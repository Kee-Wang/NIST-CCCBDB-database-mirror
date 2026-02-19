#!/usr/bin/env python3
"""
audit_gaps.py

Scan the CCCBDB mirror and classify every data issue into one of two
categories:

    Missing data   — data doesn't exist on CCCBDB or could not be fetched
    Data quality   — data exists but is wrong

Outputs:
    data/suspected_errors/data_gaps.json       — structured JSON
    data/suspected_errors/DATA_GAPS_AUDIT.md   — human-readable report

Usage:
    python scripts/audit_gaps.py               # full audit
    python scripts/audit_gaps.py --summary     # summary to stdout only
    python scripts/audit_gaps.py --verbose     # print every issue
    python scripts/audit_gaps.py --diff        # compare with previous run
"""

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from datetime import datetime

# Resolve paths relative to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from cccbdb_fetch_geometry import (
    _build_species,
    _molecule_dir,
    load_cached,
    load_species_list,
)
from db_utils import count_atoms

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Compute data paths relative to SCRIPT_DIR (not CWD, which varies)
DATA_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "data"))
SUSPECTED_DIR = os.path.join(DATA_DIR, "suspected_errors")
FETCH_STATS_DIR = os.path.join(DATA_DIR, "fetch_stats")
OUTPUT_JSON = os.path.join(SUSPECTED_DIR, "data_gaps.json")
OUTPUT_MD = os.path.join(SUSPECTED_DIR, "DATA_GAPS_AUDIT.md")
ANOMALIES_JSON = os.path.join(SUSPECTED_DIR, "anomalies.json")

# ---------------------------------------------------------------------------
# Load known anomalies
# ---------------------------------------------------------------------------

def _load_anomalies():
    """Load anomalies.json and index by species_id and category."""
    if not os.path.isfile(ANOMALIES_JSON):
        return {}
    with open(ANOMALIES_JSON, "r", encoding="utf-8") as f:
        data = json.load(f)
    index = {}
    for a in data.get("anomalies", []):
        sid = a.get("species_id", "")
        index[sid] = a
    return index

# ---------------------------------------------------------------------------
# Load fetch stats
# ---------------------------------------------------------------------------

def _load_fetch_stats():
    """Load all fetch_stats/*.json and build a species->page_type->errors index."""
    index = defaultdict(lambda: defaultdict(list))
    if not os.path.isdir(FETCH_STATS_DIR):
        return index
    for fname in sorted(os.listdir(FETCH_STATS_DIR)):
        if not fname.endswith(".json"):
            continue
        fpath = os.path.join(FETCH_STATS_DIR, fname)
        try:
            with open(fpath, "r", encoding="utf-8") as f:
                data = json.load(f)
        except (json.JSONDecodeError, IOError):
            continue
        for req in data.get("requests", []):
            if req.get("status") != "success":
                species = req.get("species", "")
                page_type = req.get("page_type", "")
                error_type = req.get("error_type", "unknown")
                if species and page_type:
                    index[species][page_type].append(error_type)
    return index

# ---------------------------------------------------------------------------
# Session failure detection
# ---------------------------------------------------------------------------

def _is_energy_session_failure(html_path):
    """Check if energy HTML file is a session failure (no real data)."""
    if not os.path.isfile(html_path):
        return False
    try:
        size = os.path.getsize(html_path)
        # Very small files are likely session failures
        if size < 500:
            return True
        with open(html_path, "r", encoding="utf-8") as f:
            content = f.read()
        # Session failure: page has no table2 (the data table)
        if 'id="table2"' not in content:
            return True
        # Page loaded but table has no data rows (bgcolor is the discriminator)
        if 'bgcolor' not in content:
            return True
    except (IOError, UnicodeDecodeError):
        return False
    return False

def _is_geom3x_session_failure(html_path):
    """Check if geom3x HTML file is a session failure (empty page with no coordinates).

    Uses BeautifulSoup to ignore HTML comments — some session failures have
    coordinate-like data inside ``<!-- ... -->`` blocks that fool raw string searches.
    """
    if not os.path.isfile(html_path):
        return False
    try:
        from bs4 import BeautifulSoup

        with open(html_path, "r", encoding="utf-8") as f:
            content = f.read()
        soup = BeautifulSoup(content, "html.parser")
        # Real geom3x pages have visible coordinate cells with class="r" or class="num"
        if soup.find(class_="r") or soup.find(class_="num"):
            return False
        return True
    except (IOError, UnicodeDecodeError):
        return False

def _is_spin_session_failure(html_path):
    """Check if spin HTML file is a session failure."""
    if not os.path.isfile(html_path):
        return False
    try:
        with open(html_path, "r", encoding="utf-8") as f:
            content = f.read()
        # The telltale sign of a session failure on spin2x.asp
        if "() is closed shell" in content:
            return True
    except (IOError, UnicodeDecodeError):
        return False
    return False

def _energy_has_molecule_name(html_path):
    """Check if energy HTML page identifies the molecule (session was established).

    A page with ``<h1>Energies for (molecule name)</h1>`` means the session
    was established but CCCBDB genuinely has no energy data for this species.
    """
    if not os.path.isfile(html_path):
        return False
    try:
        with open(html_path, "r", encoding="utf-8") as f:
            content = f.read()
        m = re.search(r"<h1>\s*Energies\s+for\s+\(([^)]+)\)", content, re.IGNORECASE)
        return bool(m)
    except (IOError, UnicodeDecodeError):
        return False

def _geom3x_has_molecule_name(html_path):
    """Check if geom3x HTML page identifies the molecule (session was established).

    A page with ``<H1> Geometry for (molecule name) ...`` means the session
    was established but CCCBDB has no calculated geometry for this species.
    """
    if not os.path.isfile(html_path):
        return False
    try:
        with open(html_path, "r", encoding="utf-8") as f:
            content = f.read()
        m = re.search(r"<[Hh]1[^>]*>\s*Geometry\s+for\s+\(([^)]+)\)", content)
        return bool(m)
    except (IOError, UnicodeDecodeError):
        return False

# ---------------------------------------------------------------------------
# Per-species gap audit
# ---------------------------------------------------------------------------

def audit_species(data, species, casno, mol_dir, anomaly_index, fetch_errors):
    """Audit one species for data gaps. Returns list of gap dicts."""
    gaps = []
    formula = data.get("formula", "")
    charge = data.get("charge", 0)
    n_atoms_formula = count_atoms(formula)
    species_id = f"{species}_{casno}"

    anomaly = anomaly_index.get(species_id)
    species_fetch_errors = fetch_errors.get(species_id, {})

    # ---------------------------------------------------------------
    # GEOMETRY gaps
    # ---------------------------------------------------------------
    geoms = data.get("geometries", {})
    has_exp = bool(geoms.get("experimental") and geoms["experimental"].get("xyz"))
    has_calc = bool(geoms.get("calculated"))
    has_any_geom = has_exp or has_calc
    is_single_atom = n_atoms_formula <= 1

    if not has_any_geom:
        if is_single_atom:
            gaps.append(_gap(species_id, "geometry", "single_atom",
                             "Single atom — no molecular geometry possible"))
        elif anomaly and anomaly.get("category") == "truncated_cartesian" and anomaly.get("status") == "upstream_unrecoverable":
            gaps.append(_gap(species_id, "geometry", "no_experimental_geometry",
                             f"Upstream: {anomaly.get('description', 'no coordinate data')}"))
        else:
            # Check if we have HTML files
            has_expgeom_html = os.path.isfile(os.path.join(mol_dir, "expgeom2x.html"))
            has_geom2x_html = os.path.isfile(os.path.join(mol_dir, "geom2x.html"))
            has_geom3x_html = os.path.isfile(os.path.join(mol_dir, "geom3x.html"))

            if not has_expgeom_html:
                # Never fetched or HTTP error
                errors = species_fetch_errors.get("expgeom2x", [])
                gaps.append(_gap(species_id, "geometry", "fetch_error",
                                 f"No expgeom2x.html archived; fetch errors: {errors or 'none logged'}"))
            elif not has_geom2x_html and not has_geom3x_html:
                # Experimental page exists but no geom summary page
                geom2x_errors = species_fetch_errors.get("geom2x", [])
                geom3x_errors = species_fetch_errors.get("geom3x", [])
                if geom2x_errors or geom3x_errors:
                    note = ""
                    if any(e == "500" for e in geom3x_errors):
                        note = " (geom3x 500 may be permanent)"
                    gaps.append(_gap(species_id, "geometry", "fetch_error",
                                     f"No calc geometry HTML; geom2x errors: {geom2x_errors}, "
                                     f"geom3x errors: {geom3x_errors}{note}",
                                     possibly_permanent=any(e == "500" for e in geom3x_errors),
                                     html_path="expgeom2x.html"))
                else:
                    # No error records either — might be geom2x showed no non-DFT methods
                    gaps.append(_gap(species_id, "geometry", "no_calc_methods",
                                     "No non-DFT calculated geometry methods available",
                                     html_path="expgeom2x.html" if has_expgeom_html else None))
            elif has_geom3x_html:
                # HTML exists but no parsed data — distinguish upstream empty vs session failure
                geom3x_path = os.path.join(mol_dir, "geom3x.html")
                if _geom3x_has_molecule_name(geom3x_path):
                    # Session established, molecule identified — CCCBDB has no geometry data
                    gaps.append(_gap(species_id, "geometry", "no_calc_geometry",
                                     "Geometry page loaded with molecule name but no coordinates",
                                     html_path="geom3x.html"))
                elif _is_geom3x_session_failure(geom3x_path):
                    gaps.append(_gap(species_id, "geometry", "session_failure",
                                     "geom3x.html exists but session not established (empty page)",
                                     html_path="geom3x.html"))
                else:
                    gaps.append(_gap(species_id, "geometry", "session_failure",
                                     "geom3x.html archived but parser found no coordinates",
                                     html_path="geom3x.html"))
            else:
                # geom2x exists but no geom3x
                geom3x_errors = species_fetch_errors.get("geom3x", [])
                if geom3x_errors:
                    gaps.append(_gap(species_id, "geometry", "fetch_error",
                                     f"geom3x fetch failed: {geom3x_errors}",
                                     possibly_permanent=any(e == "500" for e in geom3x_errors),
                                     html_path="geom2x.html"))
                else:
                    gaps.append(_gap(species_id, "geometry", "no_calc_methods",
                                     "geom2x.html exists but no non-DFT methods found",
                                     html_path="geom2x.html"))

    # ---------------------------------------------------------------
    # SPIN gaps
    # ---------------------------------------------------------------
    spin = data.get("spin")
    has_spin_html = os.path.isfile(os.path.join(mol_dir, "spin2x.html"))

    if anomaly and anomaly.get("category") in ("impossible_spin", "suspect_spin"):
        if anomaly.get("status") in ("fixed_by_parser",):
            pass  # Parser already handles it
        else:
            gaps.append(_gap(species_id, "spin", "wrong_spin",
                             anomaly.get("description", "upstream spin error"),
                             html_path="spin2x.html"))
    elif spin is None:
        # Never attempted — classify as HTTP error (fetch never succeeded)
        gaps.append(_gap(species_id, "spin", "fetch_error",
                         "Spin data never fetched (spin=None)"))
    elif spin == {}:
        # Attempted but no data
        if has_spin_html:
            if _is_spin_session_failure(os.path.join(mol_dir, "spin2x.html")):
                gaps.append(_gap(species_id, "spin", "session_failure",
                                 "spin2x.html contains session failure text '() is closed shell'",
                                 html_path="spin2x.html"))
            else:
                gaps.append(_gap(species_id, "spin", "no_spin_data",
                                 "Spin page loaded but no spin data found (genuinely empty)",
                                 html_path="spin2x.html"))
        else:
            # No HTML and empty spin — could be either
            spin_errors = species_fetch_errors.get("spin2x", [])
            if spin_errors:
                gaps.append(_gap(species_id, "spin", "fetch_error",
                                 f"spin2x fetch errors: {spin_errors}"))
            else:
                gaps.append(_gap(species_id, "spin", "no_spin_data",
                                 "No spin data and no spin HTML (assumed genuinely empty)"))

    # ---------------------------------------------------------------
    # ENERGY gaps
    # ---------------------------------------------------------------
    energy = data.get("energy")
    has_energy_html = os.path.isfile(os.path.join(mol_dir, "energy2x.html"))

    if anomaly and anomaly.get("category") == "suspect_energy":
        gaps.append(_gap(species_id, "energy", "suspect_energy",
                         anomaly.get("description", "suspect energy data"),
                         html_path="energy2x.html"))

    elif not energy:
        if has_energy_html:
            energy_html_path = os.path.join(mol_dir, "energy2x.html")
            if _energy_has_molecule_name(energy_html_path):
                # Session established, molecule identified — CCCBDB has no energy data
                gaps.append(_gap(species_id, "energy", "no_energy_data",
                                 "Energy page loaded with molecule name but no data",
                                 html_path="energy2x.html"))
            elif _is_energy_session_failure(energy_html_path):
                gaps.append(_gap(species_id, "energy", "session_failure",
                                 "energy2x.html exists but session not established (no data table)",
                                 html_path="energy2x.html"))
            else:
                gaps.append(_gap(species_id, "energy", "session_failure",
                                 "energy2x.html archived but parser found no data",
                                 html_path="energy2x.html"))
        else:
            energy_errors = species_fetch_errors.get("energy2x", [])
            if energy_errors:
                gaps.append(_gap(species_id, "energy", "fetch_error",
                                 f"No energy2x.html; fetch errors: {energy_errors}"))
            else:
                gaps.append(_gap(species_id, "energy", "fetch_error",
                                 "No energy2x.html and no fetch error records"))

    # ---------------------------------------------------------------
    # CAS typo (upstream)
    # ---------------------------------------------------------------
    if anomaly and anomaly.get("category") == "suspect_cas_typo":
        gaps.append(_gap(species_id, "metadata", "cas_typo",
                         anomaly.get("description", "suspect CAS number"),
                         html_path="expgeom2x.html"))

    return gaps

def _gap(species_id, data_type, code, description, possibly_permanent=False, html_path=None):
    """Create a gap record dict."""
    gap = {
        "species_id": species_id,
        "data_type": data_type,
        "code": code,
        "description": description,
    }
    if possibly_permanent:
        gap["possibly_permanent"] = True
    if html_path:
        gap["html_path"] = html_path
    return gap

# ---------------------------------------------------------------------------
# Full audit
# ---------------------------------------------------------------------------

def run_audit():
    """Run the full gap audit across all species. Returns structured results."""
    anomaly_index = _load_anomalies()
    fetch_errors = _load_fetch_stats()

    species_list = load_species_list()
    print(f"Auditing {len(species_list)} species...")

    all_gaps = []

    for row in species_list:
        formula = row["formula"]
        charge = int(row["charge"])
        casno = row["casno"]
        species = _build_species(formula, charge)
        mol_dir = _molecule_dir(species, casno)
        data = load_cached(species, casno)
        if not data:
            all_gaps.append(_gap(f"{species}_{casno}", "all", "fetch_error",
                                 "No molecule.json archived"))
            continue

        gaps = audit_species(data, species, casno, mol_dir, anomaly_index, fetch_errors)
        all_gaps.extend(gaps)

    # Classify into two categories
    QUALITY_CODES = {"wrong_spin", "suspect_energy"}
    missing_data = [g for g in all_gaps if g["code"] not in QUALITY_CODES]
    data_quality = [g for g in all_gaps if g["code"] in QUALITY_CODES]

    # Sub-classify by code
    by_code = defaultdict(list)
    for g in all_gaps:
        by_code[g["code"]].append(g)

    results = {
        "generated": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "total_species": len(species_list),
        "total_gaps": len(all_gaps),
        "summary": {
            "missing_data": len(missing_data),
            "data_quality": len(data_quality),
        },
        "by_code": {code: len(gaps) for code, gaps in sorted(by_code.items())},
        "gaps": {
            "missing_data": missing_data,
            "data_quality": data_quality,
        },
    }
    return results

# ---------------------------------------------------------------------------
# Output: JSON
# ---------------------------------------------------------------------------

def write_json(results, path=None):
    """Write data_gaps.json."""
    if path is None:
        path = OUTPUT_JSON
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"JSON written to: {path}")

# ---------------------------------------------------------------------------
# Output: Markdown
# ---------------------------------------------------------------------------

def write_markdown(results, path=None):
    """Write human-readable DATA_GAPS_AUDIT.md."""
    if path is None:
        path = OUTPUT_MD

    s = results["summary"]
    by_code = results["by_code"]
    gaps = results["gaps"]

    lines = [
        "# Data Gap Audit",
        "",
        f"*Generated: {results['generated']}*",
        f"*Species audited: {results['total_species']}*",
        "",
        "## Executive Summary",
        "",
        "| Category | Count | Description |",
        "|----------|------:|-------------|",
        f"| **Missing data** | {s['missing_data']} | Data doesn't exist on CCCBDB |",
        f"| **Data quality** | {s['data_quality']} | Data exists but is wrong |",
        f"| **Total** | {results['total_gaps']} | |",
        "",
        "### Breakdown by Code",
        "",
        "| Code | Count | Meaning |",
        "|------|------:|---------|",
    ]

    code_descriptions = {
        "session_failure": "HTML archived but ASP.NET session not established",
        "fetch_error": "HTML not archived — server returned HTTP 500/503/timeout",
        "single_atom": "Atoms have no molecular structure",
        "no_spin_data": "Spin page genuinely empty on CCCBDB",
        "no_experimental_geometry": "CCCBDB says 'No coordinate data'",
        "no_calc_methods": "No non-DFT calculated geometry methods available",
        "cas_typo": "Invalid CAS on species list",
        "no_energy_data": "Energy page genuinely empty on CCCBDB",
        "no_calc_geometry": "Geometry page genuinely empty on CCCBDB",
        "wrong_spin": "Parity violation or known wrong spin",
        "suspect_energy": "Identical energies for different basis sets",
    }

    for code in sorted(by_code.keys()):
        desc = code_descriptions.get(code, "")
        lines.append(f"| `{code}` | {by_code[code]} | {desc} |")
    lines.append("")

    # Missing Data details
    lines.extend([
        "## Missing Data",
        "",
        "Data genuinely does not exist on CCCBDB, or could not be fetched.",
        "",
    ])
    _write_category_tables(lines, gaps["missing_data"])

    # Data Quality details
    lines.extend([
        "## Data Quality",
        "",
        "Data exists on CCCBDB but is known to be wrong.",
        "",
    ])
    _write_category_tables(lines, gaps["data_quality"])

    # Recovery section
    lines.extend([
        "## Recovery",
        "",
        "```bash",
        "# Heal session failures (scan, clear bad HTML, re-fetch)",
        "python scripts/cccbdb_mirror.py --heal",
        "",
        "# Re-fetch all missing pages",
        "python scripts/cccbdb_mirror.py --pages all --missing-only",
        "",
        "# Re-run audit to track progress",
        "python scripts/audit_gaps.py --diff",
        "```",
        "",
    ])

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Markdown written to: {path}")

def _format_species_id(species_id, html_path=None):
    """Format a species ID, optionally as a hyperlink to its archived HTML."""
    if html_path:
        return f"[`{species_id}`](../molecules/{species_id}/{html_path})"
    return f"`{species_id}`"

def _write_compact_rows(lines, items, per_row):
    """Write items as comma-separated rows, *per_row* items each."""
    for i in range(0, len(items), per_row):
        chunk = items[i:i + per_row]
        lines.append(", ".join(chunk))

def _write_category_tables(lines, gap_list):
    """Write per-code tables for a category's gap list."""
    if not gap_list:
        lines.append("*No issues in this category.*\n")
        return

    by_code = defaultdict(list)
    for g in gap_list:
        by_code[g["code"]].append(g)

    for code in sorted(by_code.keys()):
        code_gaps = by_code[code]
        lines.append(f"### `{code}` ({len(code_gaps)} {'issue' if len(code_gaps) == 1 else 'issues'})")
        lines.append("")

        # Group by data_type for readability
        by_type = defaultdict(list)
        for g in code_gaps:
            by_type[g["data_type"]].append(g)

        for dtype in sorted(by_type.keys()):
            typed_gaps = by_type[dtype]
            lines.append(f"**{dtype}** ({len(typed_gaps)}):")
            lines.append("")

            # Check if all entries share the same description → compact format
            descriptions = set(g.get("description", "") for g in typed_gaps)
            any_html = any(g.get("html_path") for g in typed_gaps)

            if len(descriptions) == 1:
                # Shared description — print once, then list all IDs compactly
                lines.append(descriptions.pop())
                lines.append("")
                formatted = [
                    _format_species_id(g["species_id"], g.get("html_path"))
                    for g in typed_gaps
                ]
                # Hyperlinked IDs are ~60 chars each → ~4 per row;
                # plain IDs are ~15 chars → ~10 per row
                per_row = 4 if any_html else 10
                _write_compact_rows(lines, formatted, per_row)
            else:
                # Mixed descriptions — list format, hyperlink IDs when possible
                for g in typed_gaps:
                    sid_fmt = _format_species_id(g["species_id"], g.get("html_path"))
                    desc = g.get("description", "")
                    perm = " **(possibly permanent)**" if g.get("possibly_permanent") else ""
                    lines.append(f"- {sid_fmt}: {desc}{perm}")
            lines.append("")

# ---------------------------------------------------------------------------
# Diff mode
# ---------------------------------------------------------------------------

def run_diff(new_results):
    """Compare new results with previous data_gaps.json."""
    if not os.path.isfile(OUTPUT_JSON):
        print("\nNo previous data_gaps.json found — this is the first run.")
        print("Run again after retries to see progress.\n")
        return

    with open(OUTPUT_JSON, "r", encoding="utf-8") as f:
        old = json.load(f)

    old_summary = old.get("summary", {})
    new_summary = new_results.get("summary", {})
    old_by_code = old.get("by_code", {})
    new_by_code = new_results.get("by_code", {})

    print(f"\n--- Diff vs previous run ({old.get('generated', '?')}) ---")
    print(f"{'Category':<30} {'Previous':>10} {'Current':>10} {'Change':>10}")
    print("-" * 65)

    for key in ("missing_data", "data_quality"):
        o = old_summary.get(key, 0)
        n = new_summary.get(key, 0)
        delta = n - o
        sign = "+" if delta > 0 else ""
        print(f"  {key:<28} {o:>10} {n:>10} {sign}{delta:>9}")

    print()
    all_codes = sorted(set(list(old_by_code.keys()) + list(new_by_code.keys())))
    for code in all_codes:
        o = old_by_code.get(code, 0)
        n = new_by_code.get(code, 0)
        if o != n:
            delta = n - o
            sign = "+" if delta > 0 else ""
            print(f"  {code:<28} {o:>10} {n:>10} {sign}{delta:>9}")

    # Species-level diff
    old_species = set()
    for cat in old.get("gaps", {}).values():
        for g in cat:
            old_species.add((g["species_id"], g["code"]))
    new_species = set()
    for cat in new_results.get("gaps", {}).values():
        for g in cat:
            new_species.add((g["species_id"], g["code"]))

    resolved = old_species - new_species
    new_gaps = new_species - old_species

    if resolved:
        print(f"\n  Resolved since last run: {len(resolved)} issues")
        for sid, code in sorted(resolved)[:10]:
            print(f"    - {sid} ({code})")
        if len(resolved) > 10:
            print(f"    ... and {len(resolved) - 10} more")

    if new_gaps:
        print(f"\n  New issues since last run: {len(new_gaps)}")
        for sid, code in sorted(new_gaps)[:10]:
            print(f"    + {sid} ({code})")
        if len(new_gaps) > 10:
            print(f"    ... and {len(new_gaps) - 10} more")

    print()

# ---------------------------------------------------------------------------
# Summary printer
# ---------------------------------------------------------------------------

def print_summary(results):
    """Print a concise summary to stdout."""
    s = results["summary"]
    by_code = results["by_code"]

    print(f"\n=== Data Gap Audit Summary ===")
    print(f"Species audited: {results['total_species']}")
    print(f"Total issues:    {results['total_gaps']}")
    print()
    print(f"  Missing data:  {s['missing_data']}")
    print(f"  Data quality:  {s['data_quality']}")
    print()
    print(f"  {'Code':<28} {'Count':>6}")
    print(f"  {'-'*28} {'-'*6}")
    for code in sorted(by_code.keys()):
        print(f"  {code:<28} {by_code[code]:>6}")
    print()

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Audit data issues in the CCCBDB mirror archive",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/audit_gaps.py               # full audit → markdown + JSON
  python scripts/audit_gaps.py --summary     # summary to stdout only
  python scripts/audit_gaps.py --verbose     # print every issue
  python scripts/audit_gaps.py --diff        # compare with previous run
""")
    parser.add_argument("--summary", action="store_true",
                        help="Print summary to stdout only (no file output)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Print every issue")
    parser.add_argument("--diff", action="store_true",
                        help="Compare with previous data_gaps.json before overwriting")
    args = parser.parse_args()

    # Set CWD for _molecule_dir() which uses a relative CACHE_DIR path
    os.chdir(SCRIPT_DIR)

    results = run_audit()
    print_summary(results)

    if args.verbose:
        for cat_name, cat_gaps in results["gaps"].items():
            if not cat_gaps:
                continue
            print(f"\n--- {cat_name} ---")
            for g in cat_gaps:
                perm = " [possibly permanent]" if g.get("possibly_permanent") else ""
                print(f"  {g['species_id']:<30} {g['code']:<32} {g['data_type']:<10}{perm}")
                print(f"    {g['description']}")

    if args.summary:
        return

    if args.diff:
        run_diff(results)

    # Write outputs
    write_json(results)
    write_markdown(results)

if __name__ == "__main__":
    main()
