#!/usr/bin/env python3
"""
cccbdb_mirror.py

Extensible framework for mirroring CCCBDB pages as a local archive.
Raw HTML is saved as ground truth; parsed JSON provides easy access.

This file is separate from cccbdb_fetch_geometry.py (which stays as the
battle-tested, working tool). The mirror framework has a different scope:
systematic archival of ALL page types via a registry pattern.

Architecture:
    - Each CCCBDB page type is a PageType dataclass in a registry.
    - fetch_molecule() orchestrates fetching requested pages for one species.
    - reparse_molecule() re-parses saved HTML without network access.
    - New page types are added by defining a PageType + parser function.

Usage:
    python cccbdb_mirror.py --bootstrap               # init from species list HTML
    python cccbdb_mirror.py --pages spin,geometry     # fetch specific types
    python cccbdb_mirror.py --pages all               # fetch all registered
    python cccbdb_mirror.py --reparse spin            # re-parse saved HTML
    python cccbdb_mirror.py --stats                   # per-page-type coverage
    python cccbdb_mirror.py --limit 10                # cap at 10 molecules
    python cccbdb_mirror.py --missing-only            # only fetch missing data
    python cccbdb_mirror.py --heal                    # detect & heal session failures
    python cccbdb_mirror.py --verify --pages spin     # verify archived spin data
"""

import argparse
import difflib
import json
import os
import re
import time
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import date, datetime
from typing import Callable, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request

from cccbdb_fetch_geometry import (
    CACHE_DIR,
    CCCBDB_BASE,
    REQUEST_DELAY,
    _build_species,
    _molecule_dir,
    _correct_xyz_elements_from_internal_coords,
    _validate_xyz_elements,
    cache_path,
    create_session,
    fetch_with_session,
    get_best_geometry,
    is_cached,
    load_cached,
    load_species_list,
    migrate_json,
    parse_calc_geom_summary,
    parse_geometry_page,
    parse_internal_coordinates,
    parse_spin_page,
    save_cached,
    save_html,
    select_best_calc_geometry,
)
from audit_gaps import (
    _is_energy_session_failure,
    _is_geom3x_session_failure,
    _is_spin_session_failure,
    _energy_has_molecule_name,
    _geom3x_has_molecule_name,
)
from cccbdb_parse_listall import parse_species_list_html
from db_utils import _parse_formula, count_electrons

# ---------------------------------------------------------------------------
# Fetch statistics tracking
# ---------------------------------------------------------------------------

STATS_DIR = os.path.join(os.path.dirname(CACHE_DIR), "fetch_stats")

class FetchStats:
    """Track HTTP request outcomes for server health reporting."""

    def __init__(self):
        self.run_start = datetime.now().isoformat()
        self.run_end = None
        self.requests = []  # [{url, page_type, timestamp, status, error_type, duration_ms, species}]
        self._counts = defaultdict(int)  # quick summary counters

    def record(self, url, page_type, status, error_type=None,
               duration_ms=0, species=None):
        """Record a single HTTP request outcome."""
        self.requests.append({
            "url": url,
            "page_type": page_type,
            "timestamp": datetime.now().isoformat(),
            "status": status,            # "success", "http_error", "timeout", "other_error"
            "error_type": error_type,     # e.g. "500", "503", "timeout", None
            "duration_ms": round(duration_ms),
            "species": species,
        })
        self._counts[status] += 1
        if error_type:
            self._counts[f"err_{error_type}"] += 1

    def summary(self):
        """Return summary dict of all stats."""
        total = len(self.requests)
        successes = sum(1 for r in self.requests if r["status"] == "success")
        errors = total - successes

        # Per page type
        by_page = defaultdict(lambda: {"total": 0, "success": 0, "errors": {}})
        for r in self.requests:
            pt = r["page_type"]
            by_page[pt]["total"] += 1
            if r["status"] == "success":
                by_page[pt]["success"] += 1
            elif r["error_type"]:
                by_page[pt]["errors"][r["error_type"]] = \
                    by_page[pt]["errors"].get(r["error_type"], 0) + 1

        # Error type breakdown
        error_types = defaultdict(int)
        for r in self.requests:
            if r["status"] != "success" and r["error_type"]:
                error_types[r["error_type"]] += 1

        # Timing stats (successful requests only)
        durations = [r["duration_ms"] for r in self.requests
                     if r["status"] == "success" and r["duration_ms"] > 0]
        avg_ms = sum(durations) / len(durations) if durations else 0
        max_ms = max(durations) if durations else 0
        min_ms = min(durations) if durations else 0

        # Error bursts (consecutive errors)
        max_burst = 0
        current_burst = 0
        for r in self.requests:
            if r["status"] != "success":
                current_burst += 1
                max_burst = max(max_burst, current_burst)
            else:
                current_burst = 0

        return {
            "run_start": self.run_start,
            "run_end": self.run_end or datetime.now().isoformat(),
            "total_requests": total,
            "successes": successes,
            "errors": errors,
            "error_rate_pct": round(errors / total * 100, 1) if total else 0,
            "error_types": dict(error_types),
            "by_page_type": {k: dict(v) for k, v in by_page.items()},
            "timing_ms": {
                "avg": round(avg_ms),
                "min": min_ms,
                "max": max_ms,
            },
            "max_consecutive_errors": max_burst,
        }

    def save(self, path=None):
        """Save full stats to JSON file."""
        if path is None:
            os.makedirs(STATS_DIR, exist_ok=True)
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            path = os.path.join(STATS_DIR, f"fetch_stats_{ts}.json")

        self.run_end = datetime.now().isoformat()
        data = {
            "summary": self.summary(),
            "requests": self.requests,
        }
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
        print(f"  Fetch stats saved to: {path}")
        return path

    def print_summary(self):
        """Print a concise summary to stdout."""
        s = self.summary()
        print(f"\n--- Fetch Statistics ---")
        print(f"  Period: {s['run_start']} → {s['run_end']}")
        print(f"  Requests: {s['total_requests']}  "
              f"Success: {s['successes']}  "
              f"Errors: {s['errors']}  "
              f"({s['error_rate_pct']}% error rate)")

        if s["error_types"]:
            print(f"  Error breakdown:")
            for etype, count in sorted(s["error_types"].items(),
                                       key=lambda x: -x[1]):
                print(f"    {etype}: {count}")

        if s["by_page_type"]:
            print(f"  Per page type:")
            for pt, info in s["by_page_type"].items():
                errs = info.get("errors", {})
                err_str = ", ".join(f"{k}:{v}" for k, v in errs.items()) if errs else "none"
                print(f"    {pt:<14} total={info['total']:<5} "
                      f"ok={info['success']:<5} errors=[{err_str}]")

        print(f"  Response times: avg={s['timing_ms']['avg']}ms  "
              f"min={s['timing_ms']['min']}ms  max={s['timing_ms']['max']}ms")
        if s["max_consecutive_errors"] > 1:
            print(f"  Max consecutive errors: {s['max_consecutive_errors']}")

def _tracked_fetch(opener, url, page_type, stats, species=None):
    """Fetch URL with stats tracking. Replaces fetch_with_session() when stats enabled."""
    req = Request(url, headers={"User-Agent": "CCCBDB-research-fetcher/1.0"})
    t0 = time.time()
    try:
        with opener.open(req, timeout=30) as resp:
            html = resp.read().decode("utf-8", errors="replace")
        duration_ms = (time.time() - t0) * 1000
        if stats:
            stats.record(url, page_type, "success",
                         duration_ms=duration_ms, species=species)
        return html
    except HTTPError as e:
        duration_ms = (time.time() - t0) * 1000
        error_type = str(e.code)
        print(f"  ERROR fetching {url}: {e}")
        if stats:
            stats.record(url, page_type, "http_error",
                         error_type=error_type,
                         duration_ms=duration_ms, species=species)
        return None
    except (URLError, TimeoutError, OSError) as e:
        duration_ms = (time.time() - t0) * 1000
        error_str = str(e)
        if "timed out" in error_str.lower():
            error_type = "timeout"
        else:
            error_type = type(e).__name__
        print(f"  ERROR fetching {url}: {e}")
        if stats:
            stats.record(url, page_type, "timeout" if error_type == "timeout" else "other_error",
                         error_type=error_type,
                         duration_ms=duration_ms, species=species)
        return None

# ---------------------------------------------------------------------------
# PageType dataclass + registry
# ---------------------------------------------------------------------------

@dataclass
class PageType:
    """Definition of a CCCBDB page type for the mirror registry.

    Attributes:
        name:         Short identifier (e.g. "spin2x", "energy2x").
        url_template: URL template with {base}, {casno}, {charge} placeholders.
                      Session-dependent pages use "{base}/{name}.asp" (no params).
        html_filename: Filename for saved HTML (e.g. "spin2x.html").
        needs_session: True if the page relies on ASP session state.
        parser:       function(html, **ctx) -> dict | None.  May be None for stubs.
        json_key:     Top-level key in molecule.json (e.g. "spin", "energy").
        description:  Human-readable description.
        depends_on:   Page types that must be fetched first (for chained pages).
        is_session_init: True for the page that establishes the session (expgeom2x).
    """
    name: str
    url_template: str
    html_filename: str
    needs_session: bool
    parser: Optional[Callable]
    json_key: str
    description: str
    depends_on: list = field(default_factory=list)
    is_session_init: bool = False

PAGE_TYPES: dict[str, PageType] = {}

def register_page_type(pt: PageType):
    """Add a PageType to the global registry."""
    PAGE_TYPES[pt.name] = pt

# ---------------------------------------------------------------------------
# Page groups (convenience aliases for CLI --pages)
# ---------------------------------------------------------------------------

PAGE_GROUPS: dict[str, list[str]] = {
    "geometry": ["expgeom2x", "geom2x", "geom3x"],
    "spin": ["spin2x"],
    "energy": ["energy2x"],
}

def resolve_page_names(raw: str) -> list[str]:
    """Expand a comma-separated page spec into individual page names.

    Handles:
        "all"           -> all registered page types
        "spin,geometry" -> expanded via PAGE_GROUPS
        "spin2x,geom2x" -> literal page names
    """
    if raw.strip().lower() == "all":
        return list(PAGE_TYPES.keys())
    names = []
    for token in raw.split(","):
        token = token.strip().lower()
        if token in PAGE_GROUPS:
            names.extend(PAGE_GROUPS[token])
        elif token in PAGE_TYPES:
            names.append(token)
        else:
            print(f"  WARNING: unknown page type '{token}', skipping")
    # Deduplicate while preserving order
    seen = set()
    result = []
    for n in names:
        if n not in seen:
            seen.add(n)
            result.append(n)
    return result

# ---------------------------------------------------------------------------
# Topological sort for dependency ordering
# ---------------------------------------------------------------------------

def _topo_sort(page_names: list[str]) -> list[str]:
    """Sort page names respecting depends_on order.

    Also ensures the session-init page (expgeom2x) comes first when any
    session-dependent page is requested.
    """
    requested = set(page_names)

    # Add dependencies transitively
    added = True
    while added:
        added = False
        for name in list(requested):
            pt = PAGE_TYPES.get(name)
            if not pt:
                continue
            for dep in pt.depends_on:
                if dep not in requested:
                    requested.add(dep)
                    added = True

    # Ensure session-init page is present if any session-dependent page requested
    needs_session = any(PAGE_TYPES[n].needs_session for n in requested if n in PAGE_TYPES)
    session_init = [n for n in PAGE_TYPES if PAGE_TYPES[n].is_session_init]
    if needs_session and session_init:
        for si in session_init:
            requested.add(si)

    # Simple topological sort via dependency depth
    def depth(name, visited=None):
        if visited is None:
            visited = set()
        if name in visited:
            return 0
        visited.add(name)
        pt = PAGE_TYPES.get(name)
        if not pt or not pt.depends_on:
            return 0 if not pt or not pt.is_session_init else -1
        return 1 + max(depth(d, visited) for d in pt.depends_on)

    # Session-init pages get depth -1 (always first)
    def sort_key(name):
        pt = PAGE_TYPES.get(name)
        if pt and pt.is_session_init:
            return (-1, name)
        return (depth(name), name)

    return sorted(requested, key=sort_key)

# ---------------------------------------------------------------------------
# Orchestrator: fetch_molecule
# ---------------------------------------------------------------------------

def fetch_molecule(formula, name, casno, charge, pages=None, missing_only=False,
                   stats=None):
    """Fetch requested pages for one molecule, save HTML + update JSON.

    Args:
        formula:      Molecular formula (e.g. "H2O").
        name:         Human-readable name (e.g. "Water").
        casno:        CAS number string.
        charge:       Integer charge.
        pages:        List of page type names, or None for all registered.
        missing_only: If True, skip pages whose json_key already has data.
        stats:        Optional FetchStats instance for tracking request outcomes.

    Returns:
        Updated molecule data dict.
    """
    species = _build_species(formula, charge)

    # Load or initialize molecule data
    data = load_cached(species, casno)
    if data and "geometries" not in data:
        data = migrate_json(data)
    if not data:
        data = _init_molecule_data(formula, name, casno, charge)

    # Determine which pages to fetch
    if pages is None:
        page_names = list(PAGE_TYPES.keys())
    else:
        page_names = list(pages)

    ordered = _topo_sort(page_names)

    # Create session
    opener = create_session()

    # Context dict passed between pages (for chained dependencies)
    ctx = {
        "formula": formula,
        "name": name,
        "casno": casno,
        "charge": charge,
        "species": species,
        "n_electrons": count_electrons(formula, charge),
    }

    first_request = True

    mol_dir = _molecule_dir(species, casno)

    for page_name in ordered:
        pt = PAGE_TYPES.get(page_name)
        if not pt:
            continue

        # Skip if data already present and --missing-only.
        # Exception: session-init pages are ALWAYS fetched when any
        # session-dependent page will be fetched (otherwise the ASP
        # session won't know which molecule we're requesting).
        if missing_only and _has_data(data, pt, page_name, mol_dir):
            if not pt.is_session_init:
                continue
            # For session-init page, check if any later page actually needs fetching
            any_later_needed = any(
                not _has_data(data, PAGE_TYPES[pn], pn, mol_dir)
                for pn in ordered
                if pn != page_name and pn in PAGE_TYPES
                and PAGE_TYPES[pn].needs_session
            )
            if not any_later_needed:
                continue

        # Build URL
        url = _build_url(pt, ctx)

        # Rate limiting (skip delay before first request)
        if not first_request:
            time.sleep(REQUEST_DELAY)
        first_request = False

        # Fetch (use tracked version if stats enabled, otherwise plain)
        if stats:
            html = _tracked_fetch(opener, url, pt.name, stats, species=species)
        else:
            html = fetch_with_session(opener, url)
        save_html(species, casno, pt.name, html)

        if not html:
            continue

        # Parse
        if pt.parser is not None:
            parsed = pt.parser(html, **ctx)
            ctx[f"result_{page_name}"] = parsed
            _store_result(data, pt, page_name, parsed, ctx, url)

    data["fetched_date"] = str(date.today())
    save_cached(data)
    return data

# ---------------------------------------------------------------------------
# Offline re-parser
# ---------------------------------------------------------------------------

def reparse_molecule(species, casno, pages=None):
    """Re-parse saved HTML files without network access.

    For each requested page type:
        1. Check if HTML file exists in molecule folder
        2. If yes, run parser on saved HTML
        3. Update molecule.json with new parsed data

    Args:
        species: Hill formula species string (e.g. "H2O").
        casno:   CAS number string.
        pages:   List of page type names, or None for all registered.

    Returns:
        Updated molecule data dict, or None if no data archived.
    """
    data = load_cached(species, casno)
    if not data:
        return None

    if "geometries" not in data:
        data = migrate_json(data)

    if pages is None:
        page_names = list(PAGE_TYPES.keys())
    else:
        page_names = list(pages)

    ordered = _topo_sort(page_names)
    mol_dir = _molecule_dir(species, casno)

    ctx = {
        "formula": data["formula"],
        "name": data["name"],
        "casno": casno,
        "charge": data["charge"],
        "species": species,
        "n_electrons": count_electrons(data["formula"], data["charge"]),
    }

    reparsed_any = False

    for page_name in ordered:
        pt = PAGE_TYPES.get(page_name)
        if not pt or pt.parser is None:
            continue

        html_path = os.path.join(mol_dir, pt.html_filename)
        if not os.path.exists(html_path):
            continue

        with open(html_path, "r", encoding="utf-8") as f:
            html = f.read()

        parsed = pt.parser(html, **ctx)
        ctx[f"result_{page_name}"] = parsed

        url = _build_url(pt, ctx)
        _store_result(data, pt, page_name, parsed, ctx, url)
        reparsed_any = True

    # Post-processing: filter X (dummy) atoms from any calculated geometry
    # that might have been stored before the dummy atom filter was added.
    # This handles cases where geom3x.html is not archived for re-parsing.
    calc = data.get("geometries", {}).get("calculated", {})
    for method in calc:
        for basis in calc[method]:
            g = calc[method][basis]
            xyz = g.get("xyz")
            if not xyz:
                continue
            lines = [l for l in xyz.strip().split('\n')
                     if l.strip() and not l.strip().startswith('X ')]
            if len(lines) != g.get("n_atoms"):
                g["xyz"] = "\n".join(lines)
                g["n_atoms"] = len(lines)
                reparsed_any = True

    if reparsed_any:
        save_cached(data)

    return data

# ---------------------------------------------------------------------------
# Verification: re-download and compare
# ---------------------------------------------------------------------------

def _normalize_html(html):
    """Normalize HTML for comparison: collapse whitespace, strip edges."""
    if not html:
        return ""
    # Collapse runs of whitespace to a single space
    return re.sub(r'\s+', ' ', html).strip()

def verify_molecule(formula, name, casno, charge, pages=None, stats=None):
    """Re-download pages and compare against archived HTML + parsed data.

    For each requested page type:
        1. Fetch fresh HTML from CCCBDB
        2. Compare against the archived HTML on disk
        3. Re-parse the fresh HTML
        4. Compare parsed result against current molecule.json values
        5. Do NOT overwrite any files (read-only verification)

    Args:
        formula:  Molecular formula.
        name:     Human-readable name.
        casno:    CAS number string.
        charge:   Integer charge.
        pages:    List of page type names, or None for all with parsers.
        stats:    Optional FetchStats instance.

    Returns:
        dict: {page_name: {html_changed, data_changed, html_is_new, details}}
    """
    species = _build_species(formula, charge)
    data = load_cached(species, casno)
    if not data:
        return {}

    if "geometries" not in data:
        data = migrate_json(data)

    if pages is None:
        page_names = [pn for pn, pt in PAGE_TYPES.items() if pt.parser is not None]
    else:
        page_names = list(pages)

    ordered = _topo_sort(page_names)
    mol_dir = _molecule_dir(species, casno)
    opener = create_session()

    ctx = {
        "formula": formula,
        "name": name,
        "casno": casno,
        "charge": charge,
        "species": species,
        "n_electrons": count_electrons(formula, charge),
    }

    results = {}
    first_request = True

    for page_name in ordered:
        pt = PAGE_TYPES.get(page_name)
        if not pt:
            continue

        url = _build_url(pt, ctx)

        # Rate limiting
        if not first_request:
            time.sleep(REQUEST_DELAY)
        first_request = False

        # Fetch fresh HTML
        if stats:
            fresh_html = _tracked_fetch(opener, url, pt.name, stats, species=species)
        else:
            fresh_html = fetch_with_session(opener, url)

        if not fresh_html:
            results[page_name] = {
                "html_changed": False,
                "data_changed": False,
                "html_is_new": False,
                "fetch_failed": True,
                "details": "Fetch returned empty response",
            }
            continue

        # Load archived HTML
        html_path = os.path.join(mol_dir, pt.html_filename)
        old_html = None
        if os.path.isfile(html_path):
            with open(html_path, "r", encoding="utf-8") as f:
                old_html = f.read()

        # Compare HTML
        html_is_new = old_html is None
        html_changed = False
        html_details = ""
        if html_is_new:
            html_details = "No archived HTML (first fetch)"
        else:
            norm_old = _normalize_html(old_html)
            norm_new = _normalize_html(fresh_html)
            if norm_old != norm_new:
                html_changed = True
                # Generate a concise diff summary
                old_lines = old_html.splitlines(keepends=True)
                new_lines = fresh_html.splitlines(keepends=True)
                diff = list(difflib.unified_diff(
                    old_lines, new_lines,
                    fromfile="archived", tofile="fresh", n=1
                ))
                n_diff_lines = sum(1 for l in diff if l.startswith('+') or l.startswith('-'))
                html_details = f"HTML differs ({n_diff_lines} changed lines)"

        # Re-parse fresh HTML and compare with stored data
        data_changed = False
        data_details = ""
        if pt.parser is not None:
            ctx[f"result_{page_name}"] = None
            fresh_parsed = pt.parser(fresh_html, **ctx)
            ctx[f"result_{page_name}"] = fresh_parsed

            # Get stored values for comparison
            stored = _extract_stored(data, pt, page_name)
            if fresh_parsed != stored:
                data_changed = True
                data_details = _describe_data_diff(stored, fresh_parsed, page_name)

        details_parts = []
        if html_details:
            details_parts.append(html_details)
        if data_details:
            details_parts.append(data_details)

        results[page_name] = {
            "html_changed": html_changed,
            "data_changed": data_changed,
            "html_is_new": html_is_new,
            "fetch_failed": False,
            "details": "; ".join(details_parts) if details_parts else "OK",
        }

    return results

def _extract_stored(data, pt, page_name):
    """Extract the stored parsed data for a page type from molecule.json."""
    if page_name == "expgeom2x":
        return data.get("geometries", {}).get("experimental")
    elif page_name == "geom3x":
        return data.get("geometries", {}).get("calculated")
    elif page_name == "spin2x":
        return data.get("spin")
    elif page_name == "energy2x":
        return data.get("energy")
    else:
        return data.get(pt.json_key)

def _describe_data_diff(stored, fresh, page_name):
    """Produce a human-readable description of parsed data differences."""
    if stored is None and fresh is not None:
        return f"New parsed data (was empty)"
    if stored is not None and fresh is None:
        return f"Parsed data lost (fresh parse returned None)"
    if isinstance(stored, dict) and isinstance(fresh, dict):
        changed_keys = []
        all_keys = set(list(stored.keys()) + list(fresh.keys()))
        for k in sorted(all_keys):
            if stored.get(k) != fresh.get(k):
                changed_keys.append(k)
        if changed_keys:
            return f"Parsed data differs in: {', '.join(changed_keys)}"
    return f"Parsed data differs"

def _write_verify_report(all_diffs, total_verified, page_names):
    """Write VERIFY_REPORT.md summarizing verification results."""
    report_dir = os.path.normpath(os.path.join(CACHE_DIR, "..", "suspected_errors"))
    os.makedirs(report_dir, exist_ok=True)
    report_path = os.path.join(report_dir, "VERIFY_REPORT.md")

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines = [
        "# Verification Report",
        "",
        f"*Generated: {ts}*",
        f"*Species verified: {total_verified}*",
        f"*Pages checked: {', '.join(page_names)}*",
        "",
    ]

    if not all_diffs:
        lines.append("All verified pages match their archived versions. No differences found.")
    else:
        lines.extend([
            f"**{len(all_diffs)} differences found.**",
            "",
            "| Species | Page | HTML Changed | Data Changed | New HTML | Details |",
            "|---------|------|:------------:|:------------:|:--------:|---------|",
        ])
        for species_id, page_name, result in all_diffs:
            html_c = "yes" if result["html_changed"] else ""
            data_c = "yes" if result["data_changed"] else ""
            new_h = "yes" if result["html_is_new"] else ""
            details = result.get("details", "")
            lines.append(f"| `{species_id}` | `{page_name}` | {html_c} | {data_c} | {new_h} | {details} |")

    lines.append("")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Verification report written to: {report_path}")

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _init_molecule_data(formula, name, casno, charge):
    """Create a fresh molecule data dict."""
    return {
        "_schema_version": "2.1",
        "_source_database": {
            "name": "CCCBDB",
            "full_name": "Computational Chemistry Comparison and Benchmark Database",
            "version": "Release 22 (May 2022)",
            "url": "https://cccbdb.nist.gov",
            "license": "Public domain (U.S. Government work, NIST)",
        },
        "casno": casno,
        "charge": int(charge),
        "formula": formula,
        "name": name,
        "fetched_date": str(date.today()),
        "identifiers": {},
        "properties": {},
        "spin": None,
        "point_group_cccbdb": None,
        "geometries": {
            "coordinate_units": "angstrom",
            "experimental": None,
            "calculated": {},
            "best_available": None,
        },
    }

def _build_url(pt: PageType, ctx: dict) -> str:
    """Build URL for a page type, substituting context values."""
    url = pt.url_template.format(
        base=CCCBDB_BASE,
        casno=ctx["casno"],
        charge=ctx["charge"],
    )

    # geom3x needs method/basis from geom2x result
    if pt.name == "geom3x":
        geom2x_result = ctx.get("result_geom2x")
        if geom2x_result and geom2x_result.get("best"):
            best = geom2x_result["best"]
            url = (f"{CCCBDB_BASE}/geom3x.asp?"
                   f"method={best['method_id']}&basis={best['basis_id']}")
        else:
            # No calculated geometry available — return a URL that will 404
            url = f"{CCCBDB_BASE}/geom3x.asp?method=0&basis=0"

    return url

def _has_data(data: dict, pt: PageType, page_name: str, mol_dir: str = None) -> bool:
    """Check if molecule data already has content for this page type.

    Args:
        data:      Molecule data dict.
        pt:        PageType instance.
        page_name: Name of the page type.
        mol_dir:   Path to molecule directory (optional). When provided,
                   spin2x also verifies that the ground-truth HTML is archived.
    """
    if page_name == "expgeom2x":
        geoms = data.get("geometries", {})
        return (geoms.get("experimental") is not None
                and geoms["experimental"].get("xyz") is not None)
    elif page_name == "geom2x":
        # Intermediate — always re-fetch if geometry pages requested
        return False
    elif page_name == "geom3x":
        geoms = data.get("geometries", {})
        return bool(geoms.get("calculated"))
    elif page_name == "spin2x":
        spin = data.get("spin")
        # None = never attempted, {} = attempted no data, dict with keys = has data
        if spin is None:
            return False
        # Parsed data exists — also verify ground-truth HTML is archived
        if mol_dir:
            return os.path.isfile(os.path.join(mol_dir, "spin2x.html"))
        return True
    else:
        # Generic: check if json_key exists and is truthy
        return bool(data.get(pt.json_key))

def _store_result(data: dict, pt: PageType, page_name: str,
                  parsed, ctx: dict, url: str):
    """Store parsed result into the molecule data dict."""

    if page_name == "expgeom2x":
        # Store internal coordinates regardless of geometry validation
        if parsed and parsed.get("internal_coords"):
            data["geometries"]["experimental_internal_coords"] = parsed["internal_coords"]

        # Store point group
        if parsed and parsed.get("point_group"):
            data["point_group_cccbdb"] = parsed["point_group"]

        if parsed and parsed.get("xyz"):
            # Element count validation: correct or discard if elements don't match
            formula = ctx.get("formula", data.get("formula", ""))
            xyz_to_store = parsed["xyz"]
            valid = _validate_xyz_elements(xyz_to_store, formula)
            if not valid:
                # Attempt to correct element labels using internal coordinates
                ic = parsed.get("internal_coords") or data["geometries"].get(
                    "experimental_internal_coords")
                corrected = _correct_xyz_elements_from_internal_coords(
                    xyz_to_store, formula, ic)
                if corrected:
                    xyz_to_store = corrected
                    valid = True
            if valid:
                data["geometries"]["experimental"] = {
                    "xyz": xyz_to_store,
                    "point_group": parsed["point_group"],
                    "n_atoms": parsed["n_atoms"],
                    "source_url": url,
                    "geometry_type": "equilibrium",
                }
                if not data["geometries"].get("best_available"):
                    data["geometries"]["best_available"] = "experimental"
            else:
                # Bad element labels — discard experimental geometry
                data["geometries"]["experimental"] = None
                if data["geometries"].get("best_available") == "experimental":
                    # Fall back to calculated geometry if available
                    calc = data["geometries"].get("calculated", {})
                    found_calc = False
                    for method, bases in calc.items():
                        for basis in bases:
                            data["geometries"]["best_available"] = (
                                f"calculated/{method}/{basis}")
                            found_calc = True
                            break
                        if found_calc:
                            break
                    if not found_calc:
                        data["geometries"]["best_available"] = None

    elif page_name == "geom2x":
        # Intermediate result — available methods list + best selection.
        # Stored in ctx for geom3x, not directly in JSON.
        if parsed is not None:
            best = select_best_calc_geometry(parsed)
            ctx["result_geom2x"] = {"available": parsed, "best": best}
        else:
            ctx["result_geom2x"] = {"available": [], "best": None}

    elif page_name == "geom3x":
        geom2x_res = ctx.get("result_geom2x", {})
        best = geom2x_res.get("best")
        if parsed and parsed.get("xyz") and best:
            method = best["method_name"]
            basis = best["basis_name"]
            data["geometries"]["calculated"].setdefault(method, {})
            data["geometries"]["calculated"][method][basis] = {
                "xyz": parsed["xyz"],
                "point_group": (parsed["point_group"]
                                or data.get("point_group_cccbdb")),
                "n_atoms": parsed["n_atoms"],
                "method_id": best["method_id"],
                "basis_id": best["basis_id"],
                "source_url": url,
            }
            # Set as best if no experimental geometry
            if not data["geometries"].get("experimental"):
                data["geometries"]["best_available"] = (
                    f"calculated/{method}/{basis}")
            if not data.get("point_group_cccbdb") and parsed.get("point_group"):
                data["point_group_cccbdb"] = parsed["point_group"]

    elif page_name == "spin2x":
        if parsed:
            data["spin"] = parsed
        else:
            # Distinguish session failure from genuinely empty spin data.
            # Session failure: store None so --missing-only will retry.
            # Genuinely empty: store {} to avoid retrying.
            mol_dir = _molecule_dir(
                ctx.get("species", _build_species(data["formula"], data["charge"])),
                ctx.get("casno", data.get("casno", "")),
            )
            html_path = os.path.join(mol_dir, "spin2x.html")
            if os.path.isfile(html_path) and _is_spin_session_failure(html_path):
                data["spin"] = None
                os.remove(html_path)
            else:
                data["spin"] = {}

    elif page_name == "energy2x":
        if parsed:
            data["energy"] = parsed
        else:
            # Distinguish session failure from genuinely empty energy data.
            mol_dir = _molecule_dir(
                ctx.get("species", _build_species(data["formula"], data["charge"])),
                ctx.get("casno", data.get("casno", "")),
            )
            html_path = os.path.join(mol_dir, "energy2x.html")
            if os.path.isfile(html_path) and _is_energy_session_failure(html_path):
                if not _energy_has_molecule_name(html_path):
                    # True session failure — clear and allow retry
                    data["energy"] = None
                    os.remove(html_path)
                # else: upstream genuinely empty, leave energy as-is

    else:
        # Generic page types: store under json_key
        if parsed is not None:
            data[pt.json_key] = parsed

# ---------------------------------------------------------------------------
# Parser wrappers (adapt existing parsers to the **ctx signature)
# ---------------------------------------------------------------------------

def _parse_expgeom2x(html, **ctx):
    """Wrapper for parse_geometry_page (experimental geometry).

    Also parses internal coordinates and includes them in the result.
    """
    result = parse_geometry_page(html, formula=ctx.get("formula"))
    internal = parse_internal_coordinates(html)
    if internal:
        result["internal_coords"] = internal
    return result

def _parse_geom2x(html, **ctx):
    """Wrapper for parse_calc_geom_summary."""
    return parse_calc_geom_summary(html)

def _parse_geom3x(html, **ctx):
    """Wrapper for parse_geometry_page (calculated geometry)."""
    return parse_geometry_page(html, formula=ctx.get("formula"))

def _parse_spin2x(html, **ctx):
    """Wrapper for parse_spin_page with n_electrons from context."""
    return parse_spin_page(html, n_electrons=ctx.get("n_electrons"))

def _parse_energy2x(html, **ctx):
    """Parse CCCBDB energy2x.asp page to extract total energies (hartrees).

    Parses table2 ("Methods with standard basis sets"), which has the same
    structure as geom2x.asp's table2:
        - First row: 2 empty <td> + N basis <th class="nowrap"><acronym>
        - Data rows have bgcolor attribute; header repeats do not.
        - Category headers: <th rowspan="N"> (e.g. "hartree fock")
        - Method cells: <th class="nowrap"><acronym>METHOD</acronym>
        - Energy values: <td class="num"><a href="...">-76.123456</td>
        - Empty cells: <td>&nbsp;</td>

    Returns dict with:
        hf_sto3g: float or None  — HF/STO-3G energy (the reference value)
        methods: {method: {basis: energy, ...}, ...}  — all parsed energies
        basis_sets: [str, ...]  — ordered list of basis set names
    """
    if not html:
        return None

    from bs4 import BeautifulSoup

    soup = BeautifulSoup(html, "html.parser")
    table = soup.find("table", id="table2")
    if not table:
        return None

    rows = table.find_all("tr")
    if not rows:
        return None

    # Parse basis set names from first header row
    basis_names = []
    for th in rows[0].find_all("th", class_="nowrap"):
        acronym = th.find("acronym")
        if acronym:
            basis_names.append(acronym.get_text(strip=True))

    if not basis_names:
        return None

    methods = {}
    current_category = None

    for tr in rows[1:]:
        if not tr.get("bgcolor"):
            continue

        category_th = tr.find("th", attrs={"rowspan": True})
        if category_th:
            current_category = category_th.get_text(strip=True).lower()

        method_th = tr.find("th", class_="nowrap")
        if not method_th:
            continue
        acronym = method_th.find("acronym")
        if not acronym:
            continue

        method_name = acronym.get_text(strip=True)

        tds = tr.find_all("td")
        row_energies = {}
        for col, td in enumerate(tds):
            text = td.get_text(strip=True)
            if not text or text == "\xa0":
                continue
            try:
                energy = float(text)
                basis = basis_names[col] if col < len(basis_names) else f"basis_{col}"
                row_energies[basis] = energy
            except ValueError:
                continue

        if row_energies:
            # Use setdefault: first-seen value wins for each method+basis.
            # table2 can have multiple sections (total energy, atomization,
            # reaction energy, etc.) all with the same method names. The
            # first section is always total energies, which is what we want.
            existing = methods.setdefault(method_name, {})
            for basis, energy in row_energies.items():
                existing.setdefault(basis, energy)

    if not methods:
        return None

    # Extract the specific HF/STO-3G reference value
    hf_sto3g = None
    for hf_name in ("HF", "ROHF"):
        if hf_name in methods and "STO-3G" in methods[hf_name]:
            hf_sto3g = methods[hf_name]["STO-3G"]
            break

    return {
        "hf_sto3g": hf_sto3g,
        "methods": methods,
        "basis_sets": basis_names,
        "source_url": f"{CCCBDB_BASE}/energy2x.asp",
    }

def _stub_parser(html, **ctx):
    """Stub parser for page types not yet implemented. Returns None."""
    return None

# ---------------------------------------------------------------------------
# Register page types
# ---------------------------------------------------------------------------

# Priority 1: Geometry + Spin (real parsers from cccbdb_fetch_geometry.py)

register_page_type(PageType(
    name="expgeom2x",
    url_template="{base}/expgeom2x.asp?casno={casno}&charge={charge}",
    html_filename="expgeom2x.html",
    needs_session=False,
    parser=_parse_expgeom2x,
    json_key="geometries",
    description="Experimental equilibrium geometry (Cartesian coordinates)",
    is_session_init=True,
))

register_page_type(PageType(
    name="geom2x",
    url_template="{base}/geom2x.asp",
    html_filename="geom2x.html",
    needs_session=True,
    parser=_parse_geom2x,
    json_key="_geom2x_intermediate",
    description="Calculated geometry summary (available methods/basis sets)",
))

register_page_type(PageType(
    name="geom3x",
    url_template="{base}/geom3x.asp",
    html_filename="geom3x.html",
    needs_session=True,
    parser=_parse_geom3x,
    json_key="geometries",
    description="Calculated geometry (best non-DFT method, Cartesian coordinates)",
    depends_on=["geom2x"],
))

register_page_type(PageType(
    name="spin2x",
    url_template="{base}/spin2x.asp",
    html_filename="spin2x.html",
    needs_session=True,
    parser=_parse_spin2x,
    json_key="spin",
    description="Spin/multiplicity (S^2 value, closed/open shell)",
))

register_page_type(PageType(
    name="energy2x",
    url_template="{base}/energy2x.asp",
    html_filename="energy2x.html",
    needs_session=True,
    parser=_parse_energy2x,
    json_key="energy",
    description="Total energies in hartrees (HF/STO-3G reference + all methods)",
))

# Priority 2: Stubs for future page types

register_page_type(PageType(
    name="vibs2x",
    url_template="{base}/vibs2x.asp",
    html_filename="vibs2x.html",
    needs_session=True,
    parser=_stub_parser,
    json_key="vibrations",
    description="Vibrational frequencies (STUB — parser not implemented)",
))

register_page_type(PageType(
    name="dipole2x",
    url_template="{base}/dipole2x.asp",
    html_filename="dipole2x.html",
    needs_session=True,
    parser=_stub_parser,
    json_key="dipole",
    description="Dipole moments (STUB — parser not implemented)",
))

register_page_type(PageType(
    name="alldata2x",
    url_template="{base}/alldata2x.asp?casno={casno}&charge={charge}",
    html_filename="alldata2x.html",
    needs_session=False,
    parser=_stub_parser,
    json_key="overview",
    description="Overview page with InChI, names, data availability (STUB)",
))

# ---------------------------------------------------------------------------
# Coverage statistics
# ---------------------------------------------------------------------------

def print_coverage_stats(species_list):
    """Print per-page-type coverage across all cached molecules."""
    total = len(species_list)

    # Counts per page type
    html_counts = {name: 0 for name in PAGE_TYPES}
    json_counts = {name: 0 for name in PAGE_TYPES}
    cached = 0

    for row in species_list:
        species = _build_species(row["formula"], int(row["charge"]))
        casno = row["casno"]

        if not is_cached(species, casno):
            continue
        cached += 1

        mol_dir = _molecule_dir(species, casno)
        data = load_cached(species, casno)

        for name, pt in PAGE_TYPES.items():
            # Check HTML archive
            html_path = os.path.join(mol_dir, pt.html_filename)
            if os.path.exists(html_path):
                html_counts[name] += 1

            # Check parsed data in JSON
            if data and _has_data(data, pt, name):
                json_counts[name] += 1

    print(f"\nCoverage stats ({total} species, {cached} archived):")
    print(f"{'Page Type':<14} {'HTML':>6} {'Parsed':>8} {'Description'}")
    print("-" * 70)
    for name, pt in PAGE_TYPES.items():
        stub = " [STUB]" if pt.parser is _stub_parser else ""
        print(f"  {name:<12} {html_counts[name]:>6} {json_counts[name]:>8}"
              f"   {pt.description[:40]}{stub}")
    print()

# ---------------------------------------------------------------------------
# Heal session failures
# ---------------------------------------------------------------------------

# Map page types to their session failure detectors and "molecule name present"
# helpers.  A session failure means the HTML was served but the ASP.NET session
# wasn't established (so the page doesn't identify the molecule).  If the page
# *does* identify the molecule but still has no data, that's upstream-empty —
# not a session failure.
_SESSION_FAILURE_CHECKS = {
    "spin2x": {
        "is_failure": _is_spin_session_failure,
        # spin has no "molecule name present" helper — any "() is closed shell"
        # is always a session failure
        "has_molecule_name": lambda path: False,
        "json_key": "spin",
        "reset_value": None,
    },
    "energy2x": {
        "is_failure": _is_energy_session_failure,
        "has_molecule_name": _energy_has_molecule_name,
        "json_key": "energy",
        "reset_value": None,
    },
    "geom3x": {
        "is_failure": _is_geom3x_session_failure,
        "has_molecule_name": _geom3x_has_molecule_name,
        "json_key": None,  # calculated geometries are nested
        "reset_value": None,
    },
}

def heal_session_failures(species_list, max_rounds=5):
    """Detect and heal session-failed HTML files.

    Session failures occur when the CCCBDB server returns HTTP 200 but the
    ASP.NET session wasn't established, so the page body is empty/generic.
    These leave bad HTML on disk and persist {} in JSON, which --missing-only
    treats as "data exists" and skips.

    This function:
        1. Scans all species for session-failed HTML
        2. Distinguishes true session failures from upstream-empty pages
        3. Clears bad HTML + resets JSON fields
        4. Re-fetches the affected species
        5. Repeats until convergence or max_rounds

    Args:
        species_list: List of species dicts from load_species_list().
        max_rounds:   Maximum heal iterations (default 5).
    """
    print(f"=== Heal Mode ===")
    print(f"Scanning {len(species_list)} species for session failures...")

    fetch_stats = FetchStats()

    for round_num in range(1, max_rounds + 1):
        # --- Scan for session failures ---
        failures = []  # [(species, casno, formula, name, charge, page_name), ...]

        for row in species_list:
            formula = row["formula"]
            charge = int(row["charge"])
            casno = row["casno"]
            name = row["name"]
            species = _build_species(formula, charge)
            mol_dir = _molecule_dir(species, casno)

            for page_name, checks in _SESSION_FAILURE_CHECKS.items():
                html_path = os.path.join(mol_dir, f"{page_name}.html")
                if not os.path.isfile(html_path):
                    continue
                if checks["is_failure"](html_path) and not checks["has_molecule_name"](html_path):
                    failures.append((species, casno, formula, name, charge, page_name))

        if not failures:
            if round_num == 1:
                print("No session failures found. Archive is clean.")
            else:
                print(f"Round {round_num}: no remaining session failures. Done.")
            break

        # Deduplicate by species (may have multiple page failures per species)
        affected_species = {}
        for species, casno, formula, name, charge, page_name in failures:
            key = (species, casno)
            if key not in affected_species:
                affected_species[key] = {
                    "formula": formula, "name": name, "charge": charge,
                    "pages": set(),
                }
            affected_species[key]["pages"].add(page_name)

        print(f"\nRound {round_num}/{max_rounds}: found {len(failures)} session failures "
              f"across {len(affected_species)} species")

        # --- Clear bad data ---
        for (species, casno), info in affected_species.items():
            mol_dir = _molecule_dir(species, casno)
            data = load_cached(species, casno)
            if not data:
                continue

            for page_name in info["pages"]:
                checks = _SESSION_FAILURE_CHECKS[page_name]

                # Delete bad HTML
                html_path = os.path.join(mol_dir, f"{page_name}.html")
                if os.path.isfile(html_path):
                    os.remove(html_path)

                # Reset JSON field
                if page_name == "geom3x":
                    data.setdefault("geometries", {})["calculated"] = {}
                    if data.get("geometries", {}).get("best_available", "").startswith("calculated/"):
                        data["geometries"]["best_available"] = (
                            "experimental" if data["geometries"].get("experimental") else None
                        )
                else:
                    data[checks["json_key"]] = checks["reset_value"]

            save_cached(data)

        print(f"  Cleared {len(failures)} bad HTML files and reset JSON fields")

        # --- Re-fetch affected species ---
        # Determine which page groups to re-fetch per species
        recovered = 0
        still_failing = 0

        for i, ((species, casno), info) in enumerate(affected_species.items()):
            # Build the full page set needed (including session-init dependency)
            pages_needed = set()
            for page_name in info["pages"]:
                pages_needed.add(page_name)

            if i > 0:
                time.sleep(REQUEST_DELAY)

            print(f"  [{i+1}/{len(affected_species)}] Re-fetching {info['formula']} "
                  f"(CAS={casno}): {', '.join(sorted(pages_needed))}...")

            data = fetch_molecule(
                info["formula"], info["name"], casno, info["charge"],
                pages=list(pages_needed),
                stats=fetch_stats,
            )

            # Check recovery status
            for page_name in info["pages"]:
                if page_name == "spin2x":
                    if data.get("spin") and isinstance(data["spin"], dict) and data["spin"]:
                        recovered += 1
                    else:
                        still_failing += 1
                elif page_name == "energy2x":
                    if data.get("energy") and isinstance(data["energy"], dict) and data["energy"]:
                        recovered += 1
                    else:
                        still_failing += 1
                elif page_name == "geom3x":
                    if data.get("geometries", {}).get("calculated"):
                        recovered += 1
                    else:
                        still_failing += 1

        print(f"\n  Round {round_num} results: recovered={recovered}, "
              f"still failing={still_failing}")

        if still_failing == 0:
            print("All session failures recovered!")
            break
    else:
        print(f"\nReached max {max_rounds} rounds. {still_failing} failures persist.")
        print("Remaining failures may be genuine upstream issues or persistent server errors.")

    # Save fetch stats
    if fetch_stats.requests:
        fetch_stats.print_summary()
        fetch_stats.save()

# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_server_report(stats_files=None, output_path=None):
    """Generate a human-readable SERVER_STATUS.md from fetch stats JSON files.

    Args:
        stats_files: List of fetch_stats_*.json paths. If None, reads all from STATS_DIR.
        output_path: Output path for report. Defaults to data/fetch_stats/SERVER_STATUS.md.
    """
    if stats_files is None:
        if not os.path.isdir(STATS_DIR):
            print(f"No stats directory found at {STATS_DIR}")
            print("Run a fetch with stats tracking first, or use --analyze-log.")
            return
        stats_files = sorted(
            os.path.join(STATS_DIR, f) for f in os.listdir(STATS_DIR)
            if f.endswith(".json")
        )

    if not stats_files:
        print("No fetch stats files found.")
        return

    if output_path is None:
        output_path = os.path.join(STATS_DIR, "SERVER_STATUS.md")

    # Aggregate all runs
    all_requests = []
    run_summaries = []
    for sf in stats_files:
        with open(sf) as f:
            data = json.load(f)
        all_requests.extend(data.get("requests", []))
        if "summary" in data:
            run_summaries.append(data["summary"])

    total = len(all_requests)
    successes = sum(1 for r in all_requests if r["status"] == "success")
    errors = total - successes

    # Error type breakdown
    error_types = defaultdict(int)
    for r in all_requests:
        if r["status"] != "success" and r.get("error_type"):
            error_types[r["error_type"]] += 1

    # Per page type
    by_page = defaultdict(lambda: {"total": 0, "success": 0, "errors": defaultdict(int)})
    for r in all_requests:
        pt = r["page_type"]
        by_page[pt]["total"] += 1
        if r["status"] == "success":
            by_page[pt]["success"] += 1
        elif r.get("error_type"):
            by_page[pt]["errors"][r["error_type"]] += 1

    # Timing
    durations = [r["duration_ms"] for r in all_requests
                 if r["status"] == "success" and r.get("duration_ms", 0) > 0]
    avg_ms = sum(durations) / len(durations) if durations else 0

    # Species with most errors
    species_errors = defaultdict(int)
    for r in all_requests:
        if r["status"] != "success" and r.get("species"):
            species_errors[r["species"]] += 1
    worst_species = sorted(species_errors.items(), key=lambda x: -x[1])[:20]

    # Error timeline (group by hour)
    error_timeline = defaultdict(int)
    success_timeline = defaultdict(int)
    for r in all_requests:
        ts = r.get("timestamp", "")[:13]  # "2026-02-13T12" → hour bucket
        if ts:
            if r["status"] == "success":
                success_timeline[ts] += 1
            else:
                error_timeline[ts] += 1

    # Generate report
    lines = [
        "# CCCBDB Server Status Report",
        "",
        f"*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*",
        f"*Source: {len(stats_files)} fetch run(s)*",
        "",
        "## Overview",
        "",
        "This report logs all HTTP requests made during data archival from the CCCBDB",
        "website (https://cccbdb.nist.gov), a NIST-maintained ASP.NET application behind",
        "Cloudflare CDN. Transient errors (timeouts, 503s, most 500s) are expected during",
        "bulk fetching and resolve with retries — the error rate dropped from ~11% on the",
        "initial daytime pass to 0% on subsequent off-peak retries.",
        "",
        "## Aggregate Statistics",
        "",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Total HTTP requests | {total:,} |",
        f"| Successful | {successes:,} ({successes/total*100:.1f}%) |" if total else "| Successful | 0 |",
        f"| Failed | {errors:,} ({errors/total*100:.1f}%) |" if total else "| Failed | 0 |",
        f"| Avg response time (success) | {avg_ms:.0f} ms |",
        "",
    ]

    if run_summaries:
        lines.extend([
            "## Fetch Runs",
            "",
            "| Run Start | Run End | Requests | Errors | Error Rate |",
            "|-----------|---------|----------|--------|------------|",
        ])
        for s in run_summaries:
            lines.append(
                f"| {s.get('run_start', '?')[:19]} "
                f"| {s.get('run_end', '?')[:19]} "
                f"| {s.get('total_requests', 0):,} "
                f"| {s.get('errors', 0):,} "
                f"| {s.get('error_rate_pct', 0)}% |"
            )
        lines.append("")

    lines.extend([
        "## Error Types",
        "",
        "| Error Type | Count | % of Errors | Description |",
        "|------------|-------|-------------|-------------|",
    ])
    error_descriptions = {
        "500": "Internal Server Error — transient under load; permanent for some geom3x species",
        "503": "Service Unavailable — transient Cloudflare/server load, resolves on retry",
        "403": "Forbidden — Cloudflare challenge, resolves on retry",
        "timeout": "Read timeout (>30s) — transient, more common during peak hours",
        "ConnectionResetError": "TCP connection reset — transient",
        "URLError": "DNS or connection failure — transient",
    }
    for etype, count in sorted(error_types.items(), key=lambda x: -x[1]):
        pct = count / errors * 100 if errors else 0
        desc = error_descriptions.get(etype, "")
        lines.append(f"| {etype} | {count:,} | {pct:.1f}% | {desc} |")
    lines.append("")

    lines.extend([
        "## Error Rate by Page Type",
        "",
        "| Page Type | Total Requests | Successes | Errors | Error Rate | Top Error |",
        "|-----------|---------------|-----------|--------|------------|-----------|",
    ])
    for pt_name in PAGE_TYPES:
        if pt_name not in by_page:
            continue
        info = by_page[pt_name]
        t, s = info["total"], info["success"]
        e = t - s
        rate = e / t * 100 if t else 0
        top_err = max(info["errors"].items(), key=lambda x: x[1])[0] if info["errors"] else "—"
        lines.append(f"| {pt_name} | {t:,} | {s:,} | {e:,} | {rate:.1f}% | {top_err} |")
    lines.append("")

    if error_timeline:
        lines.extend([
            "## Error Timeline (by hour)",
            "",
            "| Hour | Successes | Errors | Error Rate |",
            "|------|-----------|--------|------------|",
        ])
        all_hours = sorted(set(list(error_timeline.keys()) + list(success_timeline.keys())))
        for hour in all_hours:
            s_count = success_timeline.get(hour, 0)
            e_count = error_timeline.get(hour, 0)
            total_h = s_count + e_count
            rate = e_count / total_h * 100 if total_h else 0
            lines.append(f"| {hour} | {s_count} | {e_count} | {rate:.1f}% |")
        lines.append("")

    lines.extend([
        "## Fetch Behavior Notes",
        "",
        "The error timeline above shows a clear pattern: the initial daytime fetch had",
        "~11% errors, but successive retries during off-peak hours (evenings/weekends)",
        "drove the rate to 0%. All errors except a small number of permanent `geom3x`",
        "failures were transient and caused by server load, not by bugs in the server",
        "or the client.",
        "",
        "### Transient errors resolve with retries",
        "",
        "HTTP 500, 503, and timeout errors are load-dependent. The `--missing-only` flag",
        "retries only species still missing data, and running during off-peak US Eastern",
        "hours yields near-zero error rates. Session failures (stale server-side state)",
        "are automatically detected and recovered by `--heal`.",
        "",
        "### Permanent `geom3x` 500 errors",
        "",
        "A subset of species consistently return HTTP 500 on `geom3x.asp` (calculated",
        "geometry detail pages). The server cannot compute certain method/basis",
        "combinations for these molecules. These are the only genuinely permanent",
        "failures — experimental geometry, spin, and energy data for the same species",
        "typically succeed.",
        "",
        "### Technical notes for developers",
        "",
        "- **Session requirement**: CCCBDB uses server-side ASP.NET sessions. The initial",
        "  `expgeom2x.asp?casno=X&charge=Y` request establishes the molecule context;",
        "  subsequent pages (`spin2x.asp`, `energy2x.asp`, etc.) rely on the session.",
        "  The mirror always fetches `expgeom2x.asp` first.",
        "- **Use `urllib`, not `curl`**: Cloudflare's challenge-platform token (set via",
        "  JavaScript) must be propagated. Python's `HTTPCookieProcessor` handles this;",
        "  `curl` with exported cookies does not.",
        "- **HTML tag inconsistency**: `expgeom2x.asp` uses unclosed `<TD>` tags while",
        "  `geom3x.asp` uses closed `</TD>` tags. Parsers use `(?:</TD>)?` to handle both.",
        "",
    ])

    if worst_species:
        lines.extend([
            "## Most Error-Prone Species (Top 20)",
            "",
            "| Species | Error Count |",
            "|---------|-------------|",
        ])
        for sp, count in worst_species:
            lines.append(f"| {sp} | {count} |")
        lines.append("")

    lines.extend([
        "## Re-fetch Strategy",
        "",
        "1. **Retry with `--missing-only`** — fills gaps left by transient errors.",
        "2. **Run during off-peak hours** (evenings/weekends US Eastern) — error rates",
        "   drop from ~11% to near 0%.",
        "3. **Run `--heal`** to detect and recover session failures (stale server-side",
        "   state saved as bad HTML). Scans, clears, re-fetches, loops until convergence.",
        "4. **Accept permanent `geom3x` 500s** — these species genuinely lack computed",
        "   geometry on the server. All other data types succeed for them.",
        "5. **2-second rate limit** between requests is hardcoded. Do not reduce it.",
        "",
        "---",
        f"*Report covers {len(stats_files)} fetch run(s) with {total:,} total requests.*",
        "",
    ])

    report = "\n".join(lines)
    with open(output_path, "w") as f:
        f.write(report)
    print(f"Server status report saved to: {output_path}")
    return output_path

def analyze_log_file(log_path, output_stats_path=None):
    """Parse a stdout fetch log to extract error statistics retroactively.

    Reads lines like:
        [42/2186] Si (charge=0, CAS=7440213)...
        ERROR fetching https://cccbdb.nist.gov/geom2x.asp: HTTP Error 500: Internal Server Error
        ERROR fetching https://cccbdb.nist.gov/energy2x.asp: The read operation timed out

    Returns a FetchStats-compatible dict.
    """
    species_pattern = re.compile(
        r'\[(\d+)/(\d+)\]\s+(\S+)\s+\(charge=([^,]+),\s+CAS=(\d+)\)')
    error_pattern = re.compile(
        r'ERROR fetching (https?://\S+):\s+(.*)')
    result_pattern = re.compile(
        r'->\s+best_geom=(\S+),\s+(.*)')

    requests = []
    current_species = None
    current_idx = None
    current_errors = []   # errors for current molecule
    total_species = None
    molecules_processed = 0
    molecules_with_errors = 0

    # Track which page types are being fetched (from log header)
    fetched_page_types = []
    page_types_line = re.compile(r'Fetching pages:\s+(.*)')

    def _flush_molecule(species, errors, mol_result):
        """Add success records for pages that didn't error in this molecule."""
        nonlocal molecules_with_errors
        if not species:
            return
        if errors:
            molecules_with_errors += 1

        errored_pages = {e["page_type"] for e in errors}
        # Infer successful pages: pages in fetched set that didn't error
        pages_to_check = fetched_page_types or ["expgeom2x", "energy2x", "geom2x", "spin2x", "geom3x"]
        for pt_name in pages_to_check:
            if pt_name not in errored_pages:
                requests.append({
                    "url": f"https://cccbdb.nist.gov/{pt_name}.asp",
                    "page_type": pt_name,
                    "timestamp": "",
                    "status": "success",
                    "error_type": None,
                    "duration_ms": 0,
                    "species": species,
                })

    with open(log_path) as f:
        for line in f:
            line = line.rstrip()

            # Detect which pages are being fetched
            m = page_types_line.search(line)
            if m:
                fetched_page_types = [p.strip() for p in m.group(1).split(",")]
                continue

            # New molecule
            m = species_pattern.search(line)
            if m:
                # Flush previous molecule
                _flush_molecule(current_species, current_errors, None)

                current_idx = int(m.group(1))
                total_species = int(m.group(2))
                formula = m.group(3)
                charge = m.group(4)
                casno = m.group(5)
                current_species = f"{formula}_{casno}" if charge == "0" else f"{formula}({'+' if int(charge) > 0 else ''}{charge})_{casno}"
                current_errors = []
                molecules_processed += 1
                continue

            # Error line
            m = error_pattern.search(line)
            if m:
                url = m.group(1)
                error_msg = m.group(2)

                # Determine page type from URL (most specific first)
                if "expgeom2x" in url:
                    page_type = "expgeom2x"
                elif "geom3x" in url:
                    page_type = "geom3x"
                elif "geom2x" in url:
                    page_type = "geom2x"
                elif "spin2x" in url:
                    page_type = "spin2x"
                elif "energy2x" in url:
                    page_type = "energy2x"
                else:
                    page_type = "unknown"

                # Determine error type
                if "500" in error_msg:
                    error_type = "500"
                    status = "http_error"
                elif "503" in error_msg:
                    error_type = "503"
                    status = "http_error"
                elif "403" in error_msg:
                    error_type = "403"
                    status = "http_error"
                elif "timed out" in error_msg.lower():
                    error_type = "timeout"
                    status = "timeout"
                else:
                    error_type = error_msg[:50]
                    status = "other_error"

                entry = {
                    "url": url,
                    "page_type": page_type,
                    "timestamp": "",
                    "status": status,
                    "error_type": error_type,
                    "duration_ms": 0,
                    "species": current_species,
                }
                requests.append(entry)
                current_errors.append(entry)
                continue

    # Flush last molecule
    _flush_molecule(current_species, current_errors, None)

    # Build stats
    total_requests = len(requests)
    successes = sum(1 for r in requests if r["status"] == "success")
    errors = total_requests - successes
    error_rate = errors / total_requests * 100 if total_requests else 0

    # Error type breakdown
    error_types = defaultdict(int)
    for r in requests:
        if r["status"] != "success" and r.get("error_type"):
            error_types[r["error_type"]] += 1

    # Per page type
    by_page = defaultdict(lambda: {"total": 0, "success": 0, "errors": defaultdict(int)})
    for r in requests:
        pt = r["page_type"]
        by_page[pt]["total"] += 1
        if r["status"] == "success":
            by_page[pt]["success"] += 1
        elif r.get("error_type"):
            by_page[pt]["errors"][r["error_type"]] += 1

    stats_data = {
        "summary": {
            "source": "log_analysis",
            "log_file": log_path,
            "run_start": "",
            "run_end": "",
            "molecules_processed": molecules_processed,
            "molecules_with_errors": molecules_with_errors,
            "total_species": total_species,
            "total_requests": total_requests,
            "successes": successes,
            "errors": errors,
            "error_rate_pct": round(error_rate, 1),
            "error_types": dict(error_types),
            "by_page_type": {
                k: {"total": v["total"], "success": v["success"],
                     "errors": dict(v["errors"])}
                for k, v in by_page.items()
            },
            "timing_ms": {"avg": 0, "min": 0, "max": 0},
            "max_consecutive_errors": 0,
        },
        "requests": requests,
    }

    if output_stats_path:
        os.makedirs(os.path.dirname(output_stats_path), exist_ok=True)
        with open(output_stats_path, "w") as f:
            json.dump(stats_data, f, indent=2)
        print(f"  Log analysis saved to: {output_stats_path}")

    return stats_data

# ---------------------------------------------------------------------------
# Bootstrap: create molecule directories from species list HTML
# ---------------------------------------------------------------------------

def bootstrap_from_html(html_path=None):
    """Create skeleton molecule.json files from the CCCBDB species list HTML.

    Parses data/source_pages/cccbdb_species_list.html (the archived listallx.asp
    dump) and creates a data/molecules/{species}_{casno}/ directory with a
    minimal molecule.json for each species not already on disk.

    This bridges the bootstrap gap: on a fresh clone or corrupted database,
    load_species_list() scans data/molecules/ and finds nothing. This function
    populates the directories so that subsequent --pages fetches have species
    to iterate over.

    Idempotent: skips species that already have a molecule.json.

    Args:
        html_path: Path to cccbdb_species_list.html. Defaults to
                   data/source_pages/cccbdb_species_list.html.

    Returns:
        (created, existing) tuple of counts.
    """
    species_list = parse_species_list_html(html_path)
    print(f"Parsed {len(species_list)} species from HTML")

    os.makedirs(CACHE_DIR, exist_ok=True)

    created = 0
    existing = 0

    for row in species_list:
        formula = row["formula"]
        name = row["name"]
        casno = row["casno"]
        charge = int(row["charge"])
        species = _build_species(formula, charge)

        if is_cached(species, casno):
            existing += 1
            continue

        data = _init_molecule_data(formula, name, casno, charge)
        save_cached(data)
        created += 1

    print(f"Bootstrap complete: {created} created, {existing} already exist")
    return created, existing

# ---------------------------------------------------------------------------
# CLI main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Extensible CCCBDB mirror framework",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Page types:  %(page_list)s
Groups:      %(group_list)s

Examples:
  %(prog)s --bootstrap
  %(prog)s --pages spin,geometry --limit 3
  %(prog)s --pages all --missing-only
  %(prog)s --reparse spin
  %(prog)s --stats
  %(prog)s --report
  %(prog)s --heal
  %(prog)s --verify --pages spin --limit 10
  %(prog)s --analyze-log fetch_output.log
""" % {
            "page_list": ", ".join(PAGE_TYPES.keys()),
            "group_list": ", ".join(f"{k}={v}" for k, v in PAGE_GROUPS.items()),
            "prog": "python cccbdb_mirror.py",
        })

    parser.add_argument("--pages", type=str, default=None,
                        help="Comma-separated page types or groups to fetch "
                             "(e.g. 'spin,geometry', 'all', 'energy2x')")
    parser.add_argument("--reparse", type=str, default=None,
                        help="Re-parse saved HTML for these page types "
                             "(no network access)")
    parser.add_argument("--stats", action="store_true",
                        help="Print per-page-type coverage stats")
    parser.add_argument("--limit", type=int, default=0,
                        help="Max molecules to process (0=all)")
    parser.add_argument("--missing-only", action="store_true",
                        help="Only fetch pages where data is missing")
    parser.add_argument("--list-pages", action="store_true",
                        help="List all registered page types and exit")
    parser.add_argument("--report", action="store_true",
                        help="Generate SERVER_STATUS.md from fetch stats")
    parser.add_argument("--bootstrap", action="store_true",
                        help="Initialize molecule directories from species list HTML "
                             "(data/source_pages/cccbdb_species_list.html)")
    parser.add_argument("--verify", action="store_true",
                        help="Re-download pages and compare against archived HTML "
                             "and parsed data (read-only, no files overwritten)")
    parser.add_argument("--heal", action="store_true",
                        help="Detect and heal session failures "
                             "(scan → clear → re-fetch → repeat)")
    parser.add_argument("--analyze-log", type=str, default=None,
                        metavar="LOG_FILE",
                        help="Analyze a stdout fetch log for error statistics")

    args = parser.parse_args()

    # Ensure CWD is the scripts/ directory so relative paths (../data/) resolve correctly
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Bootstrap mode
    if args.bootstrap:
        bootstrap_from_html()
        return

    # List page types
    if args.list_pages:
        print("Registered page types:")
        for name, pt in PAGE_TYPES.items():
            stub = " [STUB]" if pt.parser is _stub_parser else ""
            session = " (session)" if pt.needs_session else ""
            init = " [SESSION-INIT]" if pt.is_session_init else ""
            deps = f" depends_on={pt.depends_on}" if pt.depends_on else ""
            print(f"  {name:<14} {pt.html_filename:<18} "
                  f"{pt.description}{stub}{session}{init}{deps}")
        print(f"\nGroups:")
        for group, members in PAGE_GROUPS.items():
            print(f"  {group}: {', '.join(members)}")
        return

    # Analyze log file (offline)
    if args.analyze_log:
        print(f"Analyzing fetch log: {args.analyze_log}")
        os.makedirs(STATS_DIR, exist_ok=True)
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_path = os.path.join(STATS_DIR, f"fetch_stats_from_log_{ts}.json")
        stats_data = analyze_log_file(args.analyze_log, output_path)
        s = stats_data["summary"]
        print(f"\n--- Log Analysis ---")
        print(f"  Molecules processed: {s['molecules_processed']}")
        print(f"  Molecules with errors: {s['molecules_with_errors']}")
        print(f"  Total requests: {s.get('total_requests', 0)}  "
              f"Success: {s.get('successes', 0)}  "
              f"Errors: {s.get('errors', 0)}  "
              f"({s.get('error_rate_pct', 0)}% error rate)")
        if s["error_types"]:
            print(f"  Error breakdown:")
            for etype, count in sorted(s["error_types"].items(),
                                       key=lambda x: -x[1]):
                print(f"    {etype}: {count}")
        if s["by_page_type"]:
            print(f"  Per page type:")
            for pt, info in sorted(s["by_page_type"].items()):
                total_pt = info.get("total", 0)
                success_pt = info.get("success", 0)
                errs = info.get("errors", {})
                err_str = ", ".join(f"{k}:{v}" for k, v in errs.items()) if errs else "none"
                print(f"    {pt:<14} total={total_pt:<5} ok={success_pt:<5} errors=[{err_str}]")

        if args.report:
            generate_server_report(stats_files=[output_path])
        return

    # Report mode
    if args.report:
        generate_server_report()
        return

    # Load species list
    species_list = load_species_list()
    print(f"Loaded {len(species_list)} species from {CACHE_DIR}")

    os.makedirs(CACHE_DIR, exist_ok=True)

    # Stats mode
    if args.stats:
        print_coverage_stats(species_list)
        return

    # Heal mode
    if args.heal:
        heal_session_failures(species_list)
        return

    # Re-parse mode (offline)
    if args.reparse is not None:
        page_names = resolve_page_names(args.reparse)
        if not page_names:
            print("ERROR: no valid page types specified for --reparse")
            return

        print(f"Re-parsing: {', '.join(page_names)}")
        reparsed = 0

        for i, row in enumerate(species_list):
            if args.limit and reparsed >= args.limit:
                print(f"\nReached --limit={args.limit}, stopping.")
                break

            species = _build_species(row["formula"], int(row["charge"]))
            casno = row["casno"]

            if not is_cached(species, casno):
                continue

            result = reparse_molecule(species, casno, pages=page_names)
            if result:
                reparsed += 1
                if reparsed <= 5 or reparsed % 100 == 0:
                    print(f"  [{reparsed}] Re-parsed {row['formula']} "
                          f"(CAS={casno})")

        print(f"\nDone. Re-parsed {reparsed} molecules.")
        print_coverage_stats(species_list)
        return

    # Verify mode
    if args.verify:
        if args.pages is None:
            print("ERROR: --verify requires --pages (e.g. --verify --pages spin)")
            return

        page_names = resolve_page_names(args.pages)
        if not page_names:
            print("ERROR: no valid page types resolved from --pages")
            return

        ordered = _topo_sort(page_names)
        print(f"Verifying pages: {', '.join(ordered)}")
        print(f"  (read-only — no files will be overwritten)")

        fetch_stats = FetchStats()
        verified = 0
        all_diffs = []  # (species_id, page, result_dict)

        for i, row in enumerate(species_list):
            if args.limit and verified >= args.limit:
                print(f"\nReached --limit={args.limit}, stopping.")
                break

            formula = row["formula"]
            charge = int(row["charge"])
            casno = row["casno"]
            name = row["name"]
            species = _build_species(formula, charge)
            species_id = f"{species}_{casno}"

            if not is_cached(species, casno):
                continue

            verified += 1
            if verified <= 5 or verified % 100 == 0:
                print(f"  [{verified}] Verifying {formula} (CAS={casno})...")

            results = verify_molecule(formula, name, casno, charge,
                                      pages=ordered, stats=fetch_stats)

            for page_name, result in results.items():
                if result.get("html_changed") or result.get("data_changed") or result.get("html_is_new"):
                    all_diffs.append((species_id, page_name, result))

        # Print summary
        print(f"\n=== Verification Summary ===")
        print(f"Species verified: {verified}")
        print(f"Differences found: {len(all_diffs)}")

        if all_diffs:
            html_changes = [(s, p, r) for s, p, r in all_diffs if r["html_changed"]]
            data_changes = [(s, p, r) for s, p, r in all_diffs if r["data_changed"]]
            new_htmls = [(s, p, r) for s, p, r in all_diffs if r["html_is_new"]]

            if new_htmls:
                print(f"\n  New HTML (no archived copy): {len(new_htmls)}")
                for s, p, r in new_htmls[:10]:
                    print(f"    {s} / {p}")
                if len(new_htmls) > 10:
                    print(f"    ... and {len(new_htmls) - 10} more")

            if html_changes:
                print(f"\n  HTML content changed: {len(html_changes)}")
                for s, p, r in html_changes[:10]:
                    print(f"    {s} / {p}: {r['details']}")
                if len(html_changes) > 10:
                    print(f"    ... and {len(html_changes) - 10} more")

            if data_changes:
                print(f"\n  Parsed data changed: {len(data_changes)}")
                for s, p, r in data_changes[:10]:
                    print(f"    {s} / {p}: {r['details']}")
                if len(data_changes) > 10:
                    print(f"    ... and {len(data_changes) - 10} more")
        else:
            print("  All verified pages match their archived versions.")

        # Write verification report
        _write_verify_report(all_diffs, verified, ordered)

        # Save fetch stats
        fetch_stats.save()
        return

    # Fetch mode
    if args.pages is None:
        print("ERROR: specify --bootstrap, --pages, --reparse, --stats, --heal, --report, --verify, or --list-pages")
        parser.print_help()
        return

    page_names = resolve_page_names(args.pages)
    if not page_names:
        print("ERROR: no valid page types resolved from --pages")
        return

    ordered = _topo_sort(page_names)
    print(f"Fetching pages: {', '.join(ordered)}")
    if args.missing_only:
        print("  (--missing-only: skipping pages with existing data)")

    # Initialize fetch statistics tracker
    fetch_stats = FetchStats()

    fetched = 0
    skipped = 0

    try:
        for i, row in enumerate(species_list):
            if args.limit and fetched >= args.limit:
                print(f"\nReached --limit={args.limit}, stopping.")
                break

            casno = row["casno"]
            charge = int(row["charge"])
            formula = row["formula"]
            name = row["name"]
            species = _build_species(formula, charge)

            # In missing-only mode, skip molecules that already have all requested data
            if args.missing_only:
                data = load_cached(species, casno)
                if data:
                    mol_dir = _molecule_dir(species, casno)
                    all_present = True
                    for pn in ordered:
                        pt = PAGE_TYPES.get(pn)
                        if pt and not _has_data(data, pt, pn, mol_dir):
                            all_present = False
                            break
                    if all_present:
                        skipped += 1
                        continue

            print(f"  [{i+1}/{len(species_list)}] {formula} "
                  f"(charge={charge}, CAS={casno})...")

            data = fetch_molecule(formula, name, casno, charge,
                                  pages=ordered, missing_only=args.missing_only,
                                  stats=fetch_stats)

            best_geom = data.get("geometries", {}).get("best_available") or "none"
            spin = data.get("spin")
            spin_str = (f"mult={spin['multiplicity']}"
                        if spin and isinstance(spin, dict) and spin.get("multiplicity")
                        else "no spin")
            print(f"    -> best_geom={best_geom}, {spin_str}")

            fetched += 1

    except KeyboardInterrupt:
        print("\n\nInterrupted by user. Saving stats...")

    # Save fetch statistics
    fetch_stats.print_summary()
    stats_path = fetch_stats.save()

    print(f"\nDone. Fetched: {fetched}, Skipped: {skipped}")
    print_coverage_stats(species_list)

    # Auto-generate report
    generate_server_report(stats_files=[stats_path])

if __name__ == "__main__":
    main()
