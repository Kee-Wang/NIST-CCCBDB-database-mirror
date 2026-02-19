#!/usr/bin/env python3
"""
cccbdb_enrich_metadata.py

Purpose:
    Post-processing enrichment pass over per-molecule JSON files in
    data/molecules/. Adds schema version, provenance metadata,
    computed properties (n_electrons, molecular_weight, n_atoms),
    external identifiers (InChI, SMILES, PubChem CID via PubChem API),
    and data quality annotations.

    All changes are additive — no existing fields are modified or removed.
    Downstream code (e.g. get_best_geometry()) reads only the original
    fields and is unaffected.

Input:  data/molecules/*.json (existing per-molecule cache)
Output: Same files, enriched in-place with new fields.

Usage:
    python cccbdb_enrich_metadata.py              # enrich all JSONs
    python cccbdb_enrich_metadata.py --limit 10   # enrich at most 10
    python cccbdb_enrich_metadata.py --stats       # show enrichment stats
    python cccbdb_enrich_metadata.py --skip-pubchem # skip PubChem lookups (offline)
"""

import argparse
import json
import os
import time
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError

import db_utils

CACHE_DIR = "../data/molecules"
SCHEMA_VERSION = "2.0"
PUBCHEM_DELAY = 1.0  # seconds between PubChem requests

SOURCE_DATABASE = {
    "name": "CCCBDB",
    "full_name": "Computational Chemistry Comparison and Benchmark Database",
    "version": "Release 22 (May 2022)",
    "url": "https://cccbdb.nist.gov",
    "license": "Public domain (U.S. Government work, NIST)",
}

# Standard atomic weights (IUPAC 2021) for elements present in CCCBDB (Z <= 54 + extras).
# Values are conventional atomic weights where available.
ATOMIC_MASSES = {
    'H': 1.008, 'D': 2.014, 'T': 3.016, 'He': 4.003,
    'Li': 6.941, 'Be': 9.012, 'B': 10.81, 'C': 12.011, 'N': 14.007, 'O': 15.999,
    'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974, 'S': 32.06,
    'Cl': 35.45, 'Ar': 39.948,
    'K': 39.098, 'Ca': 40.078, 'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996,
    'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se': 78.971, 'Br': 79.904, 'Kr': 83.798,
    'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224,
    'Nb': 92.906, 'Mo': 95.95, 'Ru': 101.07, 'Rh': 102.906, 'Pd': 106.42,
    'Ag': 107.868, 'Cd': 112.414, 'In': 114.818, 'Sn': 118.710,
    'Sb': 121.760, 'Te': 127.60, 'I': 126.904, 'Xe': 131.293,
}

def format_cas_number(casno):
    """Format raw CAS digit string to standard NN...N-NN-N format.

    CAS Registry Numbers have the form PREFIX-MIDDLE-CHECK where:
    - CHECK is a single check digit (last digit)
    - MIDDLE is a 2-digit group (second and third from last)
    - PREFIX is the remaining leading digits

    Example: "7732185" -> "7732-18-5" (Water)
    """
    s = str(casno)
    if len(s) < 4:
        return s  # too short to format
    check = s[-1]
    middle = s[-3:-1]
    prefix = s[:-3]
    return f"{prefix}-{middle}-{check}"

def compute_molecular_weight(formula):
    """Compute molecular weight from formula using ATOMIC_MASSES.

    Returns float rounded to 3 decimal places, or None if an element is unknown.
    """
    atoms = db_utils._parse_formula(formula)
    total = 0.0
    for el, count in atoms.items():
        if el not in ATOMIC_MASSES:
            return None
        total += count * ATOMIC_MASSES[el]
    return round(total, 3)

def lookup_pubchem(cas_number):
    """Query PubChem REST API for identifiers by CAS number.

    Returns dict with keys: pubchem_cid, inchi, inchi_key, smiles.
    Returns None if not found or on error.
    """
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
           f"name/{cas_number}/JSON")
    req = Request(url, headers={"User-Agent": "CCCBDB-research-enricher/1.0"})
    try:
        with urlopen(req, timeout=30) as resp:
            raw = json.loads(resp.read().decode("utf-8"))
    except HTTPError as e:
        if e.code == 404:
            return None  # not found
        return None
    except (URLError, TimeoutError, OSError):
        return None

    try:
        compound = raw["PC_Compounds"][0]
    except (KeyError, IndexError):
        return None

    result = {"pubchem_cid": compound.get("id", {}).get("id", {}).get("cid")}

    # Extract properties from props array
    for prop in compound.get("props", []):
        urn = prop.get("urn", {})
        label = urn.get("label", "")
        name = urn.get("name", "")
        val = prop.get("value", {})

        if label == "InChI" and name == "Standard":
            result["inchi"] = val.get("sval")
        elif label == "InChIKey" and name == "Standard":
            result["inchi_key"] = val.get("sval")
        elif label == "SMILES" and name == "Canonical":
            result["smiles"] = val.get("sval")

    return result

def enrich_json(data, pubchem_cache=None, skip_pubchem=False):
    """Add all enrichment fields to a data dict. Idempotent.

    If _schema_version is already "2.0", only fills missing fields (does not
    overwrite existing values). This makes the function safe to run multiple times.

    Args:
        data: dict loaded from a molecule JSON file
        pubchem_cache: optional dict mapping CAS number (formatted) -> pubchem result.
            Used to share lookups across species with the same CAS.
        skip_pubchem: if True, skip PubChem API lookups entirely.

    Returns:
        (data, pubchem_was_fetched) tuple. pubchem_was_fetched is True if a new
        PubChem API call was made (for rate limiting by caller).
    """
    pubchem_fetched = False

    # Schema version
    data.setdefault("_schema_version", SCHEMA_VERSION)

    # Source database provenance
    data.setdefault("_source_database", SOURCE_DATABASE)

    # Properties (computed locally, no network)
    formula = data.get("formula", "")
    charge = int(data.get("charge", 0))

    if "properties" not in data:
        data["properties"] = {}
    props = data["properties"]
    if "n_electrons" not in props:
        props["n_electrons"] = db_utils.count_electrons(formula, charge)
    if "molecular_weight" not in props:
        props["molecular_weight"] = compute_molecular_weight(formula)
    if "n_atoms" not in props:
        props["n_atoms"] = db_utils.count_atoms(formula)

    # Identifiers
    casno = data.get("casno", "")
    if "identifiers" not in data:
        data["identifiers"] = {}
    idents = data["identifiers"]

    # Always ensure formatted CAS is present
    if "cas_number" not in idents and casno:
        idents["cas_number"] = format_cas_number(casno)

    # PubChem lookup (skip if disabled or already populated)
    if not skip_pubchem and idents.get("pubchem_cid") is None:
        cas_formatted = idents.get("cas_number", format_cas_number(casno))

        if pubchem_cache is not None and cas_formatted in pubchem_cache:
            pc = pubchem_cache[cas_formatted]
        else:
            pc = lookup_pubchem(cas_formatted)
            pubchem_fetched = True
            if pubchem_cache is not None:
                pubchem_cache[cas_formatted] = pc

        if pc:
            idents.setdefault("pubchem_cid", pc.get("pubchem_cid"))
            idents.setdefault("inchi", pc.get("inchi"))
            idents.setdefault("inchi_key", pc.get("inchi_key"))
            idents.setdefault("smiles", pc.get("smiles"))

    # Spin source tag
    if data.get("spin") and isinstance(data["spin"], dict) and data["spin"]:
        data["spin"].setdefault("source", "cccbdb_spin2x")

    # Geometry annotations
    if "geometries" in data:
        geoms = data["geometries"]
        geoms.setdefault("coordinate_units", "angstrom")

        if geoms.get("experimental") and isinstance(geoms["experimental"], dict):
            geoms["experimental"].setdefault("geometry_type", "equilibrium")

    return data, pubchem_fetched

def _iter_molecule_jsons(cache_dir):
    """Iterate over all molecule.json files in the cache directory.

    Expects per-molecule folder layout: {species}_{casno}/molecule.json.

    Yields (fpath, data) tuples.
    """
    if not os.path.exists(cache_dir):
        return
    for entry in sorted(os.listdir(cache_dir)):
        entry_path = os.path.join(cache_dir, entry)
        if not os.path.isdir(entry_path):
            continue
        fpath = os.path.join(entry_path, "molecule.json")
        if not os.path.exists(fpath):
            continue
        try:
            with open(fpath, 'r', encoding='utf-8') as f:
                data = json.load(f)
            yield fpath, data
        except (json.JSONDecodeError, IOError):
            continue

def print_stats():
    """Print enrichment coverage statistics."""
    if not os.path.exists(CACHE_DIR):
        print("No cache directory found.")
        return

    total = 0
    has_schema = 0
    has_properties = 0
    has_identifiers = 0
    has_pubchem = 0
    has_spin_source = 0
    has_coord_units = 0

    for fpath, data in _iter_molecule_jsons(CACHE_DIR):
        total += 1

        if data.get("_schema_version") == SCHEMA_VERSION:
            has_schema += 1
        if data.get("properties"):
            has_properties += 1
        idents = data.get("identifiers")
        if idents and idents.get("cas_number"):
            has_identifiers += 1
        if idents and idents.get("pubchem_cid"):
            has_pubchem += 1
        spin = data.get("spin")
        if spin and isinstance(spin, dict) and spin.get("source"):
            has_spin_source += 1
        geoms = data.get("geometries")
        if geoms and geoms.get("coordinate_units"):
            has_coord_units += 1

    print(f"Enrichment stats ({total} JSON files):")
    print(f"  Schema v{SCHEMA_VERSION}:     {has_schema}/{total}")
    print(f"  Properties:        {has_properties}/{total}")
    print(f"  Identifiers (CAS): {has_identifiers}/{total}")
    print(f"  PubChem CID:       {has_pubchem}/{total}")
    print(f"  Spin source tag:   {has_spin_source}/{total}")
    print(f"  Coordinate units:  {has_coord_units}/{total}")

def main():
    parser = argparse.ArgumentParser(
        description="Enrich CCCBDB molecule JSONs with metadata")
    parser.add_argument("--limit", type=int, default=0,
                        help="Max files to enrich (0=all)")
    parser.add_argument("--stats", action="store_true",
                        help="Print enrichment stats only")
    parser.add_argument("--skip-pubchem", action="store_true",
                        help="Skip PubChem API lookups (offline mode)")
    args = parser.parse_args()

    # Ensure CWD is the scripts/ directory so relative paths (../data/) resolve correctly
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    if args.stats:
        print_stats()
        return

    if not os.path.exists(CACHE_DIR):
        print(f"Cache directory not found: {CACHE_DIR}")
        return

    # Collect all molecule JSON paths
    json_entries = list(_iter_molecule_jsons(CACHE_DIR))
    print(f"Found {len(json_entries)} molecule JSON files in {CACHE_DIR}")

    enriched = 0
    pubchem_lookups = 0
    pubchem_cache = {}  # CAS -> pubchem result (shared across charges)

    for i, (fpath, data) in enumerate(json_entries):
        if args.limit and enriched >= args.limit:
            print(f"\nReached --limit={args.limit}, stopping.")
            break

        data, pubchem_fetched = enrich_json(
            data, pubchem_cache=pubchem_cache, skip_pubchem=args.skip_pubchem)

        with open(fpath, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)

        enriched += 1

        if pubchem_fetched:
            pubchem_lookups += 1
            if pubchem_lookups % 50 == 0:
                print(f"  [{i+1}/{len(json_entries)}] PubChem lookups: {pubchem_lookups}")
            time.sleep(PUBCHEM_DELAY)

    print(f"\nDone. Enriched: {enriched}, PubChem lookups: {pubchem_lookups}")
    print()
    print_stats()

if __name__ == "__main__":
    main()
