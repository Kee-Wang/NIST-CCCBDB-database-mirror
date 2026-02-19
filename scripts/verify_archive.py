#!/usr/bin/env python3
"""
verify_archive.py

Quality control for the CCCBDB mirror archive. Validates parsed JSON
against automated sanity checks and hardcoded gold-standard values.

Tiers:
    1. Automated sanity checks across all molecules (catch gross errors)
    2. Gold-standard spot checks for ~10 molecules (catch subtle parser bugs)
    3. CAS cross-validation against PubChem (catch CAS typos / misassignments)

Usage:
    python verify_archive.py              # run all checks (Tier 1 + 2)
    python verify_archive.py --tier1      # sanity checks only
    python verify_archive.py --tier2      # gold-standard only
    python verify_archive.py --validate-cas              # CAS cross-validation (queries PubChem)
    python verify_archive.py --validate-cas --offline     # CAS validation from cache only
    python verify_archive.py --validate-cas --verbose     # show each comparison
    python verify_archive.py --verbose    # print every failure detail
"""

import argparse
import json
import os
import sys
import time
from datetime import date
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

# Resolve paths relative to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from cccbdb_fetch_geometry import (
    CACHE_DIR, _build_species, _molecule_dir, load_cached, load_species_list,
)
from cccbdb_enrich_metadata import format_cas_number
from db_utils import count_atoms, count_electrons, hill_formula

# ---------------------------------------------------------------------------
# Tier 1: Automated sanity checks
# ---------------------------------------------------------------------------

VALID_POINT_GROUPS = {
    'C1', 'Cs', 'Ci',
    'C2', 'C3', 'C4', 'C5', 'C6',
    'C2v', 'C3v', 'C4v', 'C5v', 'C6v',
    'C2h', 'C3h', 'C4h', 'C5h', 'C6h',
    'D2', 'D3', 'D4', 'D5', 'D6',
    'D2h', 'D3h', 'D4h', 'D5h', 'D6h', 'D7h', 'D8h',
    'D2d', 'D3d', 'D4d', 'D5d',
    'S4', 'S6',
    'T', 'Td', 'Th',
    'O', 'Oh',
    'I', 'Ih',
    'Kh',       # atoms
    'C\u221ev', # C∞v (linear, no inversion)
    'D\u221eh', # D∞h (linear, with inversion)
}

class SanityChecker:
    """Runs Tier 1 sanity checks across all molecules."""

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.errors = []
        self.warnings = []
        self.stats = {
            'total': 0,
            'has_geometry': 0,
            'has_spin': 0,
            'has_energy': 0,
            'checked': 0,
        }

    def error(self, species, casno, msg):
        self.errors.append(f"{species} (CAS={casno}): {msg}")
        if self.verbose:
            print(f"  FAIL: {species} (CAS={casno}): {msg}")

    def warn(self, species, casno, msg):
        self.warnings.append(f"{species} (CAS={casno}): {msg}")

    def check_molecule(self, data, species, casno):
        """Run all sanity checks on one molecule."""
        self.stats['checked'] += 1
        formula = data.get('formula', '')
        charge = data.get('charge', 0)
        n_atoms_formula = count_atoms(formula)
        n_electrons = count_electrons(formula, charge)

        self._check_geometry(data, species, casno, n_atoms_formula)
        self._check_spin(data, species, casno, n_electrons)
        self._check_energy(data, species, casno)
        self._check_point_group(data, species, casno)

    def _check_geometry(self, data, species, casno, n_atoms_formula):
        geoms = data.get('geometries', {})
        exp = geoms.get('experimental')
        calc = geoms.get('calculated', {})
        best = geoms.get('best_available')

        has_any = bool(exp and exp.get('xyz')) or bool(calc)
        if has_any:
            self.stats['has_geometry'] += 1

        # Single atoms should never have multi-line geometry
        # (CCCBDB returns a page with point group Kh but no coordinates)
        if n_atoms_formula <= 1 and exp and exp.get('xyz'):
            xyz_lines = [l for l in exp['xyz'].strip().split('\n') if l.strip()]
            if len(xyz_lines) > 1:
                self.error(species, casno,
                           f"Single atom has multi-atom geometry ({len(xyz_lines)} atoms)")

        # If experimental geometry exists, check it
        if exp and exp.get('xyz'):
            xyz = exp['xyz']
            lines = [l for l in xyz.strip().split('\n') if l.strip()]
            n_atoms_xyz = len(lines)

            # Atom count must match formula
            if n_atoms_xyz != n_atoms_formula:
                self.error(species, casno,
                           f"Geometry atom count ({n_atoms_xyz}) != "
                           f"formula atom count ({n_atoms_formula})")

            # n_atoms field must match actual XYZ lines
            if exp.get('n_atoms') != n_atoms_xyz:
                self.error(species, casno,
                           f"n_atoms field ({exp.get('n_atoms')}) != "
                           f"actual XYZ lines ({n_atoms_xyz})")

            # Coordinates should be reasonable (< 100 Angstrom)
            for line in lines:
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        coords = [abs(float(parts[i])) for i in range(1, 4)]
                        if any(c > 100 for c in coords):
                            self.error(species, casno,
                                       f"Unreasonable coordinate: {line}")
                    except ValueError:
                        self.error(species, casno,
                                   f"Non-numeric coordinate: {line}")

        # best_available pointer must be valid
        if best == 'experimental' and not (exp and exp.get('xyz')):
            self.error(species, casno,
                       "best_available='experimental' but no experimental xyz")
        if best and best.startswith('calculated/'):
            parts = best.split('/', 2)
            if len(parts) == 3:
                _, method, basis = parts
                if method not in calc or basis not in calc.get(method, {}):
                    self.error(species, casno,
                               f"best_available='{best}' but path not in calculated")

    def _check_spin(self, data, species, casno, n_electrons):
        spin = data.get('spin')
        if spin is None or spin == {}:
            return

        self.stats['has_spin'] += 1

        mult = spin.get('multiplicity')
        s2 = spin.get('S_squared')
        closed = spin.get('closed_shell')

        if mult is None:
            self.error(species, casno, "Spin data exists but multiplicity is None")
            return

        # Multiplicity must be positive integer
        if not isinstance(mult, int) or mult < 1:
            self.error(species, casno, f"Invalid multiplicity: {mult}")
            return

        # S² must be non-negative
        if s2 is not None and s2 < 0:
            self.error(species, casno, f"Negative S²: {s2}")

        # closed_shell consistency
        if closed is True and mult != 1:
            self.error(species, casno,
                       f"closed_shell=True but multiplicity={mult}")
        if closed is False and mult == 1 and s2 is not None and s2 < 0.01:
            self.warn(species, casno,
                      "closed_shell=False but multiplicity=1 and S²≈0")

        # Electron parity: odd electrons → even mult, even electrons → odd mult
        if n_electrons is not None:
            if n_electrons % 2 == 1 and mult % 2 != 0:
                self.error(species, casno,
                           f"Odd electrons ({n_electrons}) but odd multiplicity "
                           f"({mult}) — parity violation")
            if n_electrons % 2 == 0 and mult % 2 != 1:
                self.error(species, casno,
                           f"Even electrons ({n_electrons}) but even multiplicity "
                           f"({mult}) — parity violation")

    def _check_energy(self, data, species, casno):
        energy = data.get('energy')
        if not energy:
            return

        self.stats['has_energy'] += 1
        methods = energy.get('methods', {})

        if not methods:
            self.error(species, casno, "Energy dict exists but methods is empty")
            return

        # Energy values should be negative for most species.
        # Exceptions: zero-electron species (H+, He2+) have E=0,
        # and unbound anions (H2-) can have positive energies at some
        # levels of theory. Flag positive energies as warnings, not errors,
        # since they may be physically correct.
        n_electrons = count_electrons(data.get('formula', ''), data.get('charge', 0))
        for method, bases in methods.items():
            for basis, val in bases.items():
                if val > 0 and n_electrons and n_electrons > 0:
                    self.warn(species, casno,
                              f"Positive energy: {method}/{basis} = {val}")

        # HF/STO-3G should exist if HF has any data
        hf = methods.get('HF', {})
        if hf and 'STO-3G' not in hf:
            self.warn(species, casno,
                      "HF method exists but no STO-3G basis")

        # hf_sto3g field should match methods.HF.STO-3G
        hf_sto3g = energy.get('hf_sto3g')
        hf_dict_val = hf.get('STO-3G')
        if hf_sto3g is not None and hf_dict_val is not None:
            if abs(hf_sto3g - hf_dict_val) > 1e-10:
                self.error(species, casno,
                           f"hf_sto3g ({hf_sto3g}) != methods.HF.STO-3G "
                           f"({hf_dict_val})")

        # For HF: larger basis should give lower energy (variational principle)
        if hf:
            sto3g = hf.get('STO-3G')
            cc_pvqz = hf.get('cc-pVQZ')
            if sto3g is not None and cc_pvqz is not None:
                if cc_pvqz > sto3g + 1e-6:
                    self.error(species, casno,
                               f"HF variational violation: STO-3G={sto3g:.6f} < "
                               f"cc-pVQZ={cc_pvqz:.6f}")

    def _check_point_group(self, data, species, casno):
        pg = data.get('point_group_cccbdb')
        if pg and pg not in VALID_POINT_GROUPS:
            self.warn(species, casno, f"Unknown point group: {pg!r}")

    def report(self):
        s = self.stats
        print(f"\nTier 1: Automated Sanity Checks")
        print(f"  Molecules checked: {s['checked']}")
        print(f"  With geometry:     {s['has_geometry']}")
        print(f"  With spin:         {s['has_spin']}")
        print(f"  With energy:       {s['has_energy']}")
        print(f"  Errors:            {len(self.errors)}")
        print(f"  Warnings:          {len(self.warnings)}")

        if self.errors:
            print(f"\n  ERRORS ({len(self.errors)}):")
            for e in self.errors[:20]:
                print(f"    {e}")
            if len(self.errors) > 20:
                print(f"    ... and {len(self.errors) - 20} more")

        if self.warnings:
            print(f"\n  WARNINGS ({len(self.warnings)}):")
            for w in self.warnings[:10]:
                print(f"    {w}")
            if len(self.warnings) > 10:
                print(f"    ... and {len(self.warnings) - 10} more")

        return len(self.errors) == 0

# ---------------------------------------------------------------------------
# Tier 2: Gold-standard spot checks
# ---------------------------------------------------------------------------

# Hardcoded expected values manually verified against CCCBDB HTML pages.
# Each entry: (species_folder, casno, checks_dict)
# checks_dict maps dotted paths to expected values.
GOLD_STANDARD = [
    # Water — closed shell, experimental geometry, full energy
    ("H2O_7732185", "7732185", {
        "formula": "H2O",
        "charge": 0,
        "point_group_cccbdb": "C2v",
        "geometries.experimental.n_atoms": 3,
        "geometries.experimental.point_group": "C2v",
        "geometries.best_available": "experimental",
        "spin.multiplicity": 1,
        "spin.closed_shell": True,
        "spin.S_squared": 0.0,
        "energy.hf_sto3g": -74.965901,
    }),
    # Lithium hydride — linear, experimental geometry
    ("HLi_7580678", "7580678", {
        "formula": "HLi",
        "charge": 0,
        "point_group_cccbdb": "C\u221ev",
        "geometries.experimental.n_atoms": 2,
        "geometries.best_available": "experimental",
        "spin.multiplicity": 1,
        "spin.closed_shell": True,
        "energy.hf_sto3g": -7.863382,
    }),
    # Hydrogen atom — no geometry, open shell doublet
    ("H_12385136", "12385136", {
        "formula": "H",
        "charge": 0,
        "point_group_cccbdb": "Kh",
        "geometries.experimental": None,
        "geometries.best_available": None,
        "spin.multiplicity": 2,
        "spin.closed_shell": False,
        "spin.S_squared": 0.75,
        "energy.hf_sto3g": -0.466582,
    }),
    # H- anion — closed shell singlet (2 electrons)
    ("H-_12385136", "12385136", {
        "formula": "H",
        "charge": -1,
        "spin.multiplicity": 1,
        "spin.closed_shell": True,
    }),
    # O2 — open shell triplet
    ("O2_7782447", "7782447", {
        "formula": "O2",
        "charge": 0,
        "point_group_cccbdb": "D\u221eh",
        "geometries.experimental.n_atoms": 2,
        "geometries.best_available": "experimental",
        "spin.multiplicity": 3,
        "spin.closed_shell": False,
    }),
    # CO — closed shell, linear
    ("CO_630080", "630080", {
        "formula": "CO",
        "charge": 0,
        "point_group_cccbdb": "C\u221ev",
        "geometries.experimental.n_atoms": 2,
        "geometries.best_available": "experimental",
        "spin.multiplicity": 1,
        "spin.closed_shell": True,
    }),
    # CH4 — tetrahedral, closed shell
    ("CH4_74828", "74828", {
        "formula": "CH4",
        "charge": 0,
        "point_group_cccbdb": "Td",
        "geometries.experimental.n_atoms": 5,
        "geometries.best_available": "experimental",
        "spin.multiplicity": 1,
        "spin.closed_shell": True,
    }),
    # NO — open shell doublet
    ("NO_10102439", "10102439", {
        "formula": "NO",
        "charge": 0,
        "point_group_cccbdb": "C\u221ev",
        "geometries.experimental.n_atoms": 2,
        "spin.multiplicity": 2,
        "spin.closed_shell": False,
    }),
    # NH3 — pyramidal C3v
    ("H3N_7664417", "7664417", {
        "formula": "H3N",
        "charge": 0,
        "point_group_cccbdb": "C3v",
        "geometries.experimental.n_atoms": 4,
        "spin.multiplicity": 1,
        "spin.closed_shell": True,
    }),
    # N2 — linear D∞h, triple bond
    ("N2_7727379", "7727379", {
        "formula": "N2",
        "charge": 0,
        "point_group_cccbdb": "D\u221eh",
        "geometries.experimental.n_atoms": 2,
        "geometries.best_available": "experimental",
        "spin.multiplicity": 1,
        "spin.closed_shell": True,
    }),
]

def _get_nested(data, path):
    """Get a value from nested dict using dotted path."""
    parts = path.split('.')
    obj = data
    for p in parts:
        if obj is None:
            return None
        if isinstance(obj, dict):
            obj = obj.get(p)
        else:
            return None
    return obj

def run_gold_standard(verbose=False):
    """Run Tier 2 gold-standard checks. Returns (passed, failed, skipped)."""
    passed = 0
    failed = 0
    skipped = 0
    failures = []

    for folder_name, casno, checks in GOLD_STANDARD:
        mol_dir = os.path.join(CACHE_DIR, folder_name)
        json_path = os.path.join(mol_dir, "molecule.json")

        if not os.path.exists(json_path):
            skipped += 1
            if verbose:
                print(f"  SKIP: {folder_name} — not yet archived")
            continue

        with open(json_path, 'r') as f:
            data = json.load(f)

        mol_passed = True
        for path, expected in checks.items():
            actual = _get_nested(data, path)

            # Float comparison with tolerance
            if isinstance(expected, float) and isinstance(actual, (int, float)):
                if abs(actual - expected) > 1e-6:
                    mol_passed = False
                    failures.append(
                        f"{folder_name}: {path} = {actual} "
                        f"(expected {expected}, diff={abs(actual-expected):.2e})")
            elif actual != expected:
                mol_passed = False
                failures.append(
                    f"{folder_name}: {path} = {actual!r} (expected {expected!r})")

        if mol_passed:
            passed += 1
            if verbose:
                print(f"  PASS: {folder_name}")
        else:
            failed += 1

    print(f"\nTier 2: Gold-Standard Spot Checks")
    print(f"  Passed:  {passed}/{passed + failed + skipped}")
    print(f"  Failed:  {failed}")
    print(f"  Skipped: {skipped} (not yet archived)")

    if failures:
        print(f"\n  FAILURES:")
        for f in failures:
            print(f"    {f}")

    return failed == 0

# ---------------------------------------------------------------------------
# Tier 3: CAS cross-validation via PubChem
# ---------------------------------------------------------------------------

PUBCHEM_FILENAME = "pubchem.json"
PUBCHEM_RATE_LIMIT = 1.0  # seconds between requests

def _pubchem_lookup_formula(cas_formatted):
    """Query PubChem PUG REST for molecular formula and name.

    Returns {"formula": str, "name": str, "cid": int} or None if not found.
    """
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
           f"name/{cas_formatted}/property/MolecularFormula,Title/JSON")
    req = Request(url, headers={"User-Agent": "CCCBDB-mirror-validator/1.0"})
    try:
        with urlopen(req, timeout=30) as resp:
            raw = json.loads(resp.read().decode("utf-8"))
    except HTTPError as e:
        if e.code == 404:
            return None
        return None
    except (URLError, TimeoutError, OSError):
        return None

    try:
        props = raw["PropertyTable"]["Properties"][0]
    except (KeyError, IndexError):
        return None

    return {
        "formula": props.get("MolecularFormula"),
        "name": props.get("Title"),
        "cid": props.get("CID"),
    }

def _load_pubchem_cache(molecules_dir):
    """Scan molecule directories for existing pubchem.json files.

    Returns dict: cas_formatted -> pubchem result dict (or None sentinel for not-found).
    """
    cache = {}
    if not os.path.isdir(molecules_dir):
        return cache
    for entry in os.listdir(molecules_dir):
        ppath = os.path.join(molecules_dir, entry, PUBCHEM_FILENAME)
        if not os.path.isfile(ppath):
            continue
        try:
            with open(ppath, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except (json.JSONDecodeError, IOError):
            continue
        cas_fmt = data.get("cas_formatted")
        if not cas_fmt:
            continue
        # Already seen this CAS from another charge-state dir — skip
        if cas_fmt in cache:
            continue
        if data.get("pubchem_cid") is None:
            cache[cas_fmt] = None  # not-found sentinel
        else:
            cache[cas_fmt] = {
                "formula": data.get("pubchem_formula"),
                "name": data.get("pubchem_name"),
                "cid": data.get("pubchem_cid"),
            }
    return cache

def _write_pubchem_file(mol_dir, cas_formatted, pubchem_result, our_formula_hill):
    """Write pubchem.json to a molecule directory."""
    if pubchem_result is not None:
        pubchem_formula = pubchem_result.get("formula", "")
        pubchem_hill = hill_formula(pubchem_formula) if pubchem_formula else ""
        record = {
            "cas_formatted": cas_formatted,
            "pubchem_cid": pubchem_result.get("cid"),
            "pubchem_formula": pubchem_formula,
            "pubchem_name": pubchem_result.get("name"),
            "our_formula_hill": our_formula_hill,
            "pubchem_formula_hill": pubchem_hill,
            "match": (our_formula_hill == pubchem_hill),
            "queried": str(date.today()),
        }
    else:
        record = {
            "cas_formatted": cas_formatted,
            "pubchem_cid": None,
            "pubchem_formula": None,
            "pubchem_name": None,
            "our_formula_hill": our_formula_hill,
            "match": None,
            "queried": str(date.today()),
        }
    os.makedirs(mol_dir, exist_ok=True)
    with open(os.path.join(mol_dir, PUBCHEM_FILENAME), 'w', encoding='utf-8') as f:
        json.dump(record, f, indent=2, ensure_ascii=False)
        f.write('\n')

def validate_cas_numbers(molecules_dir, offline=False, verbose=False):
    """Cross-validate all CAS numbers against PubChem.

    For each species:
    1. Format CAS with dashes
    2. Look up in per-species pubchem.json or query PubChem (unless offline)
    3. Hill-normalize PubChem formula and compare to our formula
    4. Write pubchem.json to all molecule dirs sharing that CAS

    Returns (errors, warnings, passed).
    """
    # Load existing per-species pubchem.json files as cache
    cache = _load_pubchem_cache(molecules_dir)
    species_list = load_species_list()

    errors = []
    warnings = []
    matches = 0
    not_found_count = 0
    skipped = 0
    new_lookups = 0

    # Deduplicate by CAS (multiple charges share the same CAS)
    # Also collect molecule dirs per CAS for writing pubchem.json
    cas_to_species = {}
    for row in species_list:
        casno = row["casno"]
        formula = row["formula"]
        cas_fmt = format_cas_number(casno)
        charge = int(row.get("charge", 0))
        species = _build_species(formula, charge)
        if cas_fmt not in cas_to_species:
            cas_to_species[cas_fmt] = []
        cas_to_species[cas_fmt].append({
            "formula": formula,
            "charge": charge,
            "casno": casno,
            "species": species,
            "mol_dir": _molecule_dir(species, casno),
        })

    # Track pending writes for intermediate flush
    pending_writes = []

    for cas_fmt, species_group in sorted(cas_to_species.items()):
        # Use the neutral species formula for comparison (or first available)
        ref = None
        for s in species_group:
            if s["charge"] == 0:
                ref = s
                break
        if ref is None:
            ref = species_group[0]

        our_formula_hill = hill_formula(ref["formula"])

        # Check cache first (covers both found and not-found)
        if cas_fmt in cache:
            pc = cache[cas_fmt]
        elif offline:
            skipped += 1
            if verbose:
                print(f"  SKIP (offline): {ref['species']} CAS={cas_fmt}")
            continue
        else:
            # Query PubChem
            pc = _pubchem_lookup_formula(cas_fmt)
            cache[cas_fmt] = pc
            new_lookups += 1

            # Queue writes to all molecule dirs sharing this CAS
            for s in species_group:
                s_hill = hill_formula(s["formula"])
                pending_writes.append((s["mol_dir"], cas_fmt, pc, s_hill))

            # Intermediate flush every 50 lookups
            if new_lookups % 50 == 0:
                print(f"  PubChem lookups: {new_lookups}...")
                for args in pending_writes:
                    _write_pubchem_file(*args)
                pending_writes.clear()
            time.sleep(PUBCHEM_RATE_LIMIT)

        if pc is None:
            not_found_count += 1
            if verbose:
                print(f"  NOT FOUND: {ref['species']} CAS={cas_fmt}")
            continue

        # Compare formulas (Hill-normalized)
        pubchem_formula = pc.get("formula", "")
        pubchem_hill = hill_formula(pubchem_formula) if pubchem_formula else ""

        if our_formula_hill == pubchem_hill:
            matches += 1
            if verbose:
                print(f"  MATCH: {ref['species']} CAS={cas_fmt} "
                      f"formula={our_formula_hill}")
        else:
            # Check if this is a charged species vs neutral parent mismatch
            is_charged = any(s["charge"] != 0 for s in species_group)
            if is_charged and len(species_group) > 1:
                # Multiple charge states share the CAS — check if any neutral matches
                neutral_matches = any(
                    hill_formula(s["formula"]) == pubchem_hill
                    for s in species_group if s["charge"] == 0
                )
                if neutral_matches:
                    matches += 1
                    if verbose:
                        print(f"  MATCH: {ref['species']} CAS={cas_fmt} "
                              f"(neutral parent matches)")
                    continue

            # Check if difference is only due to charge (PubChem returns neutral)
            if ref["charge"] != 0:
                warnings.append(
                    f"{ref['species']} (CAS={cas_fmt}): "
                    f"our={our_formula_hill} vs PubChem={pubchem_hill} "
                    f"(charged species, PubChem may return neutral parent)")
                if verbose:
                    print(f"  WARN (charged): {ref['species']} CAS={cas_fmt} "
                          f"our={our_formula_hill} vs PubChem={pubchem_hill}")
            else:
                errors.append(
                    f"{ref['species']} (CAS={cas_fmt}): "
                    f"our={our_formula_hill} vs PubChem={pubchem_hill} "
                    f"(name: {pc.get('name', '?')})")
                if verbose:
                    print(f"  MISMATCH: {ref['species']} CAS={cas_fmt} "
                          f"our={our_formula_hill} vs PubChem={pubchem_hill} "
                          f"(name: {pc.get('name', '?')})")

    # Flush remaining pending writes
    for args in pending_writes:
        _write_pubchem_file(*args)

    # Report
    total_cas = len(cas_to_species)
    cached = total_cas - new_lookups - skipped
    print(f"\nTier 3: CAS Validation (PubChem cross-check)")
    print(f"  Unique CAS numbers:   {total_cas}")
    print(f"  PubChem matches:      {matches} (formula confirmed)")
    print(f"  PubChem not found:    {not_found_count} (CAS not in PubChem)")
    print(f"  Formula mismatches:   {len(errors)}")
    print(f"  Warnings:             {len(warnings)} (charged species / salt differences)")
    print(f"  New PubChem lookups:  {new_lookups} (cached: {cached})")
    if skipped:
        print(f"  Skipped (offline):    {skipped}")

    if errors:
        print(f"\n  MISMATCHES ({len(errors)}):")
        for e in errors[:20]:
            print(f"    {e}")
        if len(errors) > 20:
            print(f"    ... and {len(errors) - 20} more")

    if warnings and verbose:
        print(f"\n  WARNINGS ({len(warnings)}):")
        for w in warnings[:20]:
            print(f"    {w}")
        if len(warnings) > 20:
            print(f"    ... and {len(warnings) - 20} more")

    return errors, warnings, matches

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Verify CCCBDB archive integrity")
    parser.add_argument("--tier1", action="store_true",
                        help="Run only Tier 1 sanity checks")
    parser.add_argument("--tier2", action="store_true",
                        help="Run only Tier 2 gold-standard checks")
    parser.add_argument("--validate-cas", action="store_true",
                        help="Run CAS cross-validation against PubChem")
    parser.add_argument("--offline", action="store_true",
                        help="Use cached data only (no network requests)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Print detailed output for every check")
    args = parser.parse_args()

    # Ensure CWD is the scripts/ directory so relative paths (../data/) resolve correctly
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # If specific tiers or --validate-cas requested, run only those.
    # Otherwise run Tier 1 + 2 (default, excludes CAS validation since it's slow).
    run_tier1 = args.tier1 or (not args.tier2 and not args.validate_cas)
    run_tier2 = args.tier2 or (not args.tier1 and not args.validate_cas)

    species_list = load_species_list()
    print(f"Loaded {len(species_list)} species")

    all_passed = True

    # Tier 1
    if run_tier1:
        checker = SanityChecker(verbose=args.verbose)
        for row in species_list:
            species = _build_species(row["formula"], int(row["charge"]))
            casno = row["casno"]
            data = load_cached(species, casno)
            if not data:
                continue
            checker.stats['total'] += 1
            checker.check_molecule(data, species, casno)
        if not checker.report():
            all_passed = False

    # Tier 2
    if run_tier2:
        if not run_gold_standard(verbose=args.verbose):
            all_passed = False

    # CAS validation
    if args.validate_cas:
        cas_errors, cas_warnings, cas_matches = validate_cas_numbers(
            CACHE_DIR, offline=args.offline, verbose=args.verbose
        )
        if cas_errors:
            all_passed = False

    # Summary
    print()
    if all_passed:
        print("RESULT: ALL CHECKS PASSED")
    else:
        print("RESULT: SOME CHECKS FAILED")
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())
