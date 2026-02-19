#!/usr/bin/env python3
"""
export_csv.py

Generates cccbdb-selected-energy.csv directly from per-molecule JSON files.
One row per species, 20 columns: species info, geometry, cheapest energy,
and best correlated energy.

Usage:
    python export_csv.py
"""

import csv
import gzip
import json
import os
import re
import sys
import time

# ─── Path setup ───

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
REPO_ROOT = os.path.join(SCRIPT_DIR, "..")
DATA_DIR = os.path.join(REPO_ROOT, "data")
MOLECULES_DIR = os.path.join(DATA_DIR, "molecules")
OUTPUT_CSV = os.path.join(REPO_ROOT, "cccbdb-selected-energy.csv")
OUTPUT_JSON_GZ = os.path.join(REPO_ROOT, "cccbdb-molecules.json.gz")

from db_utils import ELEMENT_Z, count_atoms, count_electrons

# ─── Ranking tables ───

# Best energy: prefer highest-ranked method, then largest basis.
# =FULL variants (all-electron correlation) always rank above frozen-core.
# QCISD(TQ) includes perturbative quadruples — ranked above QCISD(T) but
# below CCSD (size-consistency advantage of coupled cluster).
METHOD_RANK = {
    "CCSD(T)=FULL": 16, "CCSD(T)": 15,
    "CCSD=FULL": 14, "CCSD": 13,
    "QCISD(TQ)=FULL": 12, "QCISD(TQ)": 11,
    "QCISD(T)=FULL": 10, "QCISD(T)": 9,
    "QCISD=FULL": 8, "QCISD": 7,
    "MP4=FULL": 6, "MP4": 5,
    "MP3=FULL": 4, "MP3": 3,
    "MP2=FULL": 2, "MP2": 1,
}

BASIS_RANK = {
    # Quadruple-zeta
    "daug-cc-pVQZ": 22,
    "aug-cc-pVQZ": 21, "aug-cc-p(Q+d)Z": 21,
    "cc-pVQZ": 20, "cc-pV(Q+d)Z": 20, "cc-pCVQZ": 20,
    # Triple-zeta
    "daug-cc-pVTZ": 19,
    "aug-cc-pVTZ": 18, "aug-cc-pCVTZ": 18, "aug-cc-pV(T+d)Z": 17,
    "cc-pVTZ": 16, "cc-pV(T+d)Z": 16, "cc-pCVTZ": 16, "Sadlej_pVTZ": 16,
    "6-311+G(3df,2pd)": 14, "6-311+G(3df,2p)": 13,
    # Double-zeta
    "daug-cc-pVDZ": 12,
    "aug-cc-pVDZ": 12,
    "cc-pVDZ": 11, "cc-pV(D+d)Z": 11, "cc-pCVDZ": 11,
    "6-31+G**": 9, "6-31G(2df,p)": 8,
    "6-311G**": 7, "6-311G*": 6, "TZVP": 7,
    "6-31G**": 5, "6-31G*": 4, "6-31G": 3,
    "3-21G*": 2, "3-21G": 1, "STO-3G": 0,
}

def make_species_id(formula, casno, charge):
    charge = int(charge)
    if charge > 0:
        suffix = "+" * charge
    elif charge < 0:
        suffix = "-" * abs(charge)
    else:
        suffix = ""
    return f"{formula}{suffix}_{casno}"


def find_best_energy(energy_data):
    """Find best correlated energy: highest-ranked method, largest basis.

    Returns (method, basis, energy) or (None, None, None).
    """
    if not energy_data or not isinstance(energy_data, dict):
        return None, None, None

    methods = energy_data.get("methods", {})
    best = None  # (method, basis, energy, m_rank, b_rank)

    for method_name, bases in methods.items():
        m_rank = METHOD_RANK.get(method_name, -1)
        if m_rank < 0:
            continue
        if not isinstance(bases, dict):
            continue
        for basis_name, value in bases.items():
            if value is None:
                continue
            b_rank = BASIS_RANK.get(basis_name, -1)
            if best is None or (m_rank, b_rank) > (best[3], best[4]):
                best = (method_name, basis_name, value, m_rank, b_rank)

    if best:
        return best[0], best[1], best[2]
    return None, None, None

# Cheapest method ordering: HF is cheapest, then MP2, MP3, etc.
CHEAPEST_METHOD_ORDER = [
    "HF",
    "MP2", "MP2=FULL", "ROMP2",
    "MP3", "MP3=FULL",
    "CID", "CISD",
    "MP4", "MP4=FULL",
    "CCD", "CCSD", "CCSD=FULL",
    "QCISD",
    "CCSD(T)", "CCSD(T)=FULL",
    "QCISD(T)", "QCISD(T)=FULL",
    "QCISD(TQ)", "QCISD(TQ)=FULL",
]

def count_atoms_from_xyz(xyz_str):
    """Count atoms from xyz coordinate string (one atom per line)."""
    if not xyz_str:
        return None
    lines = [l for l in xyz_str.strip().split('\n') if l.strip()]
    return len(lines)

def count_electrons_from_xyz(xyz_str, charge):
    """Count electrons from xyz coordinate elements minus charge."""
    if not xyz_str:
        return None
    total = 0
    for line in xyz_str.strip().split('\n'):
        parts = line.split()
        el = parts[0] if parts else None
        z = ELEMENT_Z.get(el) if el else None
        if z is None:
            return None
        total += z
    return total - int(charge)

def find_cheapest_energy(energy_data):
    """Find cheapest energy: lowest-cost method with smallest available basis.

    Tries HF first, then MP2, MP3, etc.
    Returns (method, basis, energy) or (None, None, None).
    """
    if not energy_data or not isinstance(energy_data, dict):
        return None, None, None

    methods = energy_data.get("methods", {})
    for method_name in CHEAPEST_METHOD_ORDER:
        bases = methods.get(method_name, {})
        if not isinstance(bases, dict):
            continue
        cheapest = None  # (basis, energy, rank)
        for basis_name, value in bases.items():
            if value is None:
                continue
            b_rank = BASIS_RANK.get(basis_name, 999)
            if cheapest is None or b_rank < cheapest[2]:
                cheapest = (basis_name, value, b_rank)
        if cheapest:
            return method_name, cheapest[0], cheapest[1]
    return None, None, None

def main():
    if not os.path.exists(MOLECULES_DIR):
        print(f"Error: molecules directory not found: {MOLECULES_DIR}")
        sys.exit(1)

    t0 = time.time()

    # Load species list by scanning molecule directories
    species_list = []
    for dirname in sorted(os.listdir(MOLECULES_DIR)):
        json_path = os.path.join(MOLECULES_DIR, dirname, "molecule.json")
        if not os.path.isfile(json_path):
            continue
        with open(json_path, "r", encoding="utf-8") as f:
            mol = json.load(f)
        species_list.append({
            "formula": mol["formula"],
            "name": mol.get("name", ""),
            "casno": mol["casno"],
            "charge": str(mol.get("charge", 0)),
        })

    # CSV columns (20 total)
    columns = [
        "species_id", "formula", "name", "casno", "charge",
        "point_group", "multiplicity", "s_squared", "closed_shell",
        "n_electrons", "n_atoms",
        "geometry_source", "geometry_method", "xyz",
        "cheapest_energy_method", "cheapest_energy_basis", "cheapest_energy",
        "best_energy_method", "best_energy_basis", "best_energy",
    ]

    out_rows = []
    missing = 0

    for sp in species_list:
        formula = sp["formula"]
        casno = sp["casno"]
        charge = int(sp["charge"])
        species_id = make_species_id(formula, casno, charge)
        dir_name = make_species_id(formula, casno, charge)

        # Load molecule.json
        json_path = os.path.join(MOLECULES_DIR, dir_name, "molecule.json")
        if not os.path.exists(json_path):
            missing += 1
            continue
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)

        # Point group
        point_group = data.get("point_group_cccbdb")

        # Spin
        spin = data.get("spin")
        multiplicity = None
        s_squared = None
        closed_shell = None
        if spin and isinstance(spin, dict):
            multiplicity = spin.get("multiplicity")
            s_squared = spin.get("S_squared")
            closed_shell = spin.get("closed_shell")

        # closed_shell: write 1/0 instead of True/False for CSV compatibility
        closed_shell_val = 1 if closed_shell else (0 if closed_shell is not None else None)

        # Geometry
        geoms = data.get("geometries", {})
        best_available = geoms.get("best_available")
        geom_source = ""
        geom_method = ""
        xyz = ""

        if best_available == "experimental":
            geom_source = "experimental"
            exp = geoms.get("experimental", {})
            if isinstance(exp, dict):
                xyz = exp.get("xyz", "")
        elif best_available and best_available.startswith("calculated/"):
            geom_source = "calculated"
            parts = best_available.split("/", 2)
            if len(parts) == 3:
                m, b = parts[1], parts[2]
                geom_method = f"{m}/{b}"
                calc = geoms.get("calculated", {})
                if isinstance(calc, dict):
                    method_dict = calc.get(m, {})
                    if isinstance(method_dict, dict):
                        basis_dict = method_dict.get(b, {})
                        if isinstance(basis_dict, dict):
                            xyz = basis_dict.get("xyz", "")

        # Computed properties — priority: properties block → xyz → formula
        props = data.get("properties", {})
        n_atoms_val = (props.get("n_atoms")
                       or count_atoms_from_xyz(xyz)
                       or count_atoms(formula))
        n_electrons = (props.get("n_electrons")
                       or count_electrons_from_xyz(xyz, charge)
                       or count_electrons(formula, charge))

        # Single-atom xyz fill
        if not xyz and n_atoms_val == 1:
            el_match = re.match(r'([A-Z][a-z]?)', formula.replace('+', '').replace('-', ''))
            if el_match:
                xyz = f"{el_match.group(1)} 0.000000 0.000000 0.000000"
                if not geom_source:
                    geom_source = "single atom"

        # Energy data
        energy = data.get("energy")

        # Cheapest energy (HF/smallest basis preferred)
        ce_method, ce_basis, ce_energy = find_cheapest_energy(energy)

        # Best correlated energy
        be_method, be_basis, be_energy = find_best_energy(energy)

        out_rows.append([
            species_id, formula, sp["name"], casno, charge,
            point_group, multiplicity, s_squared, closed_shell_val,
            n_electrons, n_atoms_val,
            geom_source, geom_method, xyz,
            ce_method or "", ce_basis or "", ce_energy,
            be_method or "", be_basis or "", be_energy,
        ])

    # Sort by species_id for deterministic output
    out_rows.sort(key=lambda r: r[0])

    # Write CSV
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(out_rows)

    elapsed = time.time() - t0
    csv_size = os.path.getsize(OUTPUT_CSV) / (1024 * 1024)

    print(f"cccbdb-selected-energy.csv: {len(out_rows):,} rows ({csv_size:.1f} MB)")
    if missing:
        print(f"  ({missing} species missing molecule.json)")
    print(f"Done in {elapsed:.1f}s")

def export_combined_json():
    """Export all molecule.json files as a single gzipped JSON.

    Produces cccbdb-molecules.json.gz at the repo root — a compact,
    gzip-compressed JSON mapping species_id → molecule data for all species.
    """
    if not os.path.exists(MOLECULES_DIR):
        print(f"Error: molecules directory not found: {MOLECULES_DIR}")
        return

    t0 = time.time()
    combined = {}

    for dirname in sorted(os.listdir(MOLECULES_DIR)):
        json_path = os.path.join(MOLECULES_DIR, dirname, "molecule.json")
        if not os.path.isfile(json_path):
            continue
        with open(json_path, "r", encoding="utf-8") as f:
            mol = json.load(f)
        combined[dirname] = mol

    json_bytes = json.dumps(combined, separators=(",", ":"), ensure_ascii=False).encode("utf-8")

    with gzip.open(OUTPUT_JSON_GZ, "wb", compresslevel=9) as f:
        f.write(json_bytes)

    elapsed = time.time() - t0
    gz_size = os.path.getsize(OUTPUT_JSON_GZ) / (1024 * 1024)
    raw_size = len(json_bytes) / (1024 * 1024)

    print(f"cccbdb-molecules.json.gz: {len(combined):,} species "
          f"({gz_size:.1f} MB gzip, {raw_size:.1f} MB uncompressed)")
    print(f"Done in {elapsed:.1f}s")

if __name__ == "__main__":
    main()
    export_combined_json()
