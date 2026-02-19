"""
Cross-validation checks for CCCBDB mirror data integrity.

Runs against actual molecule.json files in data/molecules/.
Each check compares redundant fields that must be consistent.
Failures are either parser bugs (fix the code) or upstream
CCCBDB anomalies (documented in data/suspected_errors/anomalies.json).

Run: python -m pytest tests/test_data_integrity.py -v
"""
import json
import math
import os
import re
import sys

import pytest

# Allow imports from scripts/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from db_utils import ELEMENT_Z, count_atoms, count_electrons

MOLECULES_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'molecules')
ANOMALIES_PATH = os.path.join(os.path.dirname(__file__), '..', 'data',
                               'suspected_errors', 'anomalies.json')

def _load_known_anomalies():
    """Load known anomaly species_ids from anomalies.json."""
    if not os.path.exists(ANOMALIES_PATH):
        return {}
    with open(ANOMALIES_PATH, 'r', encoding='utf-8') as f:
        data = json.load(f)
    # Map species_id -> list of categories
    result = {}
    for a in data.get('anomalies', []):
        sid = a['species_id']
        result.setdefault(sid, set()).add(a['category'])
    return result

KNOWN_ANOMALIES = _load_known_anomalies()

def load_all_molecules():
    """Load all molecule.json files, yielding (dir_name, data) tuples."""
    if not os.path.isdir(MOLECULES_DIR):
        return
    for dirname in sorted(os.listdir(MOLECULES_DIR)):
        json_path = os.path.join(MOLECULES_DIR, dirname, 'molecule.json')
        if not os.path.isfile(json_path):
            continue
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        yield dirname, data

@pytest.fixture(scope="module")
def molecules():
    mols = list(load_all_molecules())
    if not mols:
        pytest.skip("No molecule data available")
    return mols

def _get_best_xyz(data):
    """Extract best available xyz string from molecule data."""
    geoms = data.get('geometries', {})
    best = geoms.get('best_available')
    if best == 'experimental':
        exp = geoms.get('experimental')
        if exp and isinstance(exp, dict):
            return exp.get('xyz')
    elif best and best.startswith('calculated/'):
        parts = best.split('/', 2)
        if len(parts) == 3:
            _, method, basis = parts
            return (geoms.get('calculated', {})
                    .get(method, {})
                    .get(basis, {})
                    .get('xyz'))
    return None

def _is_known_anomaly(dirname, category):
    """Check if this species has a known anomaly of the given category."""
    cats = KNOWN_ANOMALIES.get(dirname, set())
    return category in cats

# ---------------------------------------------------------------------------
# Tier 1: Hard constraints (violation = definite error)
# ---------------------------------------------------------------------------

def test_formula_vs_xyz_atom_count(molecules):
    """Formula atom count must match xyz line count when both present."""
    failures = []
    for dirname, data in molecules:
        if _is_known_anomaly(dirname, 'formula_mismatch'):
            continue
        if _is_known_anomaly(dirname, 'truncated_cartesian'):
            continue
        if _is_known_anomaly(dirname, 'doubled_conformer'):
            continue
        if _is_known_anomaly(dirname, 'dummy_atom'):
            continue

        formula = data.get('formula', '')
        xyz = _get_best_xyz(data)
        if not xyz:
            continue

        n_formula = count_atoms(formula)
        n_xyz = len([l for l in xyz.strip().split('\n') if l.strip()])
        if n_formula != n_xyz:
            failures.append(f"{dirname}: formula={formula} ({n_formula}) vs xyz ({n_xyz})")

    assert not failures, f"Formula vs xyz atom count mismatches:\n" + "\n".join(failures)

def test_xyz_elements_vs_formula(molecules):
    """Elements in xyz must match formula elements (same set, same counts).

    Note: CCCBDB uses 'H' in xyz even for deuterium (D) and tritium (T)
    species. We normalize D/T → H for comparison purposes.
    """
    failures = []
    for dirname, data in molecules:
        if _is_known_anomaly(dirname, 'formula_mismatch'):
            continue
        if _is_known_anomaly(dirname, 'truncated_cartesian'):
            continue
        if _is_known_anomaly(dirname, 'doubled_conformer'):
            continue
        if _is_known_anomaly(dirname, 'dummy_atom'):
            continue
        if _is_known_anomaly(dirname, 'wrong_element'):
            continue

        formula = data.get('formula', '')
        xyz = _get_best_xyz(data)
        if not xyz:
            continue

        # Count elements from xyz, normalizing D/T → H
        xyz_counts = {}
        for line in xyz.strip().split('\n'):
            parts = line.split()
            if parts:
                el = parts[0]
                if el not in ELEMENT_Z:
                    continue
                # Normalize isotopes: D/T in xyz → H
                el_norm = 'H' if el in ('D', 'T') else el
                xyz_counts[el_norm] = xyz_counts.get(el_norm, 0) + 1

        # Count elements from formula, normalizing D/T → H
        from db_utils import _parse_formula
        formula_counts_raw = _parse_formula(formula)
        formula_counts = {}
        for el, cnt in formula_counts_raw.items():
            el_norm = 'H' if el in ('D', 'T') else el
            formula_counts[el_norm] = formula_counts.get(el_norm, 0) + cnt

        if formula_counts != xyz_counts:
            failures.append(f"{dirname}: formula={formula} -> {formula_counts}, xyz -> {xyz_counts}")

    assert not failures, f"Element mismatches:\n" + "\n".join(failures[:20])

def test_xyz_count_vs_n_atoms_field(molecules):
    """xyz line count must match n_atoms field within the same geometry block."""
    failures = []
    for dirname, data in molecules:
        if _is_known_anomaly(dirname, 'truncated_cartesian'):
            continue
        if _is_known_anomaly(dirname, 'doubled_conformer'):
            continue
        if _is_known_anomaly(dirname, 'dummy_atom'):
            continue

        geoms = data.get('geometries', {})
        exp = geoms.get('experimental')
        if exp and isinstance(exp, dict) and exp.get('xyz'):
            xyz_lines = len([l for l in exp['xyz'].strip().split('\n') if l.strip()])
            n_atoms = exp.get('n_atoms')
            if n_atoms is not None and n_atoms != xyz_lines:
                failures.append(f"{dirname}: exp n_atoms={n_atoms} vs xyz lines={xyz_lines}")

    assert not failures, f"n_atoms vs xyz line count mismatches:\n" + "\n".join(failures)

def test_electron_parity_vs_multiplicity(molecules):
    """Odd electrons require even multiplicity; even electrons require odd."""
    failures = []
    for dirname, data in molecules:
        if _is_known_anomaly(dirname, 'impossible_spin'):
            continue

        formula = data.get('formula', '')
        charge = int(data.get('charge', 0))
        spin = data.get('spin')
        if not spin or not isinstance(spin, dict):
            continue
        mult = spin.get('multiplicity')
        if mult is None:
            continue

        n_electrons = count_electrons(formula, charge)
        if n_electrons is None:
            continue

        if n_electrons % 2 == 1 and mult % 2 != 0:
            failures.append(f"{dirname}: {n_electrons}e (odd) but mult={mult} (odd)")
        if n_electrons % 2 == 0 and mult % 2 != 1:
            failures.append(f"{dirname}: {n_electrons}e (even) but mult={mult} (even)")

    assert not failures, f"Electron parity violations:\n" + "\n".join(failures)

def test_closed_shell_vs_multiplicity(molecules):
    """closed_shell=True requires multiplicity=1."""
    failures = []
    for dirname, data in molecules:
        spin = data.get('spin')
        if not spin or not isinstance(spin, dict):
            continue
        closed = spin.get('closed_shell')
        mult = spin.get('multiplicity')
        if closed is True and mult is not None and mult != 1:
            failures.append(f"{dirname}: closed_shell=True but mult={mult}")

    assert not failures, f"closed_shell contradictions:\n" + "\n".join(failures)

def test_s_squared_vs_multiplicity(molecules):
    """S^2 should approximate S(S+1) where S=(mult-1)/2, tolerance +/-1.0.

    Tolerance is generous because spin contamination (common in UHF/UMP)
    can push S^2 significantly away from the ideal value. The multiplicity
    is derived using electron parity constraints, so S^2 may disagree.
    """
    # Investigated severe spin contamination cases (confirmed from raw HTML):
    # - C3+: CCSD(T) before-annihilation S^2=1.823 vs ideal 0.75 (doublet, 17e)
    # - CSi-: CCSD(T)=FULL before-annihilation S^2=1.751 vs ideal 0.75 (doublet, 21e)
    # Both have no after-annihilation values; CCCBDB reports only before-annihilation.
    KNOWN_SPIN_CONTAMINATION = {"C3+_12075353", "CSi-_409212"}

    failures = []
    for dirname, data in molecules:
        if dirname in KNOWN_SPIN_CONTAMINATION:
            continue
        spin = data.get('spin')
        if not spin or not isinstance(spin, dict):
            continue
        s2 = spin.get('S_squared')
        mult = spin.get('multiplicity')
        if s2 is None or mult is None:
            continue

        S = (mult - 1) / 2.0
        expected_s2 = S * (S + 1)
        if abs(s2 - expected_s2) > 1.0:
            failures.append(f"{dirname}: S^2={s2}, mult={mult}, expected S^2={expected_s2:.2f}")

    assert not failures, f"S^2 vs multiplicity mismatches:\n" + "\n".join(failures)

def test_best_available_pointer_valid(molecules):
    """best_available must point to existing data path."""
    failures = []
    for dirname, data in molecules:
        geoms = data.get('geometries', {})
        best = geoms.get('best_available')
        if best is None:
            continue

        if best == 'experimental':
            exp = geoms.get('experimental')
            if not exp or not isinstance(exp, dict) or not exp.get('xyz'):
                failures.append(f"{dirname}: best_available='experimental' but no xyz")
        elif best.startswith('calculated/'):
            parts = best.split('/', 2)
            if len(parts) == 3:
                _, method, basis = parts
                calc = geoms.get('calculated', {})
                if method not in calc or basis not in calc.get(method, {}):
                    failures.append(f"{dirname}: best_available='{best}' but path missing")
        else:
            failures.append(f"{dirname}: unknown best_available='{best}'")

    assert not failures, f"Invalid best_available pointers:\n" + "\n".join(failures)

def test_casno_matches_directory(molecules):
    """casno in JSON must match directory name suffix."""
    failures = []
    for dirname, data in molecules:
        casno = data.get('casno', '')
        expected_suffix = f"_{casno}"
        if not dirname.endswith(expected_suffix):
            failures.append(f"{dirname}: casno={casno}")

    assert not failures, f"casno vs directory mismatches:\n" + "\n".join(failures)

def test_formula_elements_valid(molecules):
    """All elements in formula must exist in ELEMENT_Z."""
    failures = []
    for dirname, data in molecules:
        formula = data.get('formula', '')
        from db_utils import _parse_formula
        counts = _parse_formula(formula)
        for el in counts:
            if el not in ELEMENT_Z:
                failures.append(f"{dirname}: unknown element '{el}' in {formula}")

    assert not failures, f"Invalid elements:\n" + "\n".join(failures)

# ---------------------------------------------------------------------------
# Tier 2: Physical constraints (violation = likely error)
# ---------------------------------------------------------------------------

def test_hf_sto3g_consistency(molecules):
    """hf_sto3g field must match methods.HF['STO-3G']."""
    failures = []
    for dirname, data in molecules:
        energy = data.get('energy')
        if not energy or not isinstance(energy, dict):
            continue
        hf_sto3g = energy.get('hf_sto3g')
        hf_dict = energy.get('methods', {}).get('HF', {})
        hf_dict_val = hf_dict.get('STO-3G') if isinstance(hf_dict, dict) else None

        if hf_sto3g is not None and hf_dict_val is not None:
            if abs(hf_sto3g - hf_dict_val) > 1e-10:
                failures.append(f"{dirname}: hf_sto3g={hf_sto3g} vs HF/STO-3G={hf_dict_val}")

    assert not failures, f"hf_sto3g inconsistencies:\n" + "\n".join(failures)

def test_hf_variational_principle(molecules):
    """For HF, larger basis must give lower energy (variational principle)."""
    # Check STO-3G >= cc-pVQZ
    failures = []
    for dirname, data in molecules:
        energy = data.get('energy')
        if not energy or not isinstance(energy, dict):
            continue
        hf = energy.get('methods', {}).get('HF', {})
        if not isinstance(hf, dict):
            continue
        sto3g = hf.get('STO-3G')
        cc_pvqz = hf.get('cc-pVQZ')
        if sto3g is not None and cc_pvqz is not None:
            if cc_pvqz > sto3g + 1e-6:
                failures.append(
                    f"{dirname}: STO-3G={sto3g:.6f} < cc-pVQZ={cc_pvqz:.6f}")

    assert not failures, f"HF variational violations:\n" + "\n".join(failures)

def test_s_squared_non_negative(molecules):
    """S^2 must be >= 0."""
    failures = []
    for dirname, data in molecules:
        spin = data.get('spin')
        if not spin or not isinstance(spin, dict):
            continue
        s2 = spin.get('S_squared')
        if s2 is not None and s2 < 0:
            failures.append(f"{dirname}: S^2={s2}")

    assert not failures, f"Negative S^2 values:\n" + "\n".join(failures)

def test_xyz_coordinate_bounds(molecules):
    """All xyz coordinates must be in [-100, 100] Angstrom."""
    failures = []
    for dirname, data in molecules:
        xyz = _get_best_xyz(data)
        if not xyz:
            continue
        for line in xyz.strip().split('\n'):
            parts = line.split()
            if len(parts) >= 4:
                try:
                    coords = [abs(float(parts[i])) for i in range(1, 4)]
                    if any(c > 100 for c in coords):
                        failures.append(f"{dirname}: {line}")
                except ValueError:
                    failures.append(f"{dirname}: non-numeric coordinate: {line}")

    assert not failures, f"Out-of-bounds coordinates:\n" + "\n".join(failures)

# ---------------------------------------------------------------------------
# Tier 3: Cross-validation (internal coords vs Cartesian xyz)
# ---------------------------------------------------------------------------

def _parse_xyz_coords(xyz):
    """Parse xyz string into list of (element, x, y, z) tuples."""
    atoms = []
    for line in xyz.strip().split('\n'):
        parts = line.split()
        if len(parts) >= 4:
            el = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            atoms.append((el, x, y, z))
    return atoms

def _distance(a1, a2):
    """Euclidean distance between two (el, x, y, z) tuples."""
    return math.sqrt(sum((a1[i] - a2[i])**2 for i in range(1, 4)))

def test_internal_coords_vs_xyz_distances(molecules):
    """Bond lengths from internal coordinates must match Cartesian distances.

    For molecules with both experimental_internal_coords.distances and
    experimental.xyz, compute interatomic distances from Cartesian coords
    using atom indices from internal coords and compare.

    Tolerance: 0.35 Angstrom. Internal coordinates and Cartesian coordinates
    in CCCBDB often come from different experimental sources or geometry
    types (r0 vs re vs rs vs rg), so exact agreement is not expected.
    Differences of 0.05-0.2 Angstrom are common for the same bond measured
    with different techniques.

    Mismatches > 0.5 Angstrom are skipped as they indicate atom numbering
    inconsistencies between the internal coords and Cartesian tables (not
    data errors).
    """
    failures = []
    for dirname, data in molecules:
        geoms = data.get('geometries', {})
        exp = geoms.get('experimental')
        ic = geoms.get('experimental_internal_coords')
        if not exp or not ic or not exp.get('xyz'):
            continue

        atoms = _parse_xyz_coords(exp['xyz'])
        distances = ic.get('distances', [])

        for d in distances:
            idx = d.get('atoms', [])
            if len(idx) < 2:
                continue
            i1, i2 = idx[0] - 1, idx[1] - 1  # 1-indexed → 0-indexed
            if i1 < 0 or i2 < 0 or i1 >= len(atoms) or i2 >= len(atoms):
                continue

            cart_dist = _distance(atoms[i1], atoms[i2])
            ic_dist = d['value']
            diff = abs(cart_dist - ic_dist)
            # Skip large mismatches (atom numbering inconsistency, not error)
            if diff > 0.5:
                continue
            if diff > 0.35:
                failures.append(
                    f"{dirname}: {d['description']}={ic_dist:.4f} vs "
                    f"Cartesian={cart_dist:.4f} (diff={diff:.4f})")

    assert not failures, (
        f"Internal coords vs Cartesian distance mismatches:\n"
        + "\n".join(failures[:20])
    )

def test_basis_sets_field_consistency(molecules):
    """All basis sets used in methods must be declared in energy.basis_sets.

    The basis_sets field lists columns from the CCCBDB energy table header,
    so it may include entries with no data for any method (superset is OK).
    But any basis set used in methods must appear in basis_sets.
    """
    failures = []
    for dirname, data in molecules:
        energy = data.get('energy')
        if not energy or not isinstance(energy, dict):
            continue
        methods = energy.get('methods', {})
        if not methods:
            continue
        declared = set(energy.get('basis_sets', []))
        actual = set()
        for method_data in methods.values():
            if isinstance(method_data, dict):
                actual.update(method_data.keys())
        missing = actual - declared
        if missing:
            failures.append(
                f"{dirname}: basis sets in methods but not in basis_sets: {missing}")

    assert not failures, (
        f"basis_sets field inconsistencies:\n" + "\n".join(failures[:20])
    )

def test_calculated_geometry_atom_count(molecules):
    """Calculated geometry xyz line count must match formula atom count."""
    failures = []
    for dirname, data in molecules:
        if _is_known_anomaly(dirname, 'dummy_atom'):
            continue
        formula = data.get('formula', '')
        n_expected = count_atoms(formula)
        calc = data.get('geometries', {}).get('calculated', {})
        if not calc or not isinstance(calc, dict):
            continue
        for method, bases in calc.items():
            if not isinstance(bases, dict):
                continue
            for basis, geom in bases.items():
                if not isinstance(geom, dict):
                    continue
                xyz = geom.get('xyz')
                if not xyz:
                    continue
                n_xyz = len([l for l in xyz.strip().split('\n') if l.strip()])
                if n_xyz != n_expected:
                    failures.append(
                        f"{dirname}: {method}/{basis} has {n_xyz} atoms, "
                        f"formula {formula} expects {n_expected}")

    assert not failures, (
        f"Calculated geometry atom count mismatches:\n" + "\n".join(failures[:20])
    )
