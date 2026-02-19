"""
CCCBDB — Python wrapper for the NIST CCCBDB database mirror.

Usage:
    import sys; sys.path.insert(0, "scripts")
    from cccbdb import CCCBDB
    db = CCCBDB()                # auto-finds data/molecules/

    mol = db.get("H2O")         # lookup by formula
    mol = db.get("7732185")     # lookup by CAS number

    db.energy("H2O", "CCSD(T)", "cc-pVTZ")   # single energy value
    db.search(charge=0, closed_shell=True)    # filter species
    db.species()                              # all 2,186 species (lightweight)
"""

import json
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_DATA = os.path.join(_HERE, "..", "data")

def _build_species(formula, charge):
    """Build canonical species string: formula + charge symbols."""
    c = int(charge)
    if c > 0:
        return formula + "+" * c
    elif c < 0:
        return formula + "-" * abs(c)
    return formula

class CCCBDB:
    """Wrapper for CCCBDB mirror data (reads per-molecule JSON files)."""

    def __init__(self, data_dir=None):
        data_dir = data_dir or _DEFAULT_DATA
        self._data_dir = data_dir
        self._mol_dir = os.path.join(data_dir, "molecules")

        if not os.path.isdir(self._mol_dir):
            raise FileNotFoundError(
                f"Molecules directory not found: {self._mol_dir}\n"
                "Expected data/molecules/ in the repository."
            )

        # Build index by scanning molecule directories
        self._by_id = {}       # species_id -> record
        self._by_formula = {}  # formula -> [species_id, ...]
        self._by_casno = {}    # casno -> [species_id, ...]
        self._all_ids = []     # ordered list of species_ids

        for dirname in sorted(os.listdir(self._mol_dir)):
            json_path = os.path.join(self._mol_dir, dirname, "molecule.json")
            if not os.path.isfile(json_path):
                continue

            with open(json_path, encoding="utf-8") as jf:
                mol_data = json.load(jf)

            formula = mol_data["formula"]
            name = mol_data.get("name", "")
            casno = mol_data["casno"]
            charge = int(mol_data.get("charge", 0))
            species_str = _build_species(formula, charge)
            species_id = f"{species_str}_{casno}"

            # Extract spin data
            spin = mol_data.get("spin")
            multiplicity = None
            s_squared = None
            closed_shell = None
            if spin and isinstance(spin, dict) and spin:
                multiplicity = spin.get("multiplicity")
                s_squared = spin.get("S_squared")
                closed_shell = spin.get("closed_shell")

            # Extract point group
            point_group = mol_data.get("point_group_cccbdb")
            if not point_group:
                geom = mol_data.get("geometries", {})
                exp = geom.get("experimental")
                if exp and isinstance(exp, dict):
                    point_group = exp.get("point_group")

            record = {
                "species_id": species_id,
                "formula": formula,
                "name": name,
                "casno": casno,
                "charge": charge,
                "point_group": point_group,
                "multiplicity": multiplicity,
                "s_squared": s_squared,
                "closed_shell": closed_shell,
                "_dirname": dirname,
                "_mol_data": mol_data,
            }

            self._by_id[species_id] = record
            self._by_formula.setdefault(formula, []).append(species_id)
            self._by_casno.setdefault(casno, []).append(species_id)
            self._all_ids.append(species_id)

    # ------------------------------------------------------------------
    # get() — full record for one species
    # ------------------------------------------------------------------

    def get(self, identifier):
        """Look up a species by formula, CAS number, or species_id.

        Args:
            identifier: formula (e.g. "H2O"), CAS number (e.g. "7732185"),
                        or species_id (e.g. "H2O_7732185").
                        CAS numbers are detected as all-digit strings.
                        When multiple species share a formula or CAS number
                        (e.g. H2O neutral, cation, anion), prefers charge=0.

        Returns:
            dict with species info, geometry, and energies — or None if not found.
        """
        record = self._resolve(identifier)
        if record is None:
            return None

        mol_data = record["_mol_data"]

        result = {
            "species_id": record["species_id"],
            "formula": record["formula"],
            "name": record["name"],
            "casno": record["casno"],
            "charge": record["charge"],
            "point_group": record["point_group"],
            "multiplicity": record["multiplicity"],
            "s_squared": record["s_squared"],
            "closed_shell": record["closed_shell"],
        }

        # Geometry — best available
        result["geometry"] = self._best_geometry(mol_data)

        # Energies — flat dict {method/basis: value}
        result["energies"] = self._energies(mol_data)

        return result

    def _resolve(self, identifier):
        """Resolve an identifier to a record, preferring charge=0."""
        if "_" in identifier:
            return self._by_id.get(identifier)

        if identifier.isdigit():
            ids = self._by_casno.get(identifier, [])
        else:
            ids = self._by_formula.get(identifier, [])

        if not ids:
            return None

        # Prefer charge=0, then smallest |charge|
        candidates = [self._by_id[sid] for sid in ids]
        candidates.sort(key=lambda r: (abs(r["charge"]), r["charge"]))
        return candidates[0]

    def _best_geometry(self, mol_data):
        """Return best available geometry from molecule.json."""
        if not mol_data:
            return None

        geom = mol_data.get("geometries", {})
        best = geom.get("best_available")
        if not best:
            return None

        if best == "experimental":
            exp = geom.get("experimental")
            if exp and isinstance(exp, dict):
                return {
                    "source": "experimental",
                    "xyz": exp.get("xyz"),
                    "point_group": exp.get("point_group"),
                }
        elif best.startswith("calculated/"):
            # Format: "calculated/METHOD/BASIS"
            parts = best.split("/", 2)
            if len(parts) == 3:
                _, method, basis = parts
                calc = geom.get("calculated", {})
                method_data = calc.get(method, {})
                basis_data = method_data.get(basis, {})
                if basis_data:
                    return {
                        "source": f"calculated/{method}/{basis}",
                        "xyz": basis_data.get("xyz"),
                        "point_group": basis_data.get("point_group"),
                    }

        return None

    def _energies(self, mol_data):
        """Return {method/basis: energy} from molecule.json energy.methods."""
        if not mol_data:
            return {}
        energy = mol_data.get("energy")
        if not energy:
            return {}
        methods = energy.get("methods", {})
        result = {}
        for method, bases in methods.items():
            for basis, value in bases.items():
                result[f"{method}/{basis}"] = value
        return result

    # ------------------------------------------------------------------
    # energy() — single energy value
    # ------------------------------------------------------------------

    def energy(self, identifier, method, basis):
        """Return a single energy value (hartree) or None.

        Args:
            identifier: formula, CAS number, or species_id
            method: e.g. "CCSD(T)", "HF"
            basis: e.g. "cc-pVTZ", "STO-3G"
        """
        record = self._resolve(identifier)
        if record is None:
            return None
        mol_data = record["_mol_data"]
        if not mol_data:
            return None
        energy = mol_data.get("energy")
        if not energy:
            return None
        methods = energy.get("methods", {})
        return methods.get(method, {}).get(basis)

    # ------------------------------------------------------------------
    # search() — filter species
    # ------------------------------------------------------------------

    _SEARCH_FIELDS = {
        "formula", "name", "casno", "charge", "point_group",
        "multiplicity", "closed_shell",
    }

    def search(self, limit=None, **kwargs):
        """Filter species by field values.

        Supported kwargs: formula, name, casno, charge, point_group,
            multiplicity, closed_shell.

        Returns list of lightweight dicts (no energies — call get() for full data).
        """
        for key in kwargs:
            if key not in self._SEARCH_FIELDS:
                raise ValueError(f"Unknown search field: {key!r}")

        results = []
        for sid in self._all_ids:
            record = self._by_id[sid]
            match = True
            for key, val in kwargs.items():
                if record.get(key) != val:
                    match = False
                    break
            if match:
                results.append({
                    "species_id": record["species_id"],
                    "formula": record["formula"],
                    "name": record["name"],
                    "casno": record["casno"],
                    "charge": record["charge"],
                })
                if limit is not None and len(results) >= limit:
                    break

        return results

    # ------------------------------------------------------------------
    # species() — lightweight listing
    # ------------------------------------------------------------------

    def species(self):
        """Return all 2,186 species as lightweight dicts.

        Each dict has: species_id, formula, name, casno, charge.
        """
        return [
            {
                "species_id": self._by_id[sid]["species_id"],
                "formula": self._by_id[sid]["formula"],
                "name": self._by_id[sid]["name"],
                "casno": self._by_id[sid]["casno"],
                "charge": self._by_id[sid]["charge"],
            }
            for sid in self._all_ids
        ]
