#!/usr/bin/env python3
"""
cccbdb_fetch_geometry.py

Purpose:
    Fetch molecular data from CCCBDB (Computational Chemistry Comparison and
    Benchmark Database, https://cccbdb.nist.gov) for all species in the data/molecules/ directory.

    For each species, fetches:
    - Experimental geometry (from expgeom2x.asp)
    - Calculated geometry fallback (from geom2x.asp + geom3x.asp)
    - Spin/multiplicity info (from spin2x.asp)

    Saves per-molecule folders with JSON data and archived raw HTML pages.

Input:  data/molecules/*/molecule.json
Output: data/molecules/<species>_<casno>/
            molecule.json    — parsed data
            expgeom2x.html   — raw experimental geometry page
            geom2x.html      — raw calc geometry summary page
            geom3x.html      — raw calculated geometry page
            spin2x.html      — raw spin page

Usage:
    python cccbdb_fetch_geometry.py              # fetch all missing
    python cccbdb_fetch_geometry.py --limit 10   # fetch at most 10
    python cccbdb_fetch_geometry.py --stats       # just print cache stats
    python cccbdb_fetch_geometry.py --calc-fallback          # update JSONs missing geometry
    python cccbdb_fetch_geometry.py --calc-fallback --limit 5  # update at most 5

CCCBDB Session Mechanics & Pitfalls:
    CCCBDB is an ASP.NET application behind Cloudflare. Several non-obvious
    behaviors were discovered during development:

    1. **Session cookies are mandatory for calculated geometry pages.**
       `geom3x.asp` (calculated Cartesian coordinates) takes only `method` and
       `basis` params — no `casno`. It relies on the server-side ASP session to
       know WHICH molecule to return coordinates for. The session is established
       by first visiting `expgeom2x.asp?casno=X&charge=Y`. Without a prior
       visit, geom3x returns a 500 Internal Server Error.

    2. **curl fails where Python urllib succeeds.**
       Even with exported cookie files (`-b cookies.txt`), curl requests to
       `geom3x.asp` consistently return 500 errors with Cloudflare
       `challenge-platform` scripts. Python's `http.cookiejar.CookieJar` +
       `HTTPCookieProcessor` handles the session correctly. The hypothesis is
       that Cloudflare's challenge token (set via JavaScript on the initial
       page load) is required and Python's opener propagates it properly while
       curl's cookie file export does not capture it.

    3. **Pages use inconsistent HTML tag closing.**
       `expgeom2x.asp` (experimental geometry) uses UNCLOSED `<TD>` tags:
           <TR><TD>O1<TD class="r">0.000<TD class="r">0.000<TD class="r">0.117
       `geom3x.asp` (calculated geometry) uses CLOSED `</TD>` tags:
           <TR><TD>O1</TD><TD class="num">0.000</TD>...
       The coordinate regex must handle both: `(?:</TD>)?` after each value.

    4. **Single atoms never have geometry data.**
       Species with only one atom (H, He, Li, ..., including ions like H+, He-)
       have no Cartesian coordinate tables on CCCBDB. Attempting to fetch their
       geom2x.asp summary wastes a request. The --calc-fallback mode detects
       single-atom formulas and skips geometry fetching (spin-only).

    5. **Spin page coverage is incomplete.**
       CCCBDB's spin2x.asp returns useful data (closed shell text or S^2 tables)
       for many species, but NOT all. Charged atoms and some exotic species
       return pages with no spin information at all. To avoid retrying these on
       every run, spin=None (never attempted) is distinguished from spin={}
       (attempted, no data found).

    6. **Some species return transient 500 errors.**
       A small fraction of molecules (~1-2%) trigger 500 errors on geom2x.asp
       or geom3x.asp even with valid sessions. These appear to be server-side
       issues (possibly missing data or computation failures). The code handles
       these gracefully (prints ERROR, continues to next species).

    7. **geom2x.asp table has DFT methods mixed with ab initio.**
       The "Moller Plesset perturbation" category includes DFT hybrid methods
       (B2PLYP, mPW2PLYP) alongside pure MP2/MP3. These are detected by
       checking `<acronym title="">` for "DFT" or "Hybrid" keywords, plus an
       explicit blocklist. Missing this causes DFT methods to be selected as
       "best" when they should be skipped.

    8. **Table row discrimination in geom2x.asp.**
       The summary table (`id="table2"`) has both data rows and header-repeat
       rows (repeated basis set headers mid-table). Data rows ALWAYS have a
       `bgcolor` attribute; header rows do not. This is the reliable way to
       distinguish them — using tag names or content patterns fails on edge
       cases.
"""

import argparse
import csv
import http.cookiejar
import json
import math
import os
import re
import time
from datetime import date
from urllib.parse import urlencode
from urllib.request import Request, build_opener, HTTPCookieProcessor, urlopen
from urllib.error import URLError, HTTPError

from bs4 import BeautifulSoup
from db_utils import _parse_formula, count_atoms, count_electrons, hill_formula

CACHE_DIR = "../data/molecules"
CCCBDB_BASE = "https://cccbdb.nist.gov"
REQUEST_DELAY = 2.0  # seconds between requests

# Categories to skip entirely when selecting calculated geometry.
# "density functional" covers all pure DFT methods (BLYP, B3LYP, PBE, etc.).
# "semi-empirical" and "composite" are low-accuracy categories.
DFT_CATEGORIES = {'density functional', 'semi-empirical', 'composite'}

# Individual methods to skip that appear in non-DFT categories.
# B2PLYP and mPW2PLYP are double-hybrid DFT methods listed under
# "Moller Plesset perturbation" in CCCBDB's table, NOT under "density functional".
DFT_METHODS = {'B2PLYP', 'B2PLYP=FULLultrafine', 'mPW2PLYP'}

# --- Validation helpers ---

def _validate_xyz_elements(xyz, formula):
    """Check if xyz string element counts match formula element counts.

    Normalizes D/T → H for comparison (CCCBDB uses 'H' in xyz even for
    deuterium/tritium species).

    Returns True if elements match, False if mismatch.
    """
    if not xyz or not formula:
        return True  # can't validate without both

    xyz_counts = {}
    for line in xyz.strip().split('\n'):
        parts = line.split()
        if parts:
            el = parts[0]
            el_norm = 'H' if el in ('D', 'T') else el
            xyz_counts[el_norm] = xyz_counts.get(el_norm, 0) + 1

    formula_counts_raw = _parse_formula(formula)
    formula_counts = {}
    for el, cnt in formula_counts_raw.items():
        el_norm = 'H' if el in ('D', 'T') else el
        formula_counts[el_norm] = formula_counts.get(el_norm, 0) + cnt

    return xyz_counts == formula_counts

def _correct_xyz_elements_from_internal_coords(xyz, formula, internal_coords):
    """Attempt to correct wrong element labels in xyz using internal coordinates.

    When element validation fails (e.g., OTi labeled as Ti,Ti in Cartesian table),
    the internal coordinates description (e.g., "rO=Ti") encodes the correct
    element identities for each atom index. This function extracts that mapping
    and reassigns element labels in the xyz string.

    Returns corrected xyz string if successful, or None if correction is not possible.
    """
    if not xyz or not formula or not internal_coords:
        return None

    formula_counts = _parse_formula(formula)
    # Normalize D/T → H for comparison
    formula_norm = {}
    for el, cnt in formula_counts.items():
        el_n = 'H' if el in ('D', 'T') else el
        formula_norm[el_n] = formula_norm.get(el_n, 0) + cnt

    lines = xyz.strip().split('\n')
    n_atoms = len(lines)

    # Build atom-index → element mapping from internal coordinate descriptions.
    # Descriptions like "rOTi", "rO=Ti", "aHOH" encode element symbols.
    # The atoms field gives indices [1,2] for distances, [1,2,3] for angles, etc.
    atom_elements = {}  # 1-based index → element symbol

    all_entries = (internal_coords.get('distances', [])
                   + internal_coords.get('angles', [])
                   + internal_coords.get('dihedrals', []))

    for entry in all_entries:
        desc = entry.get('description', '')
        atoms = entry.get('atoms', [])
        if not desc or not atoms:
            continue

        # Strip the prefix (r, a, d) and the = bond notation
        body = desc[1:].replace('=', '')

        # Parse element symbols from the description body.
        # E.g., "OTi" → [O, Ti], "HOH" → [H, O, H], "HCCS" → [H, C, C, S]
        desc_elements = []
        i = 0
        while i < len(body):
            if i + 1 < len(body) and body[i + 1].islower():
                desc_elements.append(body[i:i + 2])
                i += 2
            elif body[i].isupper():
                desc_elements.append(body[i])
                i += 1
            else:
                i += 1  # skip unexpected characters

        if len(desc_elements) != len(atoms):
            continue

        for el, idx in zip(desc_elements, atoms):
            if idx in atom_elements and atom_elements[idx] != el:
                return None  # conflicting assignments
            atom_elements[idx] = el

    # Check if we have assignments for all atoms
    if len(atom_elements) < n_atoms:
        return None

    # Verify the corrected element counts match the formula
    corrected_counts = {}
    for idx in range(1, n_atoms + 1):
        el = atom_elements.get(idx)
        if not el:
            return None
        el_n = 'H' if el in ('D', 'T') else el
        corrected_counts[el_n] = corrected_counts.get(el_n, 0) + 1

    if corrected_counts != formula_norm:
        return None

    # Build corrected xyz string
    corrected_lines = []
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) < 4:
            return None
        new_el = atom_elements[i + 1]
        corrected_lines.append(f"{new_el} {parts[1]} {parts[2]} {parts[3]}")

    return '\n'.join(corrected_lines)

# --- Species helpers (local, no pyscf dependency) ---

def _build_species(formula, charge):
    """Build canonical species string: Hill formula + charge symbols."""
    clean = hill_formula(formula)
    c = int(charge)
    if c > 0:
        return clean + "+" * c
    elif c < 0:
        return clean + "-" * abs(c)
    return clean

# --- Cache I/O ---

def _molecule_dir(species, casno):
    """Return per-molecule directory path: {CACHE_DIR}/{species}_{casno}/"""
    return os.path.join(CACHE_DIR, f"{species}_{casno}")

def cache_path(species, casno):
    """Return JSON cache file path: {species}_{casno}/molecule.json"""
    return os.path.join(_molecule_dir(species, casno), "molecule.json")

def _old_cache_path(casno, charge):
    """Old-style cache path for migration: {casno}_{charge}.json"""
    return os.path.join(CACHE_DIR, f"{casno}_{int(charge)}.json")

def is_cached(species, casno):
    return os.path.exists(cache_path(species, casno))

def load_cached(species, casno):
    path = cache_path(species, casno)
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    return None

def save_cached(data):
    species = _build_species(data["formula"], data["charge"])
    mol_dir = _molecule_dir(species, data["casno"])
    os.makedirs(mol_dir, exist_ok=True)
    path = os.path.join(mol_dir, "molecule.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

def save_html(species, casno, page_name, html):
    """Save raw HTML page to per-molecule folder for archival.

    Args:
        species: Hill formula species string (e.g. "H2O")
        casno: CAS number string
        page_name: base filename without extension (e.g. "expgeom2x", "spin2x")
        html: raw HTML string
    """
    if not html:
        return
    mol_dir = _molecule_dir(species, casno)
    os.makedirs(mol_dir, exist_ok=True)
    path = os.path.join(mol_dir, f"{page_name}.html")
    with open(path, "w", encoding="utf-8") as f:
        f.write(html)

def _migrate_filenames():
    """Migrate old cache formats into per-molecule folder structure.

    Handles two legacy formats:
    1. Ancient: ``{casno}_{charge}.json`` (both parts numeric)
    2. Previous: ``{species}_{casno}.json`` (flat file, letters in name)

    Target format (current):
        ``{species}_{casno}/molecule.json`` (per-molecule folder)

    Detection:
    - Ancient: both parts of ``stem.rsplit('_', 1)`` are numeric
    - Previous: any ``.json`` file directly in CACHE_DIR (not in a subfolder)

    Called automatically at startup before any cache operations. Idempotent.
    """
    if not os.path.exists(CACHE_DIR):
        return 0
    migrated = 0
    for fname in list(os.listdir(CACHE_DIR)):
        if not fname.endswith('.json'):
            continue
        fpath = os.path.join(CACHE_DIR, fname)
        if not os.path.isfile(fpath):
            continue

        try:
            with open(fpath, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except (json.JSONDecodeError, IOError):
            continue

        species = _build_species(data['formula'], data['charge'])
        mol_dir = _molecule_dir(species, data['casno'])
        new_path = os.path.join(mol_dir, "molecule.json")

        if fpath != new_path and not os.path.exists(new_path):
            os.makedirs(mol_dir, exist_ok=True)
            os.rename(fpath, new_path)
            migrated += 1
    return migrated

# --- Session handling ---

def create_session():
    """Create a cookie-aware URL opener for CCCBDB session management.

    CCCBDB uses ASP.NET server-side sessions to track which molecule the user
    is viewing. Pages like geom3x.asp and spin2x.asp do NOT accept a casno
    parameter — they rely entirely on the session state established by a prior
    visit to expgeom2x.asp. A new opener (= new CookieJar = new session) must
    be created for each molecule to avoid cross-contamination.

    We use Python's http.cookiejar rather than curl because Cloudflare's
    challenge tokens are only properly propagated through Python's
    HTTPCookieProcessor. curl with exported cookie files consistently fails
    on geom3x.asp with 500 errors.
    """
    cj = http.cookiejar.CookieJar()
    return build_opener(HTTPCookieProcessor(cj))

def fetch_with_session(opener, url):
    """Fetch a URL using session opener. Returns HTML string or None."""
    req = Request(url, headers={"User-Agent": "CCCBDB-research-fetcher/1.0"})
    try:
        with opener.open(req, timeout=30) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except (URLError, HTTPError, TimeoutError, OSError) as e:
        print(f"  ERROR fetching {url}: {e}")
        return None

def fetch_page(url):
    """Fetch a URL and return the HTML content as string (no session)."""
    req = Request(url, headers={"User-Agent": "CCCBDB-research-fetcher/1.0"})
    try:
        with urlopen(req, timeout=30) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except (URLError, HTTPError, TimeoutError, OSError) as e:
        print(f"  ERROR fetching {url}: {e}")
        return None

# --- Geometry parsing ---

def parse_atom_symbol(label):
    """Strip trailing digits from atom labels: 'Li1' → 'Li', 'H2' → 'H'."""
    return re.sub(r'\d+$', '', label)

def parse_geometry_page(html, formula=None, **_kwargs):
    """Parse CCCBDB geometry page (experimental or calculated).

    Handles two CCCBDB page types with different HTML conventions:
    - expgeom2x.asp (experimental): uses UNCLOSED <TD> tags.
      Example: ``<TR><TD>O1<TD class="r">0.000<TD class="r">0.000``
    - geom3x.asp (calculated): uses CLOSED </TD> tags.
      Example: ``<TR><TD>O1</TD><TD class="num">0.000</TD>``

    The coordinate regex uses ``(?:</TD>)?`` after each value to handle both.

    Parser fixes applied:
    - Dummy atom filter: atoms with symbol 'X' are skipped (computational
      reference points, not real atoms; fixes D3N).
    - Doubled conformer detection: if an atom label repeats, only the first
      conformer's coordinates are kept (fixes C2H6S, C3H5F, C2H2F4, H2S3,
      H3NS).
    - Single-atom element validation: for single-atom species, validates the
      extracted element against the formula and corrects mismatches (fixes
      Br/Br+/Br- where CCCBDB labels the atom as 'D' instead of 'Br').
    - Diatomic bond-length fallback: if only one atom is extracted for a
      diatomic species, reconstructs the second atom's position from the
      internal coordinates table bond length (fixes BrP where Br's
      z-coordinate is ``&nbsp;``).

    Args:
        html: Raw HTML string of the geometry page.
        formula: Molecular formula for validation/fallback. Optional.

    Returns dict with:
        xyz: "Symbol x y z\\n..." string (or None if no coordinate table found)
        point_group: string (or None)
        n_atoms: int (0 if no coordinates)
    """
    soup = BeautifulSoup(html, "html.parser")
    result = {"xyz": None, "point_group": None, "n_atoms": 0}

    # Find point group — in <H2>Point Group C<sub>2v</sub></H2>
    for h2 in soup.find_all("h2"):
        text = h2.get_text(strip=True)
        if text.lower().startswith("point group"):
            pg_text = text.replace("Point Group", "").replace("point group", "").strip()
            if pg_text:
                result["point_group"] = pg_text
            break

    # Find coordinate table using regex — handles both closed and unclosed <TD> tags
    coord_lines = []
    seen_labels = set()  # Track labels for doubled conformer detection
    atom_row_pattern = re.compile(
        r'<TR>\s*<TD>\s*([A-Z][a-z]?\d+)\s*(?:</TD>)?\s*'
        r'<TD[^>]*>\s*(-?[\d.]+)\s*(?:</TD>)?\s*'
        r'<TD[^>]*>\s*(-?[\d.]+)\s*(?:</TD>)?\s*'
        r'<TD[^>]*>\s*(-?[\d.]+)',
        re.IGNORECASE
    )
    for m in atom_row_pattern.finditer(html):
        label, x_s, y_s, z_s = m.groups()

        # Doubled conformer detection: if atom label repeats (e.g. S1 seen
        # twice), the second set is a different conformer — stop collecting.
        if label in seen_labels:
            break
        seen_labels.add(label)

        try:
            symbol = parse_atom_symbol(label)

            # Dummy atom filter: skip computational reference points (e.g. "X")
            if symbol == 'X':
                continue

            x, y, z = float(x_s), float(y_s), float(z_s)
            coord_lines.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")
        except ValueError:
            continue

    # Single-atom element validation: correct wrong element labels.
    # Fixes Br/Br+/Br- where CCCBDB labels the atom as "D1" instead of "Br1".
    if len(coord_lines) == 1 and formula:
        formula_counts = _parse_formula(formula)
        if sum(formula_counts.values()) == 1:
            formula_element = list(formula_counts.keys())[0]
            formula_el_norm = 'H' if formula_element in ('D', 'T') else formula_element
            xyz_parts = coord_lines[0].split()
            xyz_el_norm = 'H' if xyz_parts[0] in ('D', 'T') else xyz_parts[0]
            if xyz_el_norm != formula_el_norm:
                coord_lines[0] = f"{formula_element} {' '.join(xyz_parts[1:])}"

    # Diatomic bond-length fallback: if only one atom extracted for a diatomic,
    # reconstruct second atom from internal coordinates bond length.
    # Fixes BrP where Br's z-coordinate is &nbsp; in the Cartesian table.
    if len(coord_lines) == 1 and formula and count_atoms(formula) == 2:
        bond = _parse_internal_bond_length(html)
        if bond:
            _, _, distance = bond
            formula_counts = _parse_formula(formula)
            xyz_el = coord_lines[0].split()[0]
            xyz_el_norm = 'H' if xyz_el in ('D', 'T') else xyz_el
            # Find the remaining element after removing the extracted one
            remaining = {}
            found = False
            for el, cnt in formula_counts.items():
                el_norm = 'H' if el in ('D', 'T') else el
                if el_norm == xyz_el_norm and not found:
                    cnt -= 1
                    found = True
                if cnt > 0:
                    remaining[el] = cnt
            if remaining:
                second_el = list(remaining.keys())[0]
                coord_lines.append(
                    f"{second_el} 0.000000 0.000000 {distance:.6f}")

    if coord_lines:
        result["xyz"] = "\n".join(coord_lines)
        result["n_atoms"] = len(coord_lines)

    return result

def _parse_internal_bond_length(html):
    """Extract first bond length from internal coordinates table.

    Used as fallback for diatomic species where the Cartesian table is
    incomplete (e.g. BrP where Br's z-coordinate is ``&nbsp;``).

    Returns (element1, element2, distance_angstrom) or None.
    """
    soup = BeautifulSoup(html, "html.parser")

    # Find "Internal coordinates" heading
    ic_header = None
    for h2 in soup.find_all("h2"):
        if "internal coordinates" in h2.get_text().lower():
            ic_header = h2
            break
    if not ic_header:
        return None

    table = ic_header.find_next("table")
    if not table:
        return None

    for tr in table.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) < 2:
            continue
        desc = tds[0].get_text(strip=True)
        if desc.startswith('r'):
            try:
                value = float(tds[1].get_text(strip=True))
                # Extract element symbols from description (e.g. "rPBr" → P, Br)
                elements_str = desc[1:]  # Remove 'r' prefix
                els = re.findall(r'[A-Z][a-z]?', elements_str)
                if len(els) >= 2:
                    return (els[0], els[1], value)
                elif len(els) == 1:
                    return (els[0], els[0], value)  # homonuclear
            except ValueError:
                continue
    return None

def parse_internal_coordinates(html):
    """Parse internal coordinates from expgeom2x.asp HTML page.

    Extracts bond distances, angles, and dihedral angles from the
    "Internal coordinates" table. Three coordinate types distinguished
    by description prefix:
        r = bond distance (Angstroms), uses atoms 1-2
        a = bond angle (degrees), uses atoms 1-2-3
        d = dihedral angle (degrees), uses atoms 1-2-3-4

    Handles doubled conformers: if a description (e.g., "rCS") is seen
    a second time, stops collecting (takes first conformer only).

    Returns dict with distances/angles/dihedrals lists, or None if
    no internal coordinates table found.
    """
    soup = BeautifulSoup(html, "html.parser")

    # Find "Internal coordinates" heading
    ic_header = None
    for h2 in soup.find_all("h2"):
        if "internal coordinates" in h2.get_text().lower():
            ic_header = h2
            break
    if not ic_header:
        return None

    table = ic_header.find_next("table")
    if not table:
        return None

    distances = []
    angles = []
    dihedrals = []
    seen_descriptions = set()

    for tr in table.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) < 6:  # Description, Value, 4 atom columns minimum
            continue

        desc = tds[0].get_text(strip=True)
        if not desc:
            continue

        # Conformer detection: stop on repeated description
        if desc in seen_descriptions:
            break
        seen_descriptions.add(desc)

        try:
            value = float(tds[1].get_text(strip=True))
        except (ValueError, IndexError):
            continue

        # Extract atom indices (columns 2-5)
        atoms = []
        for td in tds[2:6]:
            atom_text = td.get_text(strip=True)
            if atom_text:
                try:
                    atoms.append(int(atom_text))
                except ValueError:
                    pass

        # Extract reference and comment (columns 6-7, if present)
        reference = tds[6].get_text(strip=True) if len(tds) > 6 else ""
        comment = tds[7].get_text(strip=True) if len(tds) > 7 else ""

        entry = {
            "description": desc,
            "value": value,
            "atoms": atoms,
            "reference": reference,
            "comment": comment,
        }

        if desc.startswith('r'):
            distances.append(entry)
        elif desc.startswith('a'):
            angles.append(entry)
        elif desc.startswith('d'):
            dihedrals.append(entry)

    if not distances and not angles and not dihedrals:
        return None

    return {
        "distances": distances,
        "angles": angles,
        "dihedrals": dihedrals,
    }

# --- Calculated geometry summary parsing ---

def _is_dft_method(entry):
    """Check if a method entry is DFT or DFT-hybrid.

    CCCBDB's geom2x.asp table organizes methods by category, but DFT hybrids
    appear in misleading places. Specifically, B2PLYP and mPW2PLYP are listed
    under "Moller Plesset perturbation" (MP section) even though they are
    double-hybrid DFT methods. Without this check, CCSD(T) might be skipped in
    favor of these DFT methods that appear later in the MP category.

    Three-layer detection:
    1. Category-level: entire "density functional" / "semi-empirical" / "composite" sections
    2. Method blocklist: known DFT hybrids that hide in non-DFT categories
    3. Title inspection: <acronym title="..."> containing "DFT" or "Hybrid"
    """
    if entry['category'] in DFT_CATEGORIES:
        return True
    if entry['method_name'] in DFT_METHODS:
        return True
    title = entry.get('method_title', '').lower()
    if 'dft' in title or 'hybrid' in title:
        return True
    return False

def parse_calc_geom_summary(html):
    """Parse the calculated geometry summary page (geom2x.asp).

    Extracts all available method/basis combinations from ``table2``
    ("Methods with standard basis sets").

    HTML structure of geom2x.asp table2:
        - First row: 2 empty <td> + N basis <th class="nowrap"><acronym>
        - Data rows ALWAYS have a ``bgcolor`` attribute; mid-table header
          repeats do NOT. This is the only reliable discriminator — earlier
          attempts using tag types or content patterns failed on edge cases.
        - Category headers: ``<th rowspan="N">`` with category name (e.g.,
          "hartree fock", "density functional", "coupled cluster").
        - Method cells: ``<th class="nowrap"><acronym title="...">METHOD</acronym>``
        - Available combos: ``<td><a href="geom3x.asp?method=M&basis=B">geom</a>``
        - Unavailable: ``<td>&nbsp;</td>`` or ``<td><abbr>dnf</abbr></td>``

    Edge cases handled:
        - html=None: returns [] (happens when geom2x.asp returns 500 for
          some exotic species like GeCl4)
        - No table2 on page: returns [] (some species have no calculated data)
        - Column count mismatch: falls back to ``basis_{id}`` naming

    Returns list of dicts ordered as they appear in the table:
        [{'method_name': 'CCSD(T)', 'method_id': 12, 'basis_name': '6-31+G**',
          'basis_id': 4, 'category': 'coupled cluster', 'col_index': 6,
          'method_title': '...'}, ...]
    """
    if not html:
        return []
    soup = BeautifulSoup(html, "html.parser")
    table = soup.find("table", id="table2")
    if not table:
        return []

    rows = table.find_all("tr")
    if not rows:
        return []

    # Parse basis set names from first header row
    basis_names = []
    for th in rows[0].find_all("th", class_="nowrap"):
        acronym = th.find("acronym")
        if acronym:
            basis_names.append(acronym.get_text(strip=True))

    if not basis_names:
        return []

    available = []
    current_category = None

    for tr in rows[1:]:
        # Only process data rows (always have bgcolor attribute)
        if not tr.get("bgcolor"):
            continue

        # Check for category header <th rowspan="N">
        category_th = tr.find("th", attrs={"rowspan": True})
        if category_th:
            current_category = category_th.get_text(strip=True).lower()

        # Find method name in <th class="nowrap"><acronym>...</acronym></th>
        method_th = tr.find("th", class_="nowrap")
        if not method_th:
            continue
        acronym = method_th.find("acronym")
        if not acronym:
            continue

        method_name = acronym.get_text(strip=True)
        method_title = acronym.get("title", "")

        # All <td> cells map 1:1 to basis columns
        tds = tr.find_all("td")
        for col, td in enumerate(tds):
            a_tag = td.find("a")
            if not a_tag or "geom3x" not in a_tag.get("href", ""):
                continue

            href = a_tag["href"]
            m_match = re.search(r'method=(\d+)', href)
            b_match = re.search(r'basis=(\d+)', href)
            if not m_match or not b_match:
                continue

            method_id = int(m_match.group(1))
            basis_id = int(b_match.group(1))
            basis_name = basis_names[col] if col < len(basis_names) else f"basis_{basis_id}"

            available.append({
                'method_name': method_name,
                'method_title': method_title,
                'method_id': method_id,
                'basis_name': basis_name,
                'basis_id': basis_id,
                'category': current_category or '',
                'col_index': col,
            })

    return available

def select_best_calc_geometry(available):
    """Select the best non-DFT calculated geometry.

    Strategy: highest-level theory (last in table) with largest basis set
    (rightmost column).

    The CCCBDB table order is:
        HF → DFT → MP2/MP3 → (DFT hybrids) → QCISD → CCD/CCSD/CCSD(T)
    So the last non-DFT entry is typically the highest-level theory available
    (usually CCSD(T) or CCSD). Among entries for that method, the rightmost
    column corresponds to the largest basis set.

    In practice, this produces results like:
        - CCSD(T)=FULL/aug-cc-pVQZ (best case, large molecules)
        - CCSD(T)=FULL/6-31G* (common fallback)
        - CCD/3-21G* (minimal, but still non-DFT)
        - MP3=FULL/6-31+G** (when coupled cluster is unavailable)

    Returns dict with method/basis info, or None if no non-DFT entry exists.
    """
    if not available:
        return None

    filtered = [e for e in available if not _is_dft_method(e)]
    if not filtered:
        return None

    # Last method in table order = highest theory level
    best_method = filtered[-1]['method_name']

    # Among entries for this method, pick rightmost column (largest basis)
    best = max(
        (e for e in filtered if e['method_name'] == best_method),
        key=lambda e: e['col_index']
    )

    return best

# --- Spin parsing ---

def _mult_from_s2(s2, n_electrons=None):
    """Convert S² expectation value to spin multiplicity (2S+1).

    If n_electrons is provided, enforces electron parity: odd electrons require
    even multiplicity (doublet, quartet, ...) and vice versa. When the raw S²
    gives a parity-forbidden multiplicity (e.g., spin contamination giving
    S²≈1.5 → triplet for an odd-electron species), the nearest allowed
    multiplicity is returned instead.

    Bug found 2026-02-12: CP (21 electrons) had S²=1.514 from spin contamination
    → raw formula gave triplet (3), but odd electrons can only be doublet/quartet.
    """
    S = (-1 + math.sqrt(1 + 4 * s2)) / 2
    mult = max(1, round(2 * S + 1))

    if n_electrons is not None:
        # Odd electrons → even mult; even electrons → odd mult
        need_even_mult = (n_electrons % 2 == 1)
        if need_even_mult and mult % 2 != 0:
            # Round to nearest even: e.g., 3→2 (if S²<2) or 3→4 (if S²>2)
            mult_lo = mult - 1 if (mult - 1) >= 2 else 2
            mult_hi = mult + 1
            # Pick whichever is closer to the raw 2S+1
            raw_mult = 2 * S + 1
            mult = mult_lo if abs(raw_mult - mult_lo) <= abs(raw_mult - mult_hi) else mult_hi
        elif not need_even_mult and mult % 2 == 0:
            mult_lo = mult - 1 if (mult - 1) >= 1 else 1
            mult_hi = mult + 1
            raw_mult = 2 * S + 1
            mult = mult_lo if abs(raw_mult - mult_lo) <= abs(raw_mult - mult_hi) else mult_hi

    return mult

def parse_spin_page(html, n_electrons=None):
    """Parse CCCBDB spin2x.asp page to extract spin multiplicity.

    CCCBDB's spin page has three distinct formats depending on the molecule:

    1. **Closed shell** (most common):
       Plain text: "H2O (Water) is closed shell. S²=0."
       → Parse via string matching for "closed shell".

    2. **Open shell with inline value**:
       Text with HTML superscript: "S<sup>2</sup>=0.75"
       → Parse via regex on raw HTML (not soup.get_text(), which loses tags).
       Multiplicity derived from S²: S = (-1 + sqrt(1 + 4*S²)) / 2, mult = 2S+1.

    3. **Open shell with method table**:
       A table of S² values for each computational method (HF, MP2, etc.).
       → Extract all numeric values from ``<td class="num">`` cells, take the
       most common (mode) after rounding to 4 decimal places. This handles
       slight numerical variation across methods.

    Args:
        html: Raw HTML of spin2x.asp page.
        n_electrons: If provided, used for electron parity validation.
            Odd electrons → only even multiplicities allowed (doublet, quartet, ...).
            Even electrons → only odd multiplicities allowed (singlet, triplet, ...).
            Prevents spin contamination from producing parity-forbidden multiplicities.

    Coverage limitations:
        - Many charged atoms/ions (H+, He-, Be+, etc.) have no spin data on
          CCCBDB. The page loads but contains no "closed shell" text, no
          S<sup>2</sup> value, and no numeric table. Returns None in this case.
        - The caller stores None (never attempted) vs {} (attempted, no data)
          to avoid retrying species that genuinely lack spin information.

    Returns dict: {S_squared, closed_shell, multiplicity, source_url} or None.
    """
    if not html:
        return None

    soup = BeautifulSoup(html, "html.parser")
    body_text = soup.get_text()

    # Closed shell: "FORMULA (NAME) is closed shell. S²=0."
    # Guard: if page shows "() is closed shell", the session was not properly
    # established and no species was selected. Return None to avoid recording
    # a false singlet for open-shell molecules. (Bug found 2026-02-12: this
    # caused n-Propyl radical C3H7 to be incorrectly stored as singlet.)
    if '() is closed shell' in body_text:
        return None
    if 'closed shell' in body_text.lower():
        # Electron parity check: odd electrons cannot be singlet.
        # CCCBDB incorrectly reports "closed shell" for some odd-electron
        # species (e.g. CNZn with 43 electrons). Return None so the caller
        # stores {} (attempted, no valid data).
        if n_electrons is not None and n_electrons % 2 == 1:
            return None
        return {
            'S_squared': 0.0,
            'closed_shell': True,
            'multiplicity': 1,
            'source_url': f'{CCCBDB_BASE}/spin2x.asp',
        }

    # Open shell with inline S² value: "S<sup>2</sup>=VALUE"
    s2_match = re.search(r'S<sup>2</sup>\s*=\s*([\d.]+)', html)
    if s2_match:
        s2 = float(s2_match.group(1))
        mult = _mult_from_s2(s2, n_electrons)
        return {
            'S_squared': s2,
            'closed_shell': False,
            'multiplicity': mult,
            'source_url': f'{CCCBDB_BASE}/spin2x.asp',
        }

    # Open shell with table of S² values per method.
    #
    # Strategy: find the highest-quality method that has data, then pick the
    # after-annihilation S² value from the largest basis set available.
    #
    # Each table cell contains two numbers separated by <br/>:
    #   - Upper: S² before annihilation (spin-contaminated)
    #   - Lower: S² after annihilation (more reliable)
    # We always use the after-annihilation (second) value.
    #
    # Method hierarchy (best first):
    METHOD_RANK = [
        'CCSD(T)=FULL', 'CCSD(T)', 'CCSD=FULL', 'CCSD',
        'QCISD(T)=FULL', 'QCISD(T)',
        'MP4=FULL', 'MP4', 'MP3=FULL', 'MP3',
        'MP2=FULL', 'MP2',
        'ROHF', 'HF',
    ]

    best_s2 = None
    best_method = None
    best_rank = len(METHOD_RANK)  # worse than any ranked method

    for tr in soup.find_all("tr"):
        tds = tr.find_all(["td", "th"])
        if not tds:
            continue
        method_name = tds[0].get_text(strip=True)

        # Skip DFT methods (B3LYP, BLYP, PBE, etc.) and category headers
        # that don't have direct S² data cells with <br>
        rank = None
        for i, m in enumerate(METHOD_RANK):
            if method_name == m:
                rank = i
                break
        if rank is None:
            continue
        if rank >= best_rank:
            continue  # already have a better method

        # Extract after-annihilation S² values from this row's data cells.
        # Cells use <br/> to separate before (upper) and after (lower) values.
        # We want the LAST non-empty cell (largest basis set) and its SECOND value.
        row_s2 = None
        for td in tds[1:]:
            br = td.find("br")
            if not br:
                continue
            # Get children text nodes around the <br>
            children = list(td.children)
            nums = []
            for child in children:
                t = str(child).strip()
                try:
                    nums.append(float(t))
                except ValueError:
                    continue
            if len(nums) >= 2:
                # Second value = after annihilation
                row_s2 = nums[1]
            elif len(nums) == 1:
                row_s2 = nums[0]

        if row_s2 is not None:
            best_s2 = row_s2
            best_method = method_name
            best_rank = rank

    if best_s2 is not None:
        mult = _mult_from_s2(best_s2, n_electrons)
        return {
            'S_squared': best_s2,
            'closed_shell': best_s2 < 0.01,
            'multiplicity': mult,
            'method': best_method,
            'source_url': f'{CCCBDB_BASE}/spin2x.asp',
        }

    return None

# --- Migration ---

def migrate_json(data):
    """Convert old flat JSON format to new nested structure.

    Old format (from initial fetch, pre-restructure):
        {"casno": "...", "geometry_xyz": "O 0.0 0.0 0.1\\n...",
         "point_group_cccbdb": "C2v", ...}

    New format:
        {"casno": "...", "geometries": {"experimental": {...}, "calculated": {},
         "best_available": "experimental"}, "spin": null, ...}

    Idempotent: if 'geometries' key already exists, returns data unchanged.
    """
    if 'geometries' in data:
        return data  # already new format

    new = {
        '_schema_version': '2.0',
        '_source_database': {
            'name': 'CCCBDB',
            'full_name': 'Computational Chemistry Comparison and Benchmark Database',
            'version': 'Release 22 (May 2022)',
            'url': 'https://cccbdb.nist.gov',
            'license': 'Public domain (U.S. Government work, NIST)',
        },
        'casno': data['casno'],
        'charge': data['charge'],
        'formula': data['formula'],
        'name': data['name'],
        'fetched_date': data.get('fetched_date'),
        'identifiers': {},
        'properties': {},
        'spin': None,
        'point_group_cccbdb': data.get('point_group_cccbdb'),
        'geometries': {
            'coordinate_units': 'angstrom',
            'experimental': None,
            'calculated': {},
            'best_available': None,
        },
    }

    if data.get('geometry_xyz'):
        new['geometries']['experimental'] = {
            'xyz': data['geometry_xyz'],
            'point_group': data.get('point_group_cccbdb'),
            'n_atoms': data['geometry_xyz'].count('\n') + 1,
            'source_url': data.get('source_url'),
            'geometry_type': 'equilibrium',
        }
        new['geometries']['best_available'] = 'experimental'

    return new

# --- Best geometry extraction (for stats and external use) ---

def get_best_geometry(data):
    """Extract the best available geometry XYZ string from a data dict.

    Handles both old flat format (``geometry_xyz`` key) and new nested format
    (``geometries.best_available`` pointer). This allows backward compatibility
    during migration when some JSONs may still be in old format.

    Resolution order for new format:
        1. Follow ``geometries.best_available`` pointer (e.g., "experimental"
           or "calculated/CCSD(T)/6-31G*")
        2. If pointer is missing/invalid, try experimental directly
        3. If no experimental, try any calculated geometry

    Returns xyz string or None.
    """
    # New format
    if 'geometries' in data:
        geoms = data['geometries']
        best = geoms.get('best_available')
        if best == 'experimental' and geoms.get('experimental'):
            return geoms['experimental'].get('xyz')
        if best and best.startswith('calculated/'):
            parts = best.split('/', 2)
            if len(parts) == 3:
                _, method, basis = parts
                return (geoms.get('calculated', {})
                        .get(method, {})
                        .get(basis, {})
                        .get('xyz'))
        # Fallback: check any available geometry
        if geoms.get('experimental') and geoms['experimental'].get('xyz'):
            return geoms['experimental']['xyz']
        for method_dict in geoms.get('calculated', {}).values():
            for basis_dict in method_dict.values():
                if basis_dict.get('xyz'):
                    return basis_dict['xyz']
        return None

    # Old flat format
    return data.get('geometry_xyz')

# --- Fetch full species data ---

def fetch_species(formula, name, casno, charge):
    """Fetch all available data for a single species from CCCBDB.

    Request sequence (order matters due to session state):
        1. expgeom2x.asp?casno=X&charge=Y — experimental geometry.
           **Must be first**: establishes the ASP session for this molecule.
           All subsequent session-dependent pages (geom2x, geom3x, spin2x)
           will return data for THIS molecule based on the session cookie.
        2. geom2x.asp — calculated geometry summary (only if no experimental).
           Session-dependent: returns methods available for the molecule from step 1.
        3. geom3x.asp?method=M&basis=B — calculated Cartesian coordinates.
           Session-dependent: no casno param, relies on session from step 1.
        4. spin2x.asp — spin/multiplicity info.
           Session-dependent: returns S² data for the molecule from step 1.

    A new session (opener) is created per molecule to prevent cross-contamination.
    Each request is separated by REQUEST_DELAY seconds to be respectful.

    Returns a dict with the nested JSON structure (see module docstring).
    """
    opener = create_session()
    species = _build_species(formula, charge)

    exp_url = f"{CCCBDB_BASE}/expgeom2x.asp?{urlencode({'casno': casno, 'charge': charge})}"

    data = {
        'casno': casno,
        'charge': int(charge),
        'formula': formula,
        'name': name,
        'fetched_date': str(date.today()),
        'identifiers': {},
        'properties': {},
        'spin': None,
        'point_group_cccbdb': None,
        'geometries': {
            'coordinate_units': 'angstrom',
            'experimental': None,
            'calculated': {},
            'best_available': None,
        },
    }

    # 1. Fetch experimental geometry (also sets session cookie)
    html = fetch_with_session(opener, exp_url)
    save_html(species, casno, "expgeom2x", html)
    if not html:
        return data

    parsed = parse_geometry_page(html, formula=formula)
    data['point_group_cccbdb'] = parsed['point_group']

    # Store internal coordinates if available
    internal = parse_internal_coordinates(html)
    if internal:
        data['geometries']['experimental_internal_coords'] = internal

    if parsed['xyz']:
        # Element count validation: correct or discard if elements don't match formula
        xyz_to_store = parsed['xyz']
        valid = _validate_xyz_elements(xyz_to_store, formula)
        if not valid and internal:
            # Attempt to correct element labels using internal coordinate descriptions
            corrected = _correct_xyz_elements_from_internal_coords(
                xyz_to_store, formula, internal)
            if corrected:
                xyz_to_store = corrected
                valid = True
        if valid:
            data['geometries']['experimental'] = {
                'xyz': xyz_to_store,
                'point_group': parsed['point_group'],
                'n_atoms': parsed['n_atoms'],
                'source_url': exp_url,
                'geometry_type': 'equilibrium',
            }
            data['geometries']['best_available'] = 'experimental'

    # 2. If no experimental geometry, try calculated fallback
    if not parsed['xyz']:
        time.sleep(REQUEST_DELAY)
        summary_url = f"{CCCBDB_BASE}/geom2x.asp"
        summary_html = fetch_with_session(opener, summary_url)
        save_html(species, casno, "geom2x", summary_html)

        if summary_html:
            avail = parse_calc_geom_summary(summary_html)
            best = select_best_calc_geometry(avail)

            if best:
                time.sleep(REQUEST_DELAY)
                calc_url = (f"{CCCBDB_BASE}/geom3x.asp?"
                            f"method={best['method_id']}&basis={best['basis_id']}")
                calc_html = fetch_with_session(opener, calc_url)
                save_html(species, casno, "geom3x", calc_html)

                if calc_html:
                    calc_parsed = parse_geometry_page(calc_html, formula=formula)
                    if calc_parsed['xyz']:
                        method = best['method_name']
                        basis = best['basis_name']

                        data['geometries']['calculated'].setdefault(method, {})
                        data['geometries']['calculated'][method][basis] = {
                            'xyz': calc_parsed['xyz'],
                            'point_group': (calc_parsed['point_group']
                                            or data['point_group_cccbdb']),
                            'n_atoms': calc_parsed['n_atoms'],
                            'method_id': best['method_id'],
                            'basis_id': best['basis_id'],
                            'source_url': calc_url,
                        }
                        data['geometries']['best_available'] = (
                            f"calculated/{method}/{basis}")

                        if not data['point_group_cccbdb'] and calc_parsed['point_group']:
                            data['point_group_cccbdb'] = calc_parsed['point_group']

    # 3. Fetch spin info
    time.sleep(REQUEST_DELAY)
    spin_url = f"{CCCBDB_BASE}/spin2x.asp"
    spin_html = fetch_with_session(opener, spin_url)
    save_html(species, casno, "spin2x", spin_html)
    if spin_html:
        ne = count_electrons(formula, charge)
        data['spin'] = parse_spin_page(spin_html, n_electrons=ne)

    return data

# --- Species list ---

def load_species_list():
    """Load species list by scanning data/molecules/ directories.

    Each molecule.json contains formula, name, casno, charge at the top level.
    Returns list of dicts with keys: formula, name, casno, charge.
    """
    rows = []
    if not os.path.isdir(CACHE_DIR):
        return rows
    for dirname in sorted(os.listdir(CACHE_DIR)):
        json_path = os.path.join(CACHE_DIR, dirname, "molecule.json")
        if not os.path.isfile(json_path):
            continue
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        rows.append({
            "formula": data["formula"],
            "name": data.get("name", ""),
            "casno": data["casno"],
            "charge": str(data.get("charge", 0)),
        })
    return rows

# --- Stats ---

def print_stats(species_list):
    """Print cache statistics."""
    cached = 0
    has_exp_geom = 0
    has_calc_geom = 0
    no_geom = 0
    has_spin = 0
    errors = 0
    old_format = 0

    for row in species_list:
        species = _build_species(row["formula"], int(row["charge"]))
        data = load_cached(species, row["casno"])
        if not data:
            continue
        cached += 1

        if 'geometries' not in data:
            old_format += 1
            if data.get("geometry_xyz"):
                has_exp_geom += 1
            elif data.get("geometry_source") == "fetch_error":
                errors += 1
            else:
                no_geom += 1
        else:
            geoms = data['geometries']
            if geoms.get('experimental') and geoms['experimental'].get('xyz'):
                has_exp_geom += 1
            elif geoms.get('calculated'):
                has_calc_geom += 1
            else:
                no_geom += 1

            if data.get('spin') and data['spin'].get('multiplicity'):
                has_spin += 1

    total = len(species_list)
    missing = total - cached
    print(f"Archive stats ({total} species total):")
    print(f"  Archived:        {cached}")
    print(f"    Experimental:  {has_exp_geom}")
    print(f"    Calculated:    {has_calc_geom}")
    print(f"    No geom:       {no_geom}")
    print(f"    Errors:        {errors}")
    print(f"    Has spin:      {has_spin}")
    if old_format:
        print(f"    Old format:    {old_format}")
    print(f"  Not archived:    {missing}")

# --- Main ---

def main():
    parser = argparse.ArgumentParser(description="Fetch CCCBDB geometries")
    parser.add_argument("--limit", type=int, default=0,
                        help="Max species to fetch/update (0=all)")
    parser.add_argument("--stats", action="store_true",
                        help="Print cache stats only")
    parser.add_argument("--retry-errors", action="store_true",
                        help="Re-fetch entries with fetch_error")
    parser.add_argument("--calc-fallback", action="store_true",
                        help="Update cached entries: migrate format, add "
                             "calculated geometry fallback + spin data")
    parser.add_argument("--refetch-spin", action="store_true",
                        help="Re-fetch and re-parse spin data for all entries "
                             "with existing spin (corrected parser)")
    args = parser.parse_args()

    # Ensure CWD is the scripts/ directory so relative paths (../data/) resolve correctly
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    species_list = load_species_list()
    print(f"Loaded {len(species_list)} species from {CACHE_DIR}")

    os.makedirs(CACHE_DIR, exist_ok=True)

    # Auto-migrate flat files into per-molecule folders
    n_migrated = _migrate_filenames()
    if n_migrated:
        print(f"Migrated {n_migrated} cache files into per-molecule folders")

    if args.stats:
        print_stats(species_list)
        return

    if args.calc_fallback:
        _run_calc_fallback(species_list, args.limit)
        return

    if args.refetch_spin:
        _run_refetch_spin(species_list, args.limit)
        return

    # Normal mode: fetch all missing species
    fetched = 0
    skipped = 0

    for i, row in enumerate(species_list):
        casno = row["casno"]
        charge = int(row["charge"])
        formula = row["formula"]
        name = row["name"]
        species = _build_species(formula, charge)

        # Check cache
        if is_cached(species, casno):
            if args.retry_errors:
                existing = load_cached(species, casno)
                if existing and existing.get("geometry_source") != "fetch_error":
                    skipped += 1
                    continue
            else:
                skipped += 1
                continue

        if args.limit and fetched >= args.limit:
            print(f"\nReached --limit={args.limit}, stopping.")
            break

        print(f"  [{i+1}/{len(species_list)}] Fetching {formula} "
              f"(charge={charge}, CAS={casno})...")
        data = fetch_species(formula, name, casno, charge)
        save_cached(data)
        fetched += 1

        best = data['geometries'].get('best_available') or 'none'
        pg = data.get('point_group_cccbdb', '?')
        spin_str = (f"mult={data['spin']['multiplicity']}"
                    if data.get('spin') else "no spin")
        print(f"    → best={best}, PG={pg}, {spin_str}")

        time.sleep(REQUEST_DELAY)

    print(f"\nDone. Fetched: {fetched}, Skipped (cached): {skipped}")
    print_stats(species_list)

def _run_calc_fallback(species_list, limit):
    """Migrate all cached JSONs to new format, then fetch calculated
    geometry and spin data for entries that need it.

    Two-pass design:
        Pass 1 (fast, no network): Migrate any remaining old-format JSONs
        to the new nested structure. This is idempotent.

        Pass 2 (slow, network): For each species missing geometry or spin,
        create a session, visit expgeom2x (establish session), then fetch
        geom2x + geom3x (calculated geometry) and/or spin2x (spin data).

    Skip logic to avoid redundant fetches:
        - Species with both geometry and spin already populated → skip
        - Single-atom species (H, He, Li, etc.) → skip geometry fetch
          (atoms have no Cartesian coordinate tables on CCCBDB), only fetch spin
        - spin={} (empty dict) means "attempted, no data" → skip spin fetch.
          This is distinct from spin=None which means "never attempted".
          Without this distinction, species that genuinely lack CCCBDB spin
          data would be retried on every run.

    Estimated runtime: ~3.5 hours for full 2186-species database at 2s delay.
    Use --limit N to test on a subset. The process is safe to interrupt and
    resume — each species is saved immediately after processing.
    """

    # Pass 1: migrate all old-format JSONs (fast, no network)
    migrated = 0
    for row in species_list:
        species = _build_species(row["formula"], int(row["charge"]))
        existing = load_cached(species, row["casno"])
        if existing and 'geometries' not in existing:
            save_cached(migrate_json(existing))
            migrated += 1
    if migrated:
        print(f"Migrated {migrated} entries to new JSON format")

    # Pass 2: fetch calculated geometry + spin for entries that need it
    updated = 0
    skipped = 0

    for i, row in enumerate(species_list):
        if limit and updated >= limit:
            print(f"\nReached --limit={limit}, stopping fetches.")
            break

        casno = row["casno"]
        charge = int(row["charge"])
        species = _build_species(row["formula"], charge)
        data = load_cached(species, casno)
        if not data:
            skipped += 1
            continue

        has_geometry = (data.get('geometries', {}).get('best_available')
                        is not None)
        # spin=None means never attempted; spin={} means attempted but no data
        needs_spin = ('spin' not in data or data['spin'] is None)

        # Single-atom species never have geometry — only fetch spin
        formula = data['formula']
        clean = re.sub(r'[+-]', '', formula)
        atom_pairs = re.findall(r'([A-Z][a-z]?)(\d*)', clean)
        n_atoms = sum(int(n) if n else 1 for el, n in atom_pairs if el)
        is_single_atom = (n_atoms <= 1)

        needs_geometry = (not has_geometry and not is_single_atom)

        if not needs_geometry and not needs_spin:
            skipped += 1
            continue

        print(f"  [{i+1}/{len(species_list)}] Updating {formula} "
              f"(charge={charge}, CAS={casno})...")

        opener = create_session()

        # Visit expgeom page to establish session
        exp_url = (f"{CCCBDB_BASE}/expgeom2x.asp?"
                   f"{urlencode({'casno': casno, 'charge': charge})}")
        exp_html = fetch_with_session(opener, exp_url)
        save_html(species, casno, "expgeom2x", exp_html)
        time.sleep(REQUEST_DELAY)

        # Fetch calculated geometry if needed (skip for single atoms)
        if needs_geometry:
            summary_url = f"{CCCBDB_BASE}/geom2x.asp"
            summary_html = fetch_with_session(opener, summary_url)
            save_html(species, casno, "geom2x", summary_html)

            if summary_html:
                avail = parse_calc_geom_summary(summary_html)
                best = select_best_calc_geometry(avail)

                if best:
                    time.sleep(REQUEST_DELAY)
                    calc_url = (f"{CCCBDB_BASE}/geom3x.asp?"
                                f"method={best['method_id']}&basis={best['basis_id']}")
                    calc_html = fetch_with_session(opener, calc_url)
                    save_html(species, casno, "geom3x", calc_html)

                    if calc_html:
                        calc_parsed = parse_geometry_page(calc_html, formula=formula)
                        if calc_parsed['xyz']:
                            method = best['method_name']
                            basis = best['basis_name']

                            geoms = data['geometries']
                            geoms['calculated'].setdefault(method, {})
                            geoms['calculated'][method][basis] = {
                                'xyz': calc_parsed['xyz'],
                                'point_group': (calc_parsed['point_group']
                                                or data.get('point_group_cccbdb')),
                                'n_atoms': calc_parsed['n_atoms'],
                                'method_id': best['method_id'],
                                'basis_id': best['basis_id'],
                                'source_url': calc_url,
                            }
                            geoms['best_available'] = (
                                f"calculated/{method}/{basis}")

                            if (not data.get('point_group_cccbdb')
                                    and calc_parsed['point_group']):
                                data['point_group_cccbdb'] = calc_parsed['point_group']

            time.sleep(REQUEST_DELAY)

        # Fetch spin if needed
        if needs_spin:
            spin_url = f"{CCCBDB_BASE}/spin2x.asp"
            spin_html = fetch_with_session(opener, spin_url)
            save_html(species, casno, "spin2x", spin_html)
            ne = count_electrons(formula, charge)
            spin_result = parse_spin_page(spin_html, n_electrons=ne) if spin_html else None
            # Save result or empty dict (marks as attempted, prevents retry)
            data['spin'] = spin_result if spin_result else {}
            time.sleep(REQUEST_DELAY)

        data['fetched_date'] = str(date.today())
        save_cached(data)
        updated += 1

        best_str = data['geometries'].get('best_available') or 'none'
        spin_str = (f"mult={data['spin']['multiplicity']}"
                    if data.get('spin') else "no spin")
        print(f"    → best={best_str}, {spin_str}")

    print(f"\nDone. Updated: {updated}, Skipped: {skipped}")
    print_stats(species_list)

def _run_refetch_spin(species_list, limit):
    """Re-fetch and re-parse spin data for entries with existing spin.

    Targets entries where spin data was previously parsed with an older
    (potentially buggy) parser. Skips:
        - spin=None (never attempted)
        - spin={} (attempted, no data on CCCBDB)

    For each qualifying entry:
        1. Creates a new session
        2. Visits expgeom2x.asp (establishes session for the molecule)
        3. Fetches spin2x.asp
        4. Saves raw HTML to per-molecule folder
        5. Re-parses with corrected parser (METHOD_RANK + after-annihilation)
        6. Logs any multiplicity changes
        7. Saves updated JSON

    Safe to interrupt and resume — each species is saved immediately.
    """
    # Identify entries that need re-fetch
    candidates = []
    for row in species_list:
        species = _build_species(row["formula"], int(row["charge"]))
        data = load_cached(species, row["casno"])
        if not data:
            continue
        spin = data.get('spin')
        # Only re-fetch if spin has actual data (not None, not empty dict)
        if spin and isinstance(spin, dict) and spin.get('multiplicity'):
            candidates.append(row)

    print(f"Found {len(candidates)} entries with existing spin data to re-fetch")

    updated = 0
    changed = 0

    for i, row in enumerate(candidates):
        if limit and updated >= limit:
            print(f"\nReached --limit={limit}, stopping.")
            break

        casno = row["casno"]
        charge = int(row["charge"])
        formula = row["formula"]
        species = _build_species(formula, charge)
        data = load_cached(species, casno)

        old_spin = data.get('spin', {})
        old_mult = old_spin.get('multiplicity') if old_spin else None

        print(f"  [{i+1}/{len(candidates)}] Re-fetching spin for {formula} "
              f"(charge={charge}, CAS={casno})...")

        opener = create_session()

        # Establish session
        exp_url = (f"{CCCBDB_BASE}/expgeom2x.asp?"
                   f"{urlencode({'casno': casno, 'charge': charge})}")
        exp_html = fetch_with_session(opener, exp_url)
        save_html(species, casno, "expgeom2x", exp_html)
        time.sleep(REQUEST_DELAY)

        # Fetch spin page
        spin_url = f"{CCCBDB_BASE}/spin2x.asp"
        spin_html = fetch_with_session(opener, spin_url)
        save_html(species, casno, "spin2x", spin_html)

        if spin_html:
            ne = count_electrons(formula, charge)
            new_spin = parse_spin_page(spin_html, n_electrons=ne)
            if new_spin:
                new_mult = new_spin.get('multiplicity')
                data['spin'] = new_spin
                if old_mult != new_mult:
                    print(f"    CHANGED: mult {old_mult} → {new_mult}")
                    changed += 1
                else:
                    print(f"    unchanged: mult={new_mult}")
            else:
                # Parser returned None — mark as empty to avoid future retry
                data['spin'] = {}
                if old_mult:
                    print(f"    CHANGED: mult {old_mult} → None (no data)")
                    changed += 1
        else:
            print(f"    ERROR: could not fetch spin page")

        data['fetched_date'] = str(date.today())
        save_cached(data)
        updated += 1
        time.sleep(REQUEST_DELAY)

    print(f"\nDone. Re-fetched: {updated}, Multiplicity changed: {changed}")
    print_stats(species_list)

if __name__ == "__main__":
    main()
