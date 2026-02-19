# CCCBDB Data Anomaly Log

This log documents all known data anomalies, upstream suspect data, data
completeness gaps, and parser design decisions that affect the mirror's output.
Each entry includes evidence, verification instructions, and (where applicable)
links to local HTML source files so reviewers can verify independently.

**Structured data**: See `anomalies.json` in this directory (covers
Sections 1-8 and 11; Sections 9-10 are documentation-only and not in the JSON).

## How to Verify

1. Open the linked HTML file in a browser or text editor
2. Search for the evidence described (line numbers, HTML patterns)
3. Compare against the molecule.json output to confirm the resolution
4. Run verification commands shown in each entry

---

## 1. Doubled Conformer Geometries (5 species)

**Category**: `doubled_conformer`
**Detection**: `count_atoms(formula)` != xyz line count (exactly 2x)
**Root cause**: CCCBDB experimental geometry pages for these molecules contain
two complete conformer coordinate sets concatenated in a single HTML table.
Atom numbering restarts at the second conformer boundary.

**Alternative hypothesis considered**: Could these be genuine dimers rather than
duplicated conformers? No — three independent pieces of evidence rule this out:
(1) atom labels restart with the same numbering (e.g., S1 appears twice),
indicating a second copy, not a continuation; (2) the formula is for the monomer,
not a dimer (e.g., C2H6S, not C4H12S2); (3) the atom count is exactly 2x the
formula count in all 5 cases.

**Parser fix**: Track seen atom labels during extraction. When a label repeats
(e.g., "S1" seen twice), stop collecting -- take only the first conformer.

### C2H6S -- Ethanethiol (CAS 75081)
- **HTML**: [data/molecules/C2H6S_75081/expgeom2x.html](../../data/molecules/C2H6S_75081/expgeom2x.html)
- **Evidence**: Coordinate table has 18 rows. Atom labels C1,C2,S3,H4..H9
  appear in the first 9 rows, then restart as C1,C2,S3,H4..H9 in rows 10-18.
  Formula C2H6S = 9 atoms expected.
- **Before fix**: xyz has 18 atoms (both conformers concatenated)
- **After fix**: xyz has 9 atoms (first conformer only)
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/C2H6S_75081/molecule.json')); print(d['geometries']['experimental']['n_atoms'])"`
  Expected output: `9`

### C3H5F -- Allyl Fluoride (CAS 818928)
- **HTML**: [data/molecules/C3H5F_818928/expgeom2x.html](../../data/molecules/C3H5F_818928/expgeom2x.html)
- **Evidence**: Coordinate table has 18 rows with label restart.
  Formula C3H5F = 9 atoms expected.
- **Before fix**: xyz has 18 atoms
- **After fix**: xyz has 9 atoms (first conformer only)
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/C3H5F_818928/molecule.json')); print(d['geometries']['experimental']['n_atoms'])"`
  Expected output: `9`

### C2H2F4 -- 1,1,2,2-Tetrafluoroethane (CAS 359353)
- **HTML**: [data/molecules/C2H2F4_359353/expgeom2x.html](../../data/molecules/C2H2F4_359353/expgeom2x.html)
- **Evidence**: Coordinate table has 16 rows with label restart.
  Formula C2H2F4 = 8 atoms expected.
- **Before fix**: xyz has 16 atoms
- **After fix**: xyz has 8 atoms (first conformer only)
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/C2H2F4_359353/molecule.json')); print(d['geometries']['experimental']['n_atoms'])"`
  Expected output: `8`

### H2S3 -- Trisulfane (CAS 13845233)
- **HTML**: [data/molecules/H2S3_13845233/expgeom2x.html](../../data/molecules/H2S3_13845233/expgeom2x.html)
- **Evidence**: Coordinate table has 10 rows with label restart.
  Formula H2S3 = 5 atoms expected.
- **Before fix**: xyz has 10 atoms
- **After fix**: xyz has 5 atoms (first conformer only)
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/H2S3_13845233/molecule.json')); print(d['geometries']['experimental']['n_atoms'])"`
  Expected output: `5`

### H3NS -- Thiohydroxylamine (CAS 14097008)
- **HTML**: [data/molecules/H3NS_14097008/expgeom2x.html](../../data/molecules/H3NS_14097008/expgeom2x.html)
- **Evidence**: Coordinate table has 10 rows with label restart.
  Formula H3NS = 5 atoms expected.
- **Before fix**: xyz has 10 atoms
- **After fix**: xyz has 5 atoms (first conformer only)
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/H3NS_14097008/molecule.json')); print(d['geometries']['experimental']['n_atoms'])"`
  Expected output: `5`

---

## 2. Incomplete Cartesian Coordinates (2 species)

### BrP -- Phosphorus monobromide (CAS 59727161)
- **HTML**: [data/molecules/BrP_59727161/expgeom2x.html](../../data/molecules/BrP_59727161/expgeom2x.html)
- **Alternative hypothesis considered**: Could this be a genuine single-atom
  species (i.e., is the formula wrong)? No — three confirmations: (1) the
  internal coordinates table on the same page provides a P-Br bond length
  (rPBr = 2.171 angstrom, reference 1979Col:1051), confirming two atoms;
  (2) CAS 59727161 is the registered identifier for BrP; (3) the calculated
  geometry (CCSD(T)/cc-pVTZ) shows two atoms with a similar P-Br distance
  (~2.08 angstrom).
- **Evidence**: Cartesian table: P1 has full coordinates (0,0,0),
  but Br2's z-coordinate is `&nbsp;` (non-breaking space) instead of a
  numeric value. Internal coordinates table has `rPBr = 2.171` angstrom
  with reference 1979Col:1051.
- **Before fix**: xyz has 1 atom (P only); Br dropped due to missing z-coord
- **After fix**: xyz has 2 atoms -- P at (0,0,0), Br at (0,0,2.171)
  reconstructed from internal coordinate bond length
- **Cross-validation**: Calculated geometry (CCSD(T)/cc-pVTZ) shows P-Br
  distance of ~2.08 angstrom, consistent with experimental 2.171 angstrom
  (expected difference for experimental vs computed bond length)
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/BrP_59727161/molecule.json')); print(d['geometries']['experimental']['xyz'])"`
  Expected: `P 0.000000 0.000000 0.000000\nBr 0.000000 0.000000 2.171000`

### C5H4O2 -- Furfural (CAS 980101)
- **HTML**: [data/molecules/C5H4O2_980101/expgeom2x.html](../../data/molecules/C5H4O2_980101/expgeom2x.html)
- **Evidence**: Page explicitly states "No coordinate data available".
  Coordinate table has only 1 atom (O1) of 11 expected.
- **Status**: `upstream_unrecoverable` -- no parser fix possible. Calculated
  geometry fallback is used as best_available.
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/C5H4O2_980101/molecule.json')); print(d['geometries']['best_available'])"`
  Expected: `calculated/CCD/cc-pVDZ` (or similar calculated path)

---

## 3. Wrong Element Labels (4 species)

### Br, Br+, Br- -- Bromine atom (CAS 10097322)
- **HTML**: [data/molecules/Br_10097322/expgeom2x.html](../../data/molecules/Br_10097322/expgeom2x.html)
- **Alternative hypothesis considered**: Could the formula be wrong (i.e., is
  this actually a deuterium atom)? No — three independent confirmations:
  (1) CAS number 10097322 is the registered identifier for elemental Bromine;
  (2) the CCCBDB page title says "Bromine atom"; (3) Br has Z=35 while D has
  Z=1 — a deuterium atom would have entirely different energy values, and the
  HF/STO-3G energy in the data (-2572.4 Hartree) is consistent with Br (Z=35),
  not D (Z=1, which would be around -0.5 Hartree).
- **Evidence**: Coordinate table labels the atom as "D1" instead of "Br1".
  This appears to be an HTML rendering bug in CCCBDB (the label field, not
  the data itself).
- **Before fix**: xyz = "D 0.000000 0.000000 0.000000" (wrong element)
- **After fix**: xyz = "Br 0.000000 0.000000 0.000000" (element corrected
  from formula)
- **Parser fix**: For single-atom species, validate extracted element against
  formula. If mismatch, use formula element.
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/Br_10097322/molecule.json')); print(d['geometries']['experimental']['xyz'])"`
  Expected: `Br 0.000000 0.000000 0.000000`

### OTi -- Titanium monoxide (CAS 12137201)
- **HTML**: [data/molecules/OTi_12137201/expgeom2x.html](../../data/molecules/OTi_12137201/expgeom2x.html)
- **Alternative hypothesis considered**: Could the formula be wrong (i.e., is
  this actually Ti2)? No — three independent confirmations: (1) CAS number
  12137201 is the registered identifier for titanium monoxide (TiO); (2) the
  internal coordinates table says "rO=Ti" not "rTi=Ti", confirming one atom
  is oxygen; (3) the calculated geometry (CCSD(T)/cc-pVTZ) correctly contains
  O and Ti elements with a bond length (~1.62 angstrom) consistent with the
  known TiO bond length.
- **Evidence**: Coordinate table labels both atoms as "Ti1" and "Ti2".
  Internal coordinates table says "rO=Ti = 1.620 angstrom" -- confirming
  one atom is O, not Ti. The coordinates themselves (0,0,0 and 0,0,1.6202) are
  correct; only the element labels are wrong.
- **Before fix**: xyz = "Ti 0 0 0\nTi 0 0 1.6202" (wrong element labels)
- **After fix**: xyz = "O 0 0 0\nTi 0 0 1.6202" (element labels corrected from
  internal coordinate description "rO=Ti" where atom 1=O, atom 2=Ti).
  best_available = "experimental" (geometry preserved, not discarded).
- **Parser fix**: When element validation fails, `_correct_xyz_elements_from_internal_coords()`
  extracts element identities from internal coordinate descriptions (e.g., "rO=Ti"
  encodes O and Ti as the two atoms). If the corrected element counts match the formula,
  the corrected xyz is used. Cross-validated against calculated geometry
  (CCSD(T)=FULL/aug-cc-pVQZ): Ti-O distance 1.610 angstrom, consistent with experimental
  1.620 angstrom.
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/OTi_12137201/molecule.json')); print(d['geometries']['experimental']['xyz'])"`
  Expected: `O 0.000000 0.000000 0.000000\nTi 0.000000 0.000000 1.620200`

---

## 4. Dummy Atom in Geometry (1 species)

### D3N -- Trideuteroamine / ND3 (CAS 13550497)
- **HTML**: geom3x.html (not archived on disk; evidence from parsed molecule.json
  and geom2x.html which lists the calculated method/basis)
- **Alternative hypothesis considered**: Could "X" be a real atom (e.g., xenon)?
  No — three confirmations: (1) "X" is not an element symbol in standard
  notation (Xenon is "Xe"); (2) "X" is universally used in computational
  chemistry to denote dummy/ghost atoms (orientation reference points for
  C3v symmetry specification); (3) the formula D3N has exactly 4 atoms
  (1 N + 3 D), and the 4 non-X atoms in the geometry match this formula.
- **Evidence**: Calculated geometry (CCSD(T)/cc-pVTZ) includes a 5th atom
  "X" at (0, 0, 1.1195). Formula D3N = 4 atoms. "X" is a computational
  reference point, not a real atom.
- **Before fix**: 5 atoms including X dummy
- **After fix**: 4 atoms (X filtered out)
- **Parser fix**: Skip atoms with symbol "X" during coordinate extraction.
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/D3N_13550497/molecule.json')); g=d['geometries']; calc=g.get('calculated',{}); print(sum(1 for m in calc.values() for b in m.values() if b.get('n_atoms')==4))"`
  Expected: nonzero (calculated geometry has 4 atoms)

---

## 5. Impossible Spin State (1 species)

### CNZn -- Zinc isocyanide (CAS 748101419)
- **HTML**: [data/molecules/CNZn_748101419/spin2x.html](../../data/molecules/CNZn_748101419/spin2x.html)
- **Alternative hypothesis considered**: Could this be a charged species (CNZn+ or
  CNZn-) where the "closed shell" label is correct but the charge=0 is wrong?
  Three independent pieces of evidence were examined to rule this out.

- **Evidence 1 — Electron parity**: Neutral CNZn has Zn(30) + C(6) + N(7) = 43
  electrons (odd). Odd-electron systems cannot be closed shell (singlet) -- this
  violates the Pauli exclusion principle. However, CNZn+ (42e) or CNZn- (44e) would
  be even-electron and *could* be singlet. So this alone is inconclusive.

- **Evidence 2 — CCCBDB spin page claim**: spin2x.html (line 656) states
  "ZnNC (Zinc isocyanide) is closed shell. S^2=0." The h1 tag (line 653)
  shows `<S^2> for  (Zinc isocyanide)` with no formula before the name --
  which is atypical and may indicate a generic label was applied.

- **Evidence 3 — ROHF energies in CCCBDB energy page** (decisive): The
  [energy2x.html](../../data/molecules/CNZn_748101419/energy2x.html) page
  (line 781) provides **ROHF** (Restricted Open-shell Hartree-Fock) energies
  alongside standard HF. ROHF is a method exclusively designed for open-shell
  systems with unpaired electrons. CCCBDB would not compute ROHF energies for
  a truly closed-shell molecule. The HF and ROHF energies differ slightly
  (e.g., HF/3-21G = -1860.945896, ROHF/3-21G = -1860.945491), confirming
  open-shell treatment. This means CCCBDB's own energy calculations internally
  treat CNZn as open-shell, directly contradicting the spin page's "closed
  shell" claim.

- **Ruling out wrong charge**: All URL parameters in archived HTML explicitly
  use `charge=0`. No CNZn+ or CNZn- entries exist in CCCBDB's species list
  (2,186 entries checked). Both CNZn isomers (ZnCN = CAS 73963981 and
  ZnNC = CAS 748101419) are listed as neutral.

- **Conclusion**: The charge=0 is correct. The "closed shell" label on the
  spin page is a genuine upstream CCCBDB data error -- the energy page's own
  ROHF data confirms the molecule is open-shell.

- **Before fix**: spin = {S_squared: 0.0, closed_shell: true, multiplicity: 1}
- **After fix**: spin = {} (attempted, no valid data -- parity violation detected)
- **Parser fix**: Check electron parity before trusting "closed shell" text.
  If n_electrons is odd, return None.
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/CNZn_748101419/molecule.json')); print(d['spin'])"`
  Expected: `{}` (empty dict, not a singlet)

---

## 6. Suspect Upstream Energy Data (1 species)

### CH4 -- Methane (CAS 74828)
- **HTML**: [data/molecules/CH4_74828/energy2x.html](../../data/molecules/CH4_74828/energy2x.html)
- **Evidence**: CCSD(T)=FULL/STO-3G = -39.806897 and CCSD(T)=FULL/6-31G\* = -39.806897
  are **identical**. 6-31G\* is a substantially larger basis set than STO-3G and should
  produce a more negative (lower) total energy at the CCSD(T)=FULL level. The identical
  values strongly suggest one was copied from the other in the upstream CCCBDB data.
- **Impact**: The 6-31G\* CCSD(T)=FULL energy for CH4 is unreliable. Other method/basis
  combinations for CH4 appear normal.
- **Status**: `upstream_suspect` -- faithfully reproduced from HTML. No parser fix applied.
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/CH4_74828/molecule.json')); e=d['energy']['methods']['CCSD(T)=FULL']; print(e.get('STO-3G'), e.get('6-31G*'))"`
  Expected: `-39.806897 -39.806897` (identical -- suspect)

---

## 7. Upstream Wrong Spin Data (1 species)

### Fe -- Iron atom (CAS 7439896)
- **HTML**: [data/molecules/Fe_7439896/spin2x.html](../../data/molecules/Fe_7439896/spin2x.html)
- **Evidence**: molecule.json reports `{S_squared: 0.0, closed_shell: true, multiplicity: 1}`.
  Atomic iron has 26 electrons with ground-state electron configuration [Ar] 3d6 4s2,
  giving term symbol 5D4 and multiplicity **5** (four unpaired 3d electrons). A singlet
  (multiplicity=1, all electrons paired) is physically incorrect for Fe.
- **Root cause confirmed**: The archived spin2x.html contains the line
  `Fe (Iron atom) is closed shell. S<SUP>2</SUP>=0.` — CCCBDB genuinely reports Fe as
  closed-shell. This is an **upstream data error**, not a session failure.
- **Why the parser did not catch this**: The electron parity check only flags odd-electron
  species claiming closed shell. Fe has 26 electrons (even), so a singlet is not
  *parity-forbidden* -- just physically wrong. The parser has no mechanism to validate
  multiplicity against known atomic ground states.
- **Impact**: Fe's multiplicity, s_squared, and closed_shell values in the CSV are
  unreliable. Researchers studying transition metal spin states should verify against
  CCCBDB directly.
- **Status**: `upstream_confirmed` -- archived HTML proves CCCBDB genuinely claims Fe is closed-shell.
- **Verification**: `grep "closed shell" data/molecules/Fe_7439896/spin2x.html`
  Expected: `Fe (Iron atom) is closed shell. S<SUP>2</SUP>=0.`

---

## 8. Upstream Species List Anomalies (2 entries)

### BF3H3N -- Suspected CAS Number Typo
- **Directories**: `BF3H3N_137019869` and `BF3H3N_13709869`
- **Evidence**: CAS `137019869` has name "Boron trifluoride monoammoniate" but no energy
  data and no geometry data (empty shell). CAS `13709869` has name "Amminetrifluoroboron"
  with full energy data (28 methods) and calculated geometry (CCSD(T)=FULL/aug-cc-pVDZ).
  The CAS numbers differ by one digit (`137019869` has an extra `1` inserted). Both
  entries exist in CCCBDB's authoritative 2,186-species list (`listallx.asp`).
- **Impact**: `BF3H3N_137019869` is an empty entry. Researchers iterating over all species
  will encounter it as a molecule with no data.
- **Status**: `upstream_suspect` -- likely a CAS number typo in CCCBDB's species list.
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/BF3H3N_137019869/molecule.json')); print(d.get('name'), 'energy' in d and d['energy'] is not None)"`
  Expected: `Boron trifluoride monoammoniate False`

### B12H12 -- Verified Not a Duplicate
- **Directories**: `B12H12--_12356137` and `B12H12--_84912390`
- **Initial suspicion**: Two CAS numbers for the same formula and charge appeared to
  be a duplicate entry.
- **Resolution via CAS lookup**: CAS 12356-13-7 = **lithium dodecaborate** (Li2B12H12),
  CAS 84912-39-0 = **cesium dodecaborate** (Cs2B12H12). These are different salts of
  the same B12H12^2- anion. CCCBDB strips the counterion (Li, Cs) and treats both as
  the bare anion. The different data availability (4 vs 25 energy methods) is explained
  by the different counterion environments in the reference calculations.
- **Status**: `resolved_not_anomaly` -- legitimate separate entries.
- **Lesson**: Always cross-check CAS numbers against external databases (PubChem,
  CAS Common Chemistry) before concluding "duplicate."

---

## 9. Data Completeness Warnings

These are not data errors but important limitations that researchers should be aware of.

### No Upstream Energy Data (27 species)
- **Finding**: 27 molecules have `energy2x.html` on disk but no parsed energy data
  in molecule.json. All 27 HTML files contain the correct molecule name in the `<h1>`
  tag (confirming the session was properly established), but the energy data tables
  are empty. CCCBDB genuinely has no energy data for these species.
- **Impact**: These species will have blank energy columns in the CSV export.
  This is an upstream data limitation, not a recoverable gap — re-fetching will
  not yield different results.
- **Species affected**: BF3H3N (CAS 137019869), C3H5F (CAS 1184607),
  C3H6Cl2 (CAS 78999), C3H6F2 (CAS 420451), C3H7+ (CAS 2143615),
  C4H5N (CAS 109751), C4H8O2S (CAS 126330), C4H9F (CAS 359013),
  C6H12 (CAS 18931710), C6H6 (CAS 5291907), CBr3I (CAS 14349805),
  CH4NO+ (CAS 9000924), CHCl2I (CAS 594047), CHI3 (CAS 75478),
  CHS2+ (CAS 87552861), Cl5Sb (CAS 7647189), ClH2N (CAS 110599903),
  ClSb (CAS 19952126), Co (CAS 7440484), Cr (CAS 7440473),
  and 7 others (see DATA_GAPS_AUDIT.md for full list).

---

## 10. Parser Data Transformation Decisions

These are deliberate choices made by the parser that transform data from the source HTML.
They are not errors, but researchers should be aware of them.

### 10a. Coordinate Precision Formatting
- **Transformation**: HTML coordinates have 4-5 decimal places (e.g., `0.1173`).
  The parser converts to float and reformats with 6 decimal places (e.g., `0.117300`).
- **Impact**: No precision loss (IEEE 754 double can exactly represent 4-5 digit
  decimals). However, trailing zeros are added and the original format is not preserved.
  Researchers comparing JSON to HTML will see cosmetic differences.

### 10b. DFT Exclusion from CSV Selections
- **Decision**: DFT energies (LSDA, B3LYP, M06-2X, etc.) are archived in the
  per-molecule JSON but excluded from `cheapest_energy` and `best_energy` in the CSV.
- **Rationale**: DFT is not systematically improvable -- unlike ab initio methods
  (HF -> MP2 -> CCSD -> CCSD(T)), DFT accuracy depends on the functional choice
  and cannot be extrapolated to the exact answer. This makes DFT unsuitable for
  a "cheapest/best" ranking paradigm. DFT can outperform low-level ab initio for
  some systems (especially transition metals), but the improvement is not systematic.
- **Impact**: 50 species with only DFT energies have blank energy columns in the CSV.

### 10c. Single Calculated Geometry Per Molecule
- **Decision**: Only the "best" non-DFT calculated geometry is stored in molecule.json
  (highest method rank, largest basis set). CCCBDB may provide dozens to hundreds of
  method/basis geometries per species.
- **Impact**: Researchers wanting a specific method/basis geometry cannot get it from
  the JSON. The archived `geom2x.html` lists all available methods, but individual
  `geom3x.html` pages (with actual coordinates) are only archived for the selected
  "best" geometry.

### 10d. Geometry Element Validation and Correction
- **Decision**: When the elements in a parsed experimental geometry do not match the
  molecular formula, the parser first attempts to correct element labels using internal
  coordinate descriptions (e.g., "rO=Ti" tells us which atom is O and which is Ti).
  If correction succeeds, the corrected geometry is stored. If correction fails (no
  internal coordinates, or ambiguous mapping), the geometry is discarded and
  `best_available` falls through to calculated geometry.
- **Impact**: Documented for OTi (Section 3) — successfully corrected from internal
  coords. If new upstream anomalies of this type exist and internal coordinates are
  available, they will be automatically corrected. If no internal coordinates exist,
  the geometry will be discarded.

### 10e. Silent Spin Parity Validation
- **Decision**: When CCCBDB reports "closed shell" for an odd-electron species, the
  parser discards the spin data and stores `spin={}`. No warning is logged.
- **Impact**: Documented for CNZn (Section 5). Correctly prevents physically impossible
  spin assignments, but the silent discard means the anomaly is only visible through
  the empty spin dict and the anomaly log.

---

## 11. Severe Spin Contamination (2 species)

**Category**: `severe_spin_contamination`
**Detection**: S² value deviates far from ideal value for the expected multiplicity
**Root cause**: These are legitimate computational results — multi-reference character
in the wavefunction causes severe spin contamination that persists even at high levels
of theory (CCSD(T)). This is not a CCCBDB data entry error.

### C3+ -- Carbon Trimer Cation (CAS 12075353)
- **HTML**: [data/molecules/C3+_12075353/spin2x.html](../../data/molecules/C3+_12075353/spin2x.html)
- **Electron count**: C(6)×3 - 1 = 17 electrons (odd → doublet, ideal S²=0.75)
- **Evidence**: spin2x.html CCSD(T) row shows S²=1.850 (STO-3G) and S²=1.823
  (6-31+G\*\*), with no after-annihilation values available. The parser selects
  1.823 from the highest-ranked basis set. The value 1.823 is far above the
  ideal doublet value of 0.75, indicating severe spin contamination.
- **molecule.json**: `spin = {S_squared: 1.823, closed_shell: false, multiplicity: 2}`
- **Status**: `upstream_confirmed` — legitimate severe spin contamination in CCCBDB
  data. Whitelisted in `test_s_squared_vs_multiplicity`.
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/C3+_12075353/molecule.json')); print(d['spin'])"`
  Expected: `{'S_squared': 1.823, 'closed_shell': False, 'multiplicity': 2}`

### CSi- -- Silicon Monocarbide Anion (CAS 409212)
- **HTML**: [data/molecules/CSi-_409212/spin2x.html](../../data/molecules/CSi-_409212/spin2x.html)
- **Electron count**: C(6) + Si(14) + 1 = 21 electrons (odd → doublet, ideal S²=0.75)
- **Evidence**: spin2x.html CCSD(T)=FULL row shows S²=1.751 in two basis sets,
  with no after-annihilation values available. The value 1.751 is far above the
  ideal doublet value of 0.75, indicating severe spin contamination.
- **molecule.json**: `spin = {S_squared: 1.751, closed_shell: false, multiplicity: 2}`
- **Status**: `upstream_confirmed` — legitimate severe spin contamination in CCCBDB
  data. Whitelisted in `test_s_squared_vs_multiplicity`.
- **Verification**: `python -c "import json; d=json.load(open('data/molecules/CSi-_409212/molecule.json')); print(d['spin'])"`
  Expected: `{'S_squared': 1.751, 'closed_shell': False, 'multiplicity': 2}`

---

## Summary Table

| Species | Category | Status | Resolution |
|---------|----------|--------|------------|
| C2H6S Ethanethiol (CAS 75081) | doubled_conformer | fixed_by_parser | Take first conformer |
| C3H5F Allyl Fluoride (CAS 818928) | doubled_conformer | fixed_by_parser | Take first conformer |
| C2H2F4 (CAS 359353) | doubled_conformer | fixed_by_parser | Take first conformer |
| H2S3 (CAS 13845233) | doubled_conformer | fixed_by_parser | Take first conformer |
| H3NS (CAS 14097008) | doubled_conformer | fixed_by_parser | Take first conformer |
| BrP (CAS 59727161) | truncated_cartesian | fixed_by_parser | Reconstruct from internal coords |
| C5H4O2 (CAS 980101) | truncated_cartesian | upstream_unrecoverable | Calculated fallback |
| Br (CAS 10097322) | wrong_element | fixed_by_parser | Validate against formula |
| Br+ (CAS 10097322) | wrong_element | fixed_by_parser | Validate against formula |
| Br- (CAS 10097322) | wrong_element | fixed_by_parser | Validate against formula |
| OTi (CAS 12137201) | wrong_element | fixed_by_parser | Correct labels from internal coords |
| D3N (CAS 13550497) | dummy_atom | fixed_by_parser | Filter X atoms |
| CNZn (CAS 748101419) | impossible_spin | fixed_by_parser | Parity check |
| CH4 (CAS 74828) | suspect_energy | upstream_suspect | Documented |
| Fe (CAS 7439896) | suspect_spin | upstream_confirmed | CCCBDB genuinely claims closed shell (verified via spin2x.html) |
| BF3H3N (CAS 137019869) | suspect_cas_typo | upstream_suspect | Documented |
| B12H12 (CAS 12356137/84912390) | verified_not_duplicate | resolved_not_anomaly | Li vs Cs salt — different compounds |
| C3+ (CAS 12075353) | severe_spin_contamination | upstream_confirmed | S²=1.823 at CCSD(T), ideal 0.75 |
| CSi- (CAS 409212) | severe_spin_contamination | upstream_confirmed | S²=1.751 at CCSD(T)=FULL, ideal 0.75 |
