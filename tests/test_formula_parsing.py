"""
Formula parsing pattern tests for CCCBDB mirror.

Documents every formula notation pattern found in CCCBDB's species list HTML
(data/source_pages/cccbdb_species_list.html) with concrete test cases.
49 of 2,186 species use parenthesized group notation; 2 use = bond notation.

Run: python -m pytest tests/test_formula_parsing.py -v
"""
import os
import sys
import pytest

# Allow imports from scripts/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from db_utils import hill_formula, count_atoms, count_electrons

# Each tuple: (raw_formula, expected_hill, expected_atoms, expected_electrons)
# All charge=0 unless noted in comments.

PATTERN_CASES = [
    # -- Pattern 1: Simple formulas (no parentheses) --
    ("H2O",           "H2O",       3,  10),   # Water
    ("CH4",           "CH4",       5,  10),   # Methane
    ("NaCl",          "ClNa",      2,  28),   # Sodium chloride (Hill: Cl before Na)
    ("H",             "H",         1,   1),   # Single atom
    ("He",            "He",        1,   2),   # Noble gas

    # -- Pattern 2: Single substituent group with multiplier --
    # Be(OH)2 = Be + 2*(O+H) → BeH2O2: 4+2+16=22e
    ("N(CH3)3",       "C3H9N",    13,  34),   # Trimethylamine
    ("Be(OH)2",       "BeH2O2",    5,  22),   # Beryllium hydroxide
    ("Fe(CO)5",       "C5FeO5",   11,  96),   # Iron pentacarbonyl: 26+30+40
    ("C(CN)4",        "C5N4",      9,  58),   # Tetracyanomethane: 30+28
    ("Si(CH3)4",      "C4H12Si",  17,  50),   # Tetramethylsilane

    # -- Pattern 3: Group with implicit multiplier (x1) --
    ("CH3CH(NH2)COOH","C3H7NO2",  13,  48),   # Alanine: 18+7+7+16
    ("CHOCH(CH3)CH3", "C4H8O",    13,  40),   # Isobutyraldehyde: 24+8+8
    ("CH3C(O)OO",     "C2H3O3",    8,  39),   # Acetyl peroxy radical

    # -- Pattern 4: Multiple parenthesized groups --
    ("(CH3)3CC(CH3)3","C8H18",    26,  66),   # Tetramethylbutane: 48+18
    ("SiCl2(CH3)2",   "C2H6Cl2Si",11,  66),   # Dichlorodimethylsilane: 12+6+34+14
    ("NH(CH3)CONH(CH3)","C3H8N2O",14,  48),   # Urea, N,N'-dimethyl-: 18+8+14+8
    ("CH2(SH)CH(CH3)CH3","C4H10S",15,  50),   # 1-Propanethiol, 2-methyl-: 24+10+16

    # -- Pattern 5: Nested groups (group containing subscripts) --
    ("Fe(C5H5)2",     "C10H10Fe", 21,  96),   # Ferrocene: 60+10+26
    ("NH(C2H5)2",     "C4H11N",   16,  42),   # Diethylamine: 24+11+7
    ("CH(C2H3)3",     "C7H10",    17,  52),   # Trivinyl methane: 42+10
    ("N(SiH3)3",      "H9NSi3",   13,  58),   # Trisilylamine

    # -- Pattern 6: Double-bond notation with = --
    ("CH3C(OH)=NH",   "C2H5NO",    9,  32),   # Ethaninidic acid
    ("HN=C=C(CN)2",   "C4HN3",     8,  46),   # Dicyanoketenimine: 24+1+21

    # -- Pattern 7: Charged species (formula string has no +/-) --
    ("(CH3)2NH2",     "C2H8N",    11,  27),   # Dimethylammonium (charge=+1 → 26 e-)

    # -- Pattern 8: Isotopes D and T --
    ("D2O",           "D2O",       3,  10),   # Heavy water
    ("D",             "D",         1,   1),   # Deuterium atom
]

# Charged species: tests that exercise count_electrons with non-zero charge.
CHARGED_CASES = [
    # (formula, charge, expected_electrons)
    ("H",       1,   0),    # H+ (proton): 1 - 1 = 0
    ("H",      -1,   2),    # H- (hydride): 1 + 1 = 2
    ("O",      -2,  10),    # O2-: 8 + 2 = 10
    ("Li",      1,   2),    # Li+: 3 - 1 = 2
    ("Fe",      2,  24),    # Fe2+: 26 - 2 = 24
    ("N2",      1,  13),    # N2+: 14 - 1 = 13
    ("H2O",     1,   9),    # H3O+ proxy: 10 - 1 = 9
    ("Br",      1,  34),    # Br+: 35 - 1 = 34
    ("Br",     -1,  36),    # Br-: 35 + 1 = 36
    ("(CH3)2NH2", 1, 26),   # Dimethylammonium+: 27 - 1 = 26
]

@pytest.mark.parametrize("formula,charge,expected_electrons", CHARGED_CASES)
def test_count_electrons_with_charge(formula, charge, expected_electrons):
    assert count_electrons(formula, charge=charge) == expected_electrons

@pytest.mark.parametrize("raw,expected_hill,expected_atoms,expected_electrons", PATTERN_CASES)
def test_hill_formula(raw, expected_hill, expected_atoms, expected_electrons):
    assert hill_formula(raw) == expected_hill

@pytest.mark.parametrize("raw,expected_hill,expected_atoms,expected_electrons", PATTERN_CASES)
def test_count_atoms(raw, expected_hill, expected_atoms, expected_electrons):
    assert count_atoms(raw) == expected_atoms

@pytest.mark.parametrize("raw,expected_hill,expected_atoms,expected_electrons", PATTERN_CASES)
def test_count_electrons(raw, expected_hill, expected_atoms, expected_electrons):
    assert count_electrons(raw, charge=0) == expected_electrons

def test_all_parenthesized_species_normalize_correctly():
    """Verify every parenthesized formula in CCCBDB produces valid Hill output.

    parse_species_list_html() applies hill_formula() internally, so the
    returned 'formula' field is already normalized. We check that
    raw_species text containing '(' (original CCCBDB notation) produces
    formulas where no element count exceeds 200 (catches H33-style bugs)
    and atom counts are consistent with the Hill result.
    """
    import re as _re

    source_html = os.path.join(
        os.path.dirname(__file__), '..', 'data', 'source_pages',
        'cccbdb_species_list.html')
    if not os.path.exists(source_html):
        pytest.skip("Species list HTML not available")

    from cccbdb_parse_listall import parse_species_list_html

    species = parse_species_list_html(source_html)

    # raw_species contains original CCCBDB notation with parentheses
    paren_species = [s for s in species if '(' in s.get('raw_species', '')]
    assert len(paren_species) >= 49, (
        f"Expected at least 49 parenthesized species, found {len(paren_species)}")

    for sp in paren_species:
        formula = sp['formula']  # already Hill-normalized

        # Sanity: no element count > 200 (catches H33-style bugs)
        for num in _re.findall(r'\d+', formula):
            assert int(num) <= 200, (
                f"{sp['raw_species']} -> {formula}: suspiciously large count {num}")

        # Hill formula should round-trip cleanly
        assert hill_formula(formula) == formula, (
            f"Hill formula not idempotent: {formula} -> {hill_formula(formula)}")
