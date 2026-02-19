"""
db_utils.py — Shared utilities for CCCBDB mirror scripts.

Provides: ELEMENT_Z, hill_formula, count_atoms, count_electrons.
"""

import re

# --- Constants ---

ELEMENT_Z = {
    'H': 1, 'D': 1, 'T': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
    'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
    'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40,
    'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60,
    'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70,
    'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80,
    'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90,
    'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100,
    'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109,
    'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
}

# --- Functions ---

def _parse_formula(formula):
    """Parse formula with parenthesized group expansion.

    Handles: N(CH3)3 → {N:1, C:3, H:9}, Fe(C5H5)2 → {Fe:1, C:10, H:10},
    CH3C(OH)=NH → {C:2, H:5, N:1, O:1}. Uses a stack for nested groups.
    Strips charge signs (+/-) and bond notation (=) before parsing.
    """
    clean = formula.replace('+', '').replace('-', '').replace('=', '')
    stack = [{}]
    i = 0
    while i < len(clean):
        if clean[i] == '(':
            stack.append({})
            i += 1
        elif clean[i] == ')':
            i += 1
            m = re.match(r'(\d+)', clean[i:])
            mult = int(m.group(1)) if m else 1
            if m:
                i += len(m.group(1))
            if len(stack) < 2:
                # Unbalanced closing paren — skip it
                continue
            group = stack.pop()
            for el, cnt in group.items():
                stack[-1][el] = stack[-1].get(el, 0) + cnt * mult
        else:
            m = re.match(r'([A-Z][a-z]?)(\d*)', clean[i:])
            if m and m.group(1):
                el = m.group(1)
                cnt = int(m.group(2)) if m.group(2) else 1
                stack[-1][el] = stack[-1].get(el, 0) + cnt
                i += len(m.group(0))
            else:
                i += 1
    return stack[0]

def hill_formula(formula):
    """Normalize a formula string to Hill convention.

    Organic (contains C): C, H, then remaining elements alphabetical.
    Inorganic (no C): all elements alphabetical.
    Handles parenthesized groups: N(CH3)3 → C3H9N.
    """
    counts = _parse_formula(formula)
    result = ""
    if 'C' in counts:
        result += f"C{counts['C'] if counts['C'] > 1 else ''}"
        del counts['C']
        if 'H' in counts:
            result += f"H{counts['H'] if counts['H'] > 1 else ''}"
            del counts['H']
    for el in sorted(counts.keys()):
        result += f"{el}{counts[el] if counts[el] > 1 else ''}"
    return result

def count_atoms(formula):
    """Count total atoms in formula, expanding parenthesized groups."""
    counts = _parse_formula(formula)
    return sum(counts.values())

def count_electrons(formula, charge):
    """Count total electrons from formula elements minus charge."""
    counts = _parse_formula(formula)
    total = 0
    for el, cnt in counts.items():
        if el in ELEMENT_Z:
            total += cnt * ELEMENT_Z[el]
        else:
            return None
    return total - int(charge)
