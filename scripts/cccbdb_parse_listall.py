#!/usr/bin/env python3
"""
cccbdb_parse_listall.py

Purpose:
    Parse the local CCCBDB "list of species" HTML (listallx.asp dump)
    into a structured species list.

    The core function parse_species_list_html() is imported by
    cccbdb_mirror.py for bootstrapping the molecule directory structure.

Input:  data/source_pages/cccbdb_species_list.html
Output: data/cas_lookup.csv (when run as CLI)
"""

import csv
import os
import re
from urllib.parse import urlparse, parse_qs
from bs4 import BeautifulSoup, NavigableString, Tag

import db_utils

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_HTML = os.path.join(SCRIPT_DIR, "..", "data", "source_pages",
                            "cccbdb_species_list.html")
DEFAULT_OUTPUT = os.path.join(SCRIPT_DIR, "..", "data", "cas_lookup.csv")

FIELDNAMES = ["formula", "name", "casno", "charge", "raw_species"]

def parse_species_cell(cell):
    """Extract raw formula (without charge) and charge from a species <td>.

    Species cells contain element symbols as text, stoichiometry in <sub>,
    and charge signs in <sup>.  Example: BH<sub>3</sub><sup>-</sup>
    """
    parts = []
    charge = 0
    for c in cell.contents:
        if isinstance(c, Tag):
            t = c.get_text(strip=True)
            if c.name == "sub":
                parts.append(t)
            elif c.name == "sup":
                # Charge sign(s)
                if t == "+":
                    charge = 1
                elif t == "-":
                    charge = -1
                elif t == "++":
                    charge = 2
                elif t == "--":
                    charge = -2
                elif t.endswith("+"):
                    try:
                        charge = int(t[:-1])
                    except ValueError:
                        charge = 1
                elif t.endswith("-"):
                    try:
                        charge = -int(t[:-1])
                    except ValueError:
                        charge = -1
            else:
                parts.append(t)
        elif isinstance(c, NavigableString):
            t = str(c).strip()
            if t:
                parts.append(t)
    raw_formula = "".join(parts)
    return raw_formula, charge

def parse_name_cell(cell):
    """Extract molecule name and CAS/charge from the <a> tag in a name <td>.

    Name cells look like:
        <a href="alldata2x.asp?casno=12385136&charge=0">Hydrogen atom</a>

    Returns (name, casno, charge) or (None, None, None) if empty.
    """
    a_tag = cell.find("a")
    if not a_tag:
        return None, None, None

    name = a_tag.get_text(strip=True)
    href = a_tag.get("href", "")

    # Parse casno and charge from URL query params
    # Handle both relative and absolute URLs
    parsed = urlparse(href)
    params = parse_qs(parsed.query)

    casno = params.get("casno", [None])[0]
    charge_str = params.get("charge", ["0"])[0]
    try:
        charge = int(charge_str)
    except ValueError:
        charge = 0

    return name, casno, charge

def _is_data_table(table):
    """Check if a table is a data table (has Species/Name headers)."""
    th_tags = table.find_all("th")
    return any("Species" in th.get_text() for th in th_tags)

def parse_species_list_html(html_path=None):
    """Parse CCCBDB species list HTML into a list of species dicts.

    Args:
        html_path: Path to cccbdb_species_list.html. Defaults to
                   data/source_pages/cccbdb_species_list.html.

    Returns:
        List of dicts, each with keys: formula, name, casno, charge, raw_species.
        Formulas are Hill-normalized. Deduplicated by (formula, charge, casno).
    """
    if html_path is None:
        html_path = DEFAULT_HTML

    if not os.path.exists(html_path):
        raise FileNotFoundError(
            f"Species list HTML not found: {html_path}\n"
            "Expected data/source_pages/cccbdb_species_list.html"
        )

    with open(html_path, "r", encoding="utf-8") as f:
        soup = BeautifulSoup(f, "html.parser")

    rows_out = []
    seen = set()  # Track (formula, charge, casno) to avoid exact duplicates

    # Find all data tables (those with Species/Name headers)
    for table in soup.find_all("table", attrs={"border": "1"}):
        if not _is_data_table(table):
            continue

        for tr in table.find_all("tr"):
            # Skip header rows
            if tr.find("th"):
                continue

            tds = tr.find_all("td")
            if len(tds) < 2:
                continue

            # Process pairs: (species, name) at indices (0,1), (2,3), (4,5)
            for pair_idx in range(0, len(tds) - 1, 2):
                species_td = tds[pair_idx]
                name_td = tds[pair_idx + 1]

                # Skip empty pairs
                if not species_td.get_text(strip=True) and not name_td.find("a"):
                    continue

                raw_formula, species_charge = parse_species_cell(species_td)
                name, casno, url_charge = parse_name_cell(name_td)

                if not name or not raw_formula:
                    continue

                # Use charge from URL (authoritative), fall back to species cell
                charge = url_charge if url_charge is not None else species_charge

                # Hill-normalize the formula
                formula = db_utils.hill_formula(raw_formula)

                # Build raw_species string for debugging
                raw_species = species_td.get_text(strip=True)

                # Dedup key
                key = (formula, charge, casno)
                if key in seen:
                    continue
                seen.add(key)

                rows_out.append({
                    "formula": formula,
                    "name": name,
                    "casno": casno,
                    "charge": charge,
                    "raw_species": raw_species,
                })

    return rows_out

def main():
    species = parse_species_list_html()

    # Write output
    os.makedirs(os.path.dirname(DEFAULT_OUTPUT), exist_ok=True)
    with open(DEFAULT_OUTPUT, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(species)

    print(f"Parsed {len(species)} entries → {DEFAULT_OUTPUT}")

    # Summary stats
    charges = {}
    for r in species:
        c = r["charge"]
        charges[c] = charges.get(c, 0) + 1
    for c in sorted(charges):
        label = "neutral" if c == 0 else f"charge={c:+d}"
        print(f"  {label}: {charges[c]}")

if __name__ == "__main__":
    main()
