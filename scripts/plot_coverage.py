#!/usr/bin/env python3
"""
plot_coverage.py

Generates a horizontal stacked bar chart showing data coverage across 2,186
species (geometry, spin, energy). Reads from the CSV and per-molecule JSON.

Usage:
    python scripts/plot_coverage.py
    python scripts/plot_coverage.py --output docs/data-coverage.png --dpi 180
"""

import argparse
import csv
import json
import os
import sys

# ─── Path setup (same pattern as export_csv.py) ───

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
REPO_ROOT = os.path.join(SCRIPT_DIR, "..")
DATA_DIR = os.path.join(REPO_ROOT, "data")
MOLECULES_DIR = os.path.join(DATA_DIR, "molecules")
CSV_PATH = os.path.join(REPO_ROOT, "cccbdb-selected-energy.csv")

# Ab initio method prefixes — anything else is DFT
AB_INITIO_METHODS = {
    "HF", "ROHF",
    "MP2", "MP2=FULL", "PMP2", "ROMP2",
    "MP3", "MP3=FULL",
    "MP4", "MP4=FULL",
    "CID", "CID=FULL", "CISD",
    "CCD",
    "CCSD", "CCSD=FULL",
    "CCSD(T)", "CCSD(T)=FULL",
    "QCISD", "QCISD=FULL",
    "QCISD(T)", "QCISD(T)=FULL",
    "QCISD(TQ)", "QCISD(TQ)=FULL",
}

TOTAL_SPECIES = 2186

def compute_counts():
    """Compute coverage partition counts from CSV + JSON.

    Returns dict with geometry, spin, and energy breakdowns.
    """
    # ── Read CSV ──
    with open(CSV_PATH, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # ── Geometry ──
    geom_exp_and_calc = 0
    geom_exp_only = 0
    geom_calculated = 0
    geom_single_atom = 0
    geom_missing = 0
    for row in rows:
        src = row["geometry_source"].strip()
        if src == "experimental":
            # Check molecule.json for calculated geometry overlap
            formula = row["formula"]
            casno = row["casno"]
            charge = int(row["charge"])
            if charge > 0:
                suffix = "+" * charge
            elif charge < 0:
                suffix = "-" * abs(charge)
            else:
                suffix = ""
            dirname = f"{formula}{suffix}_{casno}"
            json_path = os.path.join(MOLECULES_DIR, dirname, "molecule.json")
            has_calc = False
            if os.path.isfile(json_path):
                with open(json_path, "r", encoding="utf-8") as f:
                    mol = json.load(f)
                calc = mol.get("geometries", {}).get("calculated", {})
                if calc:
                    has_calc = True
            if has_calc:
                geom_exp_and_calc += 1
            else:
                geom_exp_only += 1
        elif src == "calculated":
            geom_calculated += 1
        elif src == "single atom":
            geom_single_atom += 1
        else:
            geom_missing += 1

    # ── Spin ──
    spin_has = sum(1 for r in rows if r["multiplicity"].strip())
    spin_missing = len(rows) - spin_has

    # ── Energy: ab initio vs DFT-only vs missing ──
    # Species with cheapest_energy in CSV have ab initio data
    has_ab_initio = set()
    no_ab_initio = set()
    for row in rows:
        sid = row["species_id"]
        if row["cheapest_energy"].strip():
            has_ab_initio.add(sid)
        else:
            no_ab_initio.add(sid)

    # For species without ab initio energy, check molecule.json for DFT
    dft_only = 0
    energy_missing = 0
    for row in rows:
        sid = row["species_id"]
        if sid in has_ab_initio:
            continue

        # Find molecule directory
        formula = row["formula"]
        casno = row["casno"]
        charge = int(row["charge"])
        # Inline species_id builder (avoids dependency on fetcher scripts)
        if charge > 0:
            suffix = "+" * charge
        elif charge < 0:
            suffix = "-" * abs(charge)
        else:
            suffix = ""
        dirname = f"{formula}{suffix}_{casno}"
        json_path = os.path.join(MOLECULES_DIR, dirname, "molecule.json")

        has_any_energy = False
        if os.path.isfile(json_path):
            with open(json_path, "r", encoding="utf-8") as f:
                mol = json.load(f)
            energy = mol.get("energy")
            if energy and isinstance(energy, dict):
                methods = energy.get("methods", {})
                if methods:
                    has_any_energy = True

        if has_any_energy:
            dft_only += 1
        else:
            energy_missing += 1

    energy_ab_initio = len(has_ab_initio)

    # ── Validate partitions ──
    geom_total = geom_exp_and_calc + geom_exp_only + geom_calculated + geom_single_atom + geom_missing
    spin_total = spin_has + spin_missing
    energy_total = energy_ab_initio + dft_only + energy_missing

    assert geom_total == TOTAL_SPECIES, \
        f"Geometry partition {geom_total} != {TOTAL_SPECIES}"
    assert spin_total == TOTAL_SPECIES, \
        f"Spin partition {spin_total} != {TOTAL_SPECIES}"
    assert energy_total == TOTAL_SPECIES, \
        f"Energy partition {energy_total} != {TOTAL_SPECIES}"

    return {
        "geometry": [
            ("Exp + Calc", geom_exp_and_calc),
            ("Exp only", geom_exp_only),
            ("Calc only", geom_calculated),
            ("Single atoms", geom_single_atom),
            ("Missing", geom_missing),
        ],
        "spin": [
            ("Has data", spin_has),
            ("Missing", spin_missing),
        ],
        "energy": [
            ("Ab initio", energy_ab_initio),
            ("DFT-only", dft_only),
            ("Missing", energy_missing),
        ],
    }

def plot_coverage(counts, output_path, dpi):
    """Generate horizontal stacked bar chart."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Okabe-Ito colorblind-friendly palette
    BLUE = "#0072B2"
    SKY = "#56B4E9"
    GREEN = "#009E73"
    VERMILLION = "#D55E00"

    AMBER = "#E69F00"

    color_map = {
        "geometry": [BLUE, SKY, GREEN, AMBER, VERMILLION],
        "spin": [BLUE, VERMILLION],
        "energy": [BLUE, SKY, VERMILLION],
    }

    bar_labels = {
        "geometry": "Geometry",
        "spin": "Spin",
        "energy": "Energy",
    }

    # Per-bar annotation placement: list of "above"/"below" per small segment.
    # Geometry: "Exp only" below, "Single atoms" above, "Missing" above
    # Energy: "DFT-only" above, "Missing" below
    annotation_side = {
        "geometry": ["below", "above", "above"],
        "spin": ["above"],
        "energy": ["above", "below"],
    }

    fig, ax = plt.subplots(figsize=(9, 3.6))

    bar_order = ["geometry", "spin", "energy"]
    y_positions = [2.5, 1.25, 0]
    bar_height = 0.55

    for bar_key, y_pos in zip(bar_order, y_positions):
        segments = counts[bar_key]
        colors = color_map[bar_key]
        left = 0
        small_annotations = []  # collect (cx, label_text, side)
        sides = annotation_side.get(bar_key, [])

        for i, (label, count) in enumerate(segments):
            pct = count / TOTAL_SPECIES * 100
            width = pct
            ax.barh(y_pos, width, left=left, height=bar_height,
                    color=colors[i], edgecolor="white", linewidth=0.5)

            # Label inside bar if segment > 8%, otherwise defer annotation
            if pct > 8:
                cx = left + width / 2
                text = f"{label}\n{count:,} ({pct:.1f}%)"
                text_color = "white" if colors[i] in (BLUE, VERMILLION, GREEN) else "#333333"
                ax.text(cx, y_pos, text, ha="center", va="center",
                        fontsize=7.5, fontweight="medium", color=text_color,
                        linespacing=1.2)
            else:
                cx = left + width / 2
                side_idx = len(small_annotations)
                side = sides[side_idx] if side_idx < len(sides) else "below"
                small_annotations.append(
                    (cx, f"{label}: {count} ({pct:.1f}%)", side))

            left += width

        # Place small annotations above or below the bar
        if small_annotations:
            above_idx = 0
            below_idx = 0
            for anchor_x, text, side in small_annotations:
                text_ha = "right" if anchor_x > 90 else "center"
                text_x = anchor_x + 1.5 if anchor_x > 90 else anchor_x
                if side == "above":
                    y_off = 0.30 + 0.27 * above_idx
                    above_idx += 1
                    ax.annotate(
                        text,
                        xy=(anchor_x, y_pos + bar_height / 2),
                        xytext=(text_x, y_pos + bar_height / 2 + y_off),
                        fontsize=6.5, ha=text_ha, va="bottom", color="#555555",
                        arrowprops=dict(arrowstyle="-", color="#999999", lw=0.7),
                    )
                else:
                    y_off = 0.30 + 0.27 * below_idx
                    below_idx += 1
                    ax.annotate(
                        text,
                        xy=(anchor_x, y_pos - bar_height / 2),
                        xytext=(text_x, y_pos - bar_height / 2 - y_off),
                        fontsize=6.5, ha=text_ha, va="top", color="#555555",
                        arrowprops=dict(arrowstyle="-", color="#999999", lw=0.7),
                    )

    # Y-axis labels
    ax.set_yticks(y_positions)
    ax.set_yticklabels([bar_labels[k] for k in bar_order],
                       fontsize=10, fontweight="bold")

    # X-axis
    ax.set_xlim(0, 102)
    ax.set_xlabel("% of 2,186 species", fontsize=9)
    ax.set_xticks([0, 20, 40, 60, 80])
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.0f}%"))

    # Clean up spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(left=False)

    # Title
    ax.set_title("Data Coverage: CCCBDB Mirror (2,186 species)",
                 fontsize=11, fontweight="bold", pad=10)

    plt.tight_layout()

    # Ensure output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"Saved: {output_path} ({dpi} DPI)")

def main():
    parser = argparse.ArgumentParser(
        description="Generate data coverage chart for README")
    parser.add_argument("--output", default=os.path.join(REPO_ROOT, "docs", "data-coverage.png"),
                        help="Output image path (default: docs/data-coverage.png)")
    parser.add_argument("--dpi", type=int, default=180,
                        help="Image DPI (default: 180)")
    args = parser.parse_args()

    counts = compute_counts()

    # Print summary
    for key in ["geometry", "spin", "energy"]:
        parts = ", ".join(f"{label}: {count}" for label, count in counts[key])
        print(f"  {key.capitalize()}: {parts}")

    plot_coverage(counts, args.output, args.dpi)

if __name__ == "__main__":
    main()
