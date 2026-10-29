#!/usr/bin/env python3
"""Summarise defence-system prevalence across standard and edge MGE hotspots."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


EMPTY = "no defence detected"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split_systems(value: str) -> list[str]:
    return sorted({part.strip() for part in value.split(";")
                   if part.strip() and part.strip().lower() != EMPTY})


def hotspot_key(value: str):
    parts = re.split(r"(\d+)$", value)
    return parts[0], int(parts[1]) if len(parts) > 1 and parts[1] else -1


def write_tsv(path: Path, header: list[str], rows) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot defence-system prevalence for MGE hotspots")
    parser.add_argument("--results-dir", type=Path, required=True,
                        help="merged MGE annotation-input directory containing hotspots.tsv")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--title", required=True)
    args = parser.parse_args()
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(args.results_dir.resolve() / "hotspots.tsv")

    total = Counter()
    positive = Counter()
    empty = Counter()
    multi = Counter()
    counts = Counter()
    types = defaultdict(set)
    for row in rows:
        hotspot = row["hotspot_id"]
        systems = split_systems(row.get("defence_system", ""))
        total[hotspot] += 1
        if systems:
            positive[hotspot] += 1
        else:
            empty[hotspot] += 1
        if len(systems) > 1:
            multi[hotspot] += 1
        for system in systems:
            counts[(hotspot, system)] += 1
            types[hotspot].add(system)

    hotspots = sorted(total, key=hotspot_key)
    systems = sorted({system for _, system in counts},
                     key=lambda system: (-sum(counts[(hotspot, system)] for hotspot in hotspots), system))
    summary = []
    distribution = []
    for hotspot in hotspots:
        occurrences = sum(counts[(hotspot, system)] for system in types[hotspot])
        summary.append([
            hotspot, total[hotspot], positive[hotspot], empty[hotspot],
            f"{positive[hotspot] / total[hotspot]:.6f}", len(types[hotspot]),
            occurrences, multi[hotspot],
        ])
        for system in systems:
            count = counts[(hotspot, system)]
            if not count:
                continue
            distribution.append([
                hotspot, system, count, total[hotspot], positive[hotspot],
                f"{count / total[hotspot]:.6f}",
                f"{count / positive[hotspot]:.6f}" if positive[hotspot] else "0.000000",
            ])
    write_tsv(out / "hotspot_defence_summary.tsv", [
        "hotspot_id", "total_loci", "defence_positive_loci", "boundary_only_loci",
        "defence_positive_fraction", "distinct_system_types", "system_occurrences", "multi_system_loci",
    ], summary)
    write_tsv(out / "defence_distribution_by_hotspot.tsv", [
        "hotspot_id", "system_type", "loci_with_system", "total_hotspot_loci",
        "defence_positive_loci", "fraction_of_all_loci", "fraction_of_defence_loci",
    ], distribution)

    matrix = np.array([[100 * counts[(hotspot, system)] / total[hotspot] for system in systems]
                       for hotspot in hotspots])
    sns.set_theme(style="white", context="paper", font_scale=1.05)
    cmap = LinearSegmentedColormap.from_list(
        "phorific_prevalence", ["#fffaf5", "#f4c49b", "#d77a55", "#26395f"]
    )
    # Keep the full defence matrix continuous in one horizontal panel.
    width = max(15.0, 5.0 + 0.52 * len(systems))
    height = max(5.5, 3.0 + 0.55 * len(hotspots))
    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)
    ylabels = [f"{h}  (n={total[h]})" for h in hotspots]
    annotations = np.array([
        [f"{value:.1f}%" if value > 0 else "" for value in row]
        for row in matrix
    ])
    sns.heatmap(
        matrix, ax=ax, cmap=cmap, vmin=0, vmax=100,
        annot=annotations, fmt="", annot_kws={"fontsize": 7.5},
        linewidths=0.6, linecolor="white",
        cbar_kws={"label": "Hotspot member loci containing system (%)", "shrink": 0.76},
        xticklabels=systems, yticklabels=ylabels,
    )
    ax.set_title(args.title, loc="left", fontsize=15, fontweight="bold", pad=14)
    ax.set_xlabel("Defence system")
    ax.set_ylabel("Hotspot architecture (member loci)")
    ax.tick_params(axis="x", rotation=42, labelsize=9)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    ax.tick_params(axis="y", rotation=0, labelsize=9, length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.savefig(out / "defence_system_prevalence_heatmap.pdf", bbox_inches="tight")
    fig.savefig(out / "defence_system_prevalence_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Summarised {len(rows)} loci across {len(hotspots)} hotspots and {len(systems)} systems.")


if __name__ == "__main__":
    main()
