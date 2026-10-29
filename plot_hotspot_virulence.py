#!/usr/bin/env python3
"""Plot selected protein-VFDB matches by hotspot and category."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


CATEGORIES = ["adhesin", "exotoxin", "endotoxin", "T7SS"]
COLORS = {
    "adhesin": "#F3A261",
    "exotoxin": "#C96A4A",
    "endotoxin": "#E6B566",
    "T7SS": "#5B6F91",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = ["hotspot_id", "category", "positive_genes", "positive_loci", "total_loci", "positive_percent"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_virulence_plot(results_dir: Path, annotation_dir: Path) -> None:
    hits = read_tsv(annotation_dir / "diamond_vfdb_selected.tsv")
    hotspots = read_tsv(results_dir / "hotspots.tsv")
    total_loci: dict[str, set[str]] = defaultdict(set)
    for row in hotspots:
        total_loci[row["hotspot_id"]].add(row["locus_id"])
    hotspot_ids = sorted(total_loci)

    loci: dict[tuple[str, str], set[str]] = defaultdict(set)
    genes: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in hits:
        key = (row["hotspot_id"], row["category"])
        loci[key].add(row["locus_id"])
        genes[key].add(row["gene_key"])

    rows: list[dict[str, object]] = []
    for hotspot in hotspot_ids:
        for category in CATEGORIES:
            count = len(loci[(hotspot, category)])
            total = len(total_loci[hotspot])
            rows.append(
                {
                    "hotspot_id": hotspot,
                    "category": category,
                    "positive_genes": len(genes[(hotspot, category)]),
                    "positive_loci": count,
                    "total_loci": total,
                    "positive_percent": round(100 * count / total, 4) if total else 0,
                }
            )
    write_tsv(annotation_dir / "hotspot_virulence_summary_data.tsv", rows)

    fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharey=True, constrained_layout=True)
    y = np.arange(len(CATEGORIES))
    for ax, hotspot in zip(axes.flat, hotspot_ids):
        values = [len(loci[(hotspot, category)]) for category in CATEGORIES]
        colors = [COLORS[category] for category in CATEGORIES]
        ax.barh(y, values, color=colors, height=0.58)
        ax.set_title(hotspot.replace("hotspot_", "Hotspot "), loc="left", weight="bold")
        ax.set_yticks(y, CATEGORIES)
        ax.invert_yaxis()
        ax.grid(False)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.set_xlabel("Loci with match")
        xmax = max(values) if max(values) else 1
        ax.set_xlim(0, max(1.2, xmax * 1.3))
        for yi, value in zip(y, values):
            if value:
                ax.text(value + xmax * 0.04, yi, str(value), va="center", fontsize=10)
    for ax in axes.flat[len(hotspot_ids):]:
        ax.set_visible(False)
    fig.suptitle("Selected VFDB protein matches by hotspot", fontsize=18, weight="bold")
    fig.savefig(annotation_dir / "hotspot_virulence_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(annotation_dir / "hotspot_virulence_summary.pdf", bbox_inches="tight")
    plt.close(fig)
