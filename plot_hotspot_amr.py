#!/usr/bin/env python3
"""Plot AMRFinderPlus hits by defence-hotspot architecture."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, Normalize


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def natural_hotspot_key(value: str) -> tuple[str, int]:
    prefix, _, number = value.rpartition("_")
    return prefix, int(number) if number.isdigit() else 0


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_plots(results_dir: Path, annotation_dir: Path, output_prefix: Path) -> None:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    hotspots = read_tsv(results_dir / "hotspots.tsv")
    gene_map = read_tsv(annotation_dir / "hotspot_gene_map.tsv")
    amr_hits = read_tsv(annotation_dir / "amrfinder.tsv")

    gene_context = {row["gene_key"]: row for row in gene_map}
    loci_by_hotspot: dict[str, set[str]] = defaultdict(set)
    for row in hotspots:
        loci_by_hotspot[row["hotspot_id"]].add(row["locus_id"])

    hits_by_locus: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    unmatched: list[str] = []
    for hit in amr_hits:
        context = gene_context.get(hit["Protein id"])
        if not context:
            unmatched.append(hit["Protein id"])
            continue
        key = (context["hotspot_id"], context["locus_id"], context["genome"])
        hits_by_locus[key].add(hit["Element symbol"])

    hotspot_ids = sorted(loci_by_hotspot, key=natural_hotspot_key)
    positive_by_hotspot: dict[str, set[str]] = defaultdict(set)
    hit_count_by_hotspot: dict[str, int] = defaultdict(int)
    for (hotspot_id, locus_id, _genome), genes in hits_by_locus.items():
        positive_by_hotspot[hotspot_id].add(locus_id)
        hit_count_by_hotspot[hotspot_id] += len(genes)

    summary_rows: list[dict[str, object]] = []
    for hotspot_id in hotspot_ids:
        total = len(loci_by_hotspot[hotspot_id])
        positive = len(positive_by_hotspot[hotspot_id])
        summary_rows.append(
            {
                "hotspot_id": hotspot_id,
                "total_loci": total,
                "amr_positive_loci": positive,
                "amr_positive_percent": round(100 * positive / total, 4) if total else 0,
                "amr_gene_hits": hit_count_by_hotspot[hotspot_id],
            }
        )
    write_tsv(
        output_prefix.with_name(output_prefix.name + "_data.tsv"),
        summary_rows,
        ["hotspot_id", "total_loci", "amr_positive_loci", "amr_positive_percent", "amr_gene_hits"],
    )

    all_genes = {gene for values in hits_by_locus.values() for gene in values}
    preferred_order = ["vanR-B", "vanH-B", "vanB", "vanX-B", "vanY-B"]
    genes = [gene for gene in preferred_order if gene in all_genes]
    genes.extend(sorted(all_genes - set(genes)))
    # Keep the matrix continuous: all AMR genes occupy one horizontal panel.
    gene_chunks = [genes] if genes else []

    sns.set_theme(style="white", context="talk", font_scale=0.9)
    figure_height = 5.2 if not gene_chunks else 7.0
    figure_width = max(16.5, 5.0 + 0.55 * len(genes))
    fig = plt.figure(figsize=(figure_width, figure_height), constrained_layout=True)
    grid = fig.add_gridspec(
        1 + max(1, len(gene_chunks)), 1,
        height_ratios=[1.15] + [2.0] * max(1, len(gene_chunks)),
    )

    ax = fig.add_subplot(grid[0, 0])
    percentages = [float(row["amr_positive_percent"]) for row in summary_rows]
    y = np.arange(len(hotspot_ids))
    point_colours = ["#173F5F" if value else "#BCC5CC" for value in percentages]
    ax.hlines(y, 0, percentages, color="#173F5F", linewidth=2.2, zorder=1)
    ax.scatter(percentages, y, s=90, color=point_colours, edgecolor="white", linewidth=1.2, zorder=2)
    ax.set_title("AMR-positive loci", loc="left", weight="bold", pad=16)
    ax.set_xlabel("Loci with ≥1 AMR gene (%)")
    ax.set_ylabel("")
    ax.set_yticks(y, labels=hotspot_ids)
    ax.invert_yaxis()
    xmax = max(percentages + [0.1])
    ax.set_xlim(-0.04 * xmax, xmax * 1.48)
    ax.grid(False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#343A40")
    ax.tick_params(axis="y", length=0)
    for yi, row in zip(y, summary_rows):
        positive = int(row["amr_positive_loci"])
        total = int(row["total_loci"])
        value = float(row["amr_positive_percent"])
        ax.text(
            value + xmax * 0.045,
            yi,
            f"{positive}/{total}",
            ha="left",
            va="center",
            fontsize=10,
            color="#343A40",
        )

    if genes:
        locus_gene_counts: dict[tuple[str, str], int] = defaultdict(int)
        for (hotspot_id, _locus_id, _genome), locus_genes in hits_by_locus.items():
            for gene in locus_genes:
                locus_gene_counts[(hotspot_id, gene)] += 1
        matrix = np.array(
            [
                [100 * locus_gene_counts[(hotspot_id, gene)] / len(loci_by_hotspot[hotspot_id]) for gene in genes]
                for hotspot_id in hotspot_ids
            ]
        )
        cmap = LinearSegmentedColormap.from_list(
            "phorific_amr", ["#fffaf5", "#f4c49b", "#d77a55", "#26395f"]
        )
        vmax = max(float(matrix.max()), 0.1)
        heat_axes = []
        start = 0
        for panel, chunk in enumerate(gene_chunks, 1):
            ax2 = fig.add_subplot(grid[panel, 0])
            heat_axes.append(ax2)
            stop = start + len(chunk)
            panel_matrix = matrix[:, start:stop]
            annotations = np.full(panel_matrix.shape, "", dtype=object)
            for i, hotspot_id in enumerate(hotspot_ids):
                total = len(loci_by_hotspot[hotspot_id])
                for j, gene in enumerate(chunk):
                    count = locus_gene_counts[(hotspot_id, gene)]
                    if count:
                        annotations[i, j] = f"{count}/{total}\n{panel_matrix[i, j]:.1f}%"
            sns.heatmap(
                panel_matrix, ax=ax2, cmap=cmap, cbar=False,
                linewidths=0.7, linecolor="#FFFFFF",
                xticklabels=chunk, yticklabels=hotspot_ids,
                vmin=0, vmax=vmax, annot=annotations, fmt="",
                annot_kws={"fontsize": 8},
            )
            ax2.set_title("AMR gene frequency by hotspot", loc="left", weight="bold", pad=10)
            ax2.set_xlabel("AMRFinderPlus gene")
            ax2.set_ylabel("")
            ax2.tick_params(axis="x", rotation=42, labelsize=9, length=0)
            for label in ax2.get_xticklabels():
                label.set_horizontalalignment("right")
            ax2.tick_params(axis="y", rotation=0, labelsize=10, length=0)
            start = stop
        colourbar = fig.colorbar(
            plt.cm.ScalarMappable(norm=Normalize(0, vmax), cmap=cmap),
            ax=heat_axes, shrink=0.78, pad=0.015,
        )
        colourbar.set_label("Loci carrying gene (%)")
    else:
        ax2 = fig.add_subplot(grid[1, 0])
        ax2.text(0.5, 0.5, "No AMRFinderPlus hits mapped to hotspot loci", ha="center", va="center")
        ax2.set_axis_off()

    fig.suptitle("Antimicrobial resistance across defence hotspots", fontsize=20, weight="bold")
    fig.savefig(output_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    print(f"Mapped AMR hits: {sum(len(v) for v in hits_by_locus.values())}")
    print(f"AMR-positive loci: {len(hits_by_locus)}")
    print(f"Unmatched AMR protein IDs: {len(unmatched)}")
    print(output_prefix.with_suffix(".png").resolve())
    print(output_prefix.with_suffix(".pdf").resolve())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--annotation-dir", type=Path)
    parser.add_argument("--output-prefix", type=Path)
    args = parser.parse_args()
    results_dir = args.results_dir.resolve()
    annotation_dir = (args.annotation_dir or results_dir / "hotspot_annotations").resolve()
    output_prefix = (args.output_prefix or annotation_dir / "hotspot_amr_summary").resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    make_plots(results_dir, annotation_dir, output_prefix)


if __name__ == "__main__":
    main()
