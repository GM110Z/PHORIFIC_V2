#!/usr/bin/env python3
"""Summarise fallback Pfam annotations across hotspot architectures."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from hotspot_annotation_pipeline import parse_amr, parse_pfam, read_tsv, truthy
from pfam_priority import excluded_reason, interesting_category


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def hotspot_key(value: str) -> tuple[str, int]:
    prefix, _, number = value.rpartition("_")
    return prefix, int(number) if number.isdigit() else 0


def make_pfam_plots(results_dir: Path, annotation_dir: Path, output_prefix: Path, top_n: int = 15) -> None:
    genes = read_tsv(annotation_dir / "hotspot_gene_map.tsv")
    hotspots = read_tsv(results_dir / "hotspots.tsv")
    pfam = parse_pfam(annotation_dir / "pfam.domtblout")
    amr_keys = set(parse_amr(annotation_dir / "amrfinder.tsv"))

    contexts: dict[str, list[dict[str, str]]] = defaultdict(list)
    for gene in genes:
        contexts[gene["gene_key"]].append(gene)

    loci_by_hotspot: dict[str, set[str]] = defaultdict(set)
    for row in hotspots:
        loci_by_hotspot[row["hotspot_id"]].add(row["locus_id"])
    hotspot_ids = sorted(loci_by_hotspot, key=hotspot_key)

    filtered_rows: list[dict[str, object]] = []
    generic_excluded_rows: list[dict[str, object]] = []
    locus_domains: set[tuple[str, str, str]] = set()
    names: dict[str, str] = {}
    descriptions: dict[str, str] = {}
    categories: dict[str, str] = {}
    for gene_key, hits in pfam.items():
        for context in contexts.get(gene_key, []):
            if truthy(context["is_defence"]) or gene_key in amr_keys:
                continue
            for hit in hits:
                accession = hit["accession"].split(".")[0]
                row = {
                    "hotspot_id": context["hotspot_id"],
                    "locus_id": context["locus_id"],
                    "genome": context["genome"],
                    "gene_id": context["gene_id"],
                    "gene_key": gene_key,
                    "pfam_accession": accession,
                    "pfam_name": hit["name"],
                    "pfam_description": hit["description"],
                    "i_evalue": hit["i_evalue"],
                    "score": hit["score"],
                }
                reason = excluded_reason(hit)
                category = interesting_category(hit) if not reason else ""
                if reason or not category:
                    row["exclusion_reason"] = reason or "not_prioritised_accessory_function"
                    generic_excluded_rows.append(row)
                    continue
                row["candidate_category"] = category
                names[accession] = hit["name"]
                descriptions[accession] = hit["description"]
                categories[accession] = category
                locus_domains.add((context["hotspot_id"], context["locus_id"], accession))
                filtered_rows.append(row)
    write_tsv(
        annotation_dir / "pfam_fallback_annotations.tsv",
        filtered_rows,
        [
            "hotspot_id", "locus_id", "genome", "gene_id", "gene_key", "pfam_accession",
            "pfam_name", "pfam_description", "candidate_category", "i_evalue", "score",
        ],
    )
    write_tsv(
        annotation_dir / "pfam_domains_excluded_from_priority.tsv",
        generic_excluded_rows,
        [
            "hotspot_id", "locus_id", "genome", "gene_id", "gene_key", "pfam_accession",
            "pfam_name", "pfam_description", "exclusion_reason", "i_evalue", "score",
        ],
    )
    (annotation_dir / "pfam_generic_domains_excluded.tsv").unlink(missing_ok=True)

    overall = Counter(accession for _hotspot, _locus, accession in locus_domains)
    by_hotspot = Counter((hotspot, accession) for hotspot, _locus, accession in locus_domains)
    frequency_rows: list[dict[str, object]] = []
    for hotspot_id in hotspot_ids:
        total = len(loci_by_hotspot[hotspot_id])
        for accession in sorted(overall):
            count = by_hotspot[(hotspot_id, accession)]
            frequency_rows.append(
                {
                    "hotspot_id": hotspot_id,
                    "pfam_accession": accession,
                    "pfam_name": names[accession],
                    "pfam_description": descriptions[accession],
                    "candidate_category": categories[accession],
                    "positive_loci": count,
                    "total_loci": total,
                    "positive_percent": round(100 * count / total, 4) if total else 0,
                }
            )
    write_tsv(
        annotation_dir / "pfam_by_hotspot.tsv",
        frequency_rows,
        [
            "hotspot_id", "pfam_accession", "pfam_name", "pfam_description", "candidate_category",
            "positive_loci", "total_loci", "positive_percent",
        ],
    )

    ranked = sorted(overall, key=lambda accession: (-overall[accession], names[accession], accession))
    selected = ranked[: max(1, top_n)]
    matrix = np.array(
        [
            [100 * by_hotspot[(hotspot_id, accession)] / len(loci_by_hotspot[hotspot_id]) for accession in selected]
            for hotspot_id in hotspot_ids
        ]
    )
    labels = [f"{names[accession]}\n{accession}" for accession in selected]

    sns.set_theme(style="white", context="talk", font_scale=0.82)
    fig = plt.figure(figsize=(18, 9), constrained_layout=True)
    grid = fig.add_gridspec(1, 2, width_ratios=[0.72, 1.55])

    ax = fig.add_subplot(grid[0, 0])
    order = list(reversed(selected))
    values = [overall[accession] for accession in order]
    y = np.arange(len(order))
    ax.hlines(y, 0, values, color="#26395F", linewidth=2.2)
    ax.scatter(values, y, s=72, color="#26395F", edgecolor="white", linewidth=1)
    ax.set_yticks(y, [f"{names[a]} ({a})" for a in order], fontsize=10)
    ax.set_xlabel("Hotspot loci carrying domain")
    ax.set_title("Most frequent prioritised Pfam candidates", loc="left", weight="bold", pad=14)
    ax.grid(False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    xmax = max(values) if values else 1
    ax.set_xlim(0, xmax * 1.18)
    for yi, value in zip(y, values):
        ax.text(value + xmax * 0.018, yi, str(value), va="center", fontsize=9, color="#343A40")

    ax2 = fig.add_subplot(grid[0, 1])
    annotations = np.full(matrix.shape, "", dtype=object)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i, j] >= 1:
                annotations[i, j] = f"{matrix[i, j]:.0f}"
    sns.heatmap(
        matrix,
        ax=ax2,
        cmap=sns.blend_palette(["#FFFAF5", "#F4C49B", "#D77A55", "#26395F"], as_cmap=True),
        vmin=0,
        vmax=100,
        xticklabels=labels,
        yticklabels=hotspot_ids,
        linewidths=0,
        annot=annotations,
        fmt="",
        annot_kws={"fontsize": 8},
        cbar_kws={"label": "Loci carrying domain (%)", "shrink": 0.66, "pad": 0.02},
    )
    ax2.set_title("Prioritised Pfam candidates by hotspot", loc="left", weight="bold", pad=14)
    ax2.set_xlabel("Pfam domain")
    ax2.set_ylabel("")
    ax2.tick_params(axis="x", rotation=60, length=0, labelsize=9)
    ax2.tick_params(axis="y", rotation=0, length=0)

    fig.suptitle("Interesting Pfam candidates after defence, AMR and mobility filtering", fontsize=20, weight="bold")
    fig.savefig(output_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    findings = [
        f"Prioritised Pfam gene-domain rows: {len(filtered_rows)}",
        f"Unique prioritised Pfam domains: {len(overall)}",
        f"Unique hotspot-locus/domain combinations: {len(locus_domains)}",
        "",
        "Highest-frequency domain in each hotspot:",
    ]
    for hotspot_id in hotspot_ids:
        candidates = [
            (100 * by_hotspot[(hotspot_id, accession)] / len(loci_by_hotspot[hotspot_id]), accession)
            for accession in overall
            if by_hotspot[(hotspot_id, accession)]
        ]
        if candidates:
            percentage, accession = max(candidates)
            findings.append(
                f"{hotspot_id}: {names[accession]} ({accession}), "
                f"{by_hotspot[(hotspot_id, accession)]}/{len(loci_by_hotspot[hotspot_id])} loci "
                f"({percentage:.1f}%)"
            )
        else:
            findings.append(f"{hotspot_id}: no retained Pfam domains")
    findings.extend(["", "Prioritised categories:"])
    for category in sorted(set(categories.values())):
        members = [accession for accession in overall if categories[accession] == category]
        combinations = sum(overall[accession] for accession in members)
        top = sorted(members, key=lambda accession: (-overall[accession], names[accession]))[:5]
        findings.append(
            f"{category}: {combinations} hotspot-locus/domain combinations; top domains "
            + ", ".join(f"{names[a]} ({a}, n={overall[a]})" for a in top)
        )
    findings.extend(["", "Selected biologically informative patterns:"])
    informative = [
        "PF07275",  # ArdA anti-restriction
        "PF12648",  # TcpE
        "PF02604",  # Phd/YeFM antitoxin
        "PF02452",  # PemK toxin
        "PF05016",  # ParE toxin
        "PF15781",  # ParE-like toxin
        "PF02661",  # Fic toxin-associated
        "PF20613",  # HipA toxin
        "PF07282",  # Cas12f1-like TNB
        "PF07751",  # Abi
        "PF02534",  # T4SS DNA transfer
        "PF17802",  # SpaA
        "PF04203",  # Sortase
        "PF16555",  # pilin
    ]
    for accession in informative:
        if accession not in overall:
            continue
        profiles = []
        for hotspot_id in hotspot_ids:
            count = by_hotspot[(hotspot_id, accession)]
            if count:
                total = len(loci_by_hotspot[hotspot_id])
                profiles.append(f"{hotspot_id} {count}/{total} ({100 * count / total:.1f}%)")
        findings.append(f"{names[accession]} ({accession}): " + "; ".join(profiles))
    (annotation_dir / "pfam_findings.txt").write_text("\n".join(findings) + "\n")
    print("\n".join(findings))
    print(output_prefix.with_suffix(".png").resolve())
    print(output_prefix.with_suffix(".pdf").resolve())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--annotation-dir", type=Path)
    parser.add_argument("--output-prefix", type=Path)
    parser.add_argument("--top-domains", type=int, default=15)
    args = parser.parse_args()
    results_dir = args.results_dir.resolve()
    annotation_dir = (args.annotation_dir or results_dir / "hotspot_annotations").resolve()
    output_prefix = (args.output_prefix or annotation_dir / "hotspot_pfam_summary").resolve()
    make_pfam_plots(results_dir, annotation_dir, output_prefix, args.top_domains)


if __name__ == "__main__":
    main()
