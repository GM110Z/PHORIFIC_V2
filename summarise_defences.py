#!/usr/bin/env python3
"""Summarise defence-system abundance and distribution across hotspots."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

from hotspot_pipeline import pdf_escape, write_vector_pdf


EMPTY_LABEL = "no defence detected"


def fail(message: str) -> None:
    raise SystemExit(f"error: {message}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        fail(f"missing required table: {path}")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def split_systems(value: str) -> list[str]:
    """Split multiple systems joined by the hotspot pipeline, but retain caller aliases joined by '/'."""
    return sorted({part.strip() for part in value.split(";")
                   if part.strip() and part.strip().lower() != EMPTY_LABEL})


def pdf_text(commands: list[str], height: float, x: float, y: float, value: str,
             size: int = 10, colour=(0.1, 0.1, 0.1), centred: bool = False) -> None:
    safe = value.encode("latin-1", errors="replace").decode("latin-1")
    if centred:
        x -= len(safe) * size * 0.25
    commands.append(
        f"{colour[0]} {colour[1]} {colour[2]} rg BT /F1 {size} Tf "
        f"{x:.1f} {height-y:.1f} Td ({pdf_escape(safe)}) Tj ET"
    )


def pdf_rect(commands: list[str], height: float, x: float, y: float,
             width: float, rect_height: float, fill, stroke=None) -> None:
    stroke = stroke or fill
    commands.append(
        f"{fill[0]} {fill[1]} {fill[2]} rg {stroke[0]} {stroke[1]} {stroke[2]} RG "
        f"{x:.1f} {height-y-rect_height:.1f} {width:.1f} {rect_height:.1f} re B"
    )


def shorten(value: str, limit: int = 49) -> str:
    return value if len(value) <= limit else value[:limit - 3] + "..."


def draw_type_counts(path: Path, rows: list[dict[str, object]], top_n: int,
                     total_calls: int, total_types: int) -> None:
    shown = rows[:top_n]
    width = 1250
    row_height = 25
    height = max(500, 135 + row_height * len(shown))
    label_x, plot_x, plot_width = 30, 390, 600
    commands = ["0.5 w"]
    pdf_text(commands, height, 30, 35, "Defence systems detected across all genomes", 20)
    pdf_text(commands, height, 30, 60,
             f"{total_calls:,} dereplicated calls; {total_types:,} system labels; showing top {len(shown)}", 11,
             (0.35, 0.35, 0.35))
    maximum = max((int(row["system_calls"]) for row in shown), default=1)
    for index, row in enumerate(shown):
        y = 95 + index * row_height
        count = int(row["system_calls"])
        bar_width = plot_width * count / maximum
        pdf_text(commands, height, label_x, y + 12, shorten(str(row["system_type"])), 9)
        pdf_rect(commands, height, plot_x, y, bar_width, 16, (0.18, 0.48, 0.66))
        pdf_text(commands, height, plot_x + bar_width + 8, y + 12,
                 f"{count:,} calls | {int(row['genomes']):,} genomes", 9)
    pdf_text(commands, height, plot_x, height - 22, "Number of dereplicated system calls", 10,
             (0.35, 0.35, 0.35))
    write_vector_pdf(path, width, height, commands)


def heat_colour(value: int, maximum: int) -> tuple[float, float, float]:
    if value <= 0:
        return (0.95, 0.95, 0.95)
    strength = math.log1p(value) / math.log1p(maximum)
    return (0.94 - 0.72 * strength, 0.96 - 0.48 * strength, 1.0 - 0.25 * strength)


def draw_hotspot_heatmap(path: Path, hotspots: list[str], systems: list[str],
                         counts: Counter[tuple[str, str]], top_n: int) -> None:
    ranked = sorted(systems, key=lambda s: (-sum(counts[(h, s)] for h in hotspots), s))[:top_n]
    cell_width, cell_height = 105, 27
    left, top = 355, 145
    width = max(1000, left + cell_width * len(hotspots) + 45)
    height = max(500, top + cell_height * len(ranked) + 65)
    commands = ["0.35 w"]
    pdf_text(commands, height, 30, 35, "Defence-system distribution across hotspots", 20)
    pdf_text(commands, height, 30, 60,
             "Cells show the number of hotspot loci containing each caller-detected system", 11,
             (0.35, 0.35, 0.35))
    maximum = max((counts[(h, s)] for h in hotspots for s in ranked), default=1)
    for column, hotspot in enumerate(hotspots):
        x = left + column * cell_width
        pdf_text(commands, height, x + cell_width / 2, 112, hotspot, 9, centred=True)
    for row_index, system in enumerate(ranked):
        y = top + row_index * cell_height
        pdf_text(commands, height, 30, y + 18, shorten(system, 52), 9)
        for column, hotspot in enumerate(hotspots):
            x = left + column * cell_width
            value = counts[(hotspot, system)]
            colour = heat_colour(value, maximum)
            pdf_rect(commands, height, x, y, cell_width - 3, cell_height - 3,
                     colour, (0.82, 0.82, 0.82))
            text_colour = (1, 1, 1) if value and math.log1p(value) / math.log1p(maximum) > 0.62 else (0.12, 0.12, 0.12)
            pdf_text(commands, height, x + (cell_width - 3) / 2, y + 17,
                     str(value) if value else "-", 9, text_colour, centred=True)
    write_vector_pdf(path, width, height, commands)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarise defence-system counts and their distribution across inferred hotspots")
    parser.add_argument("--results-dir", type=Path, required=True,
                        help="hotspot pipeline output directory")
    parser.add_argument("--outdir", type=Path, default=None,
                        help="report directory (default: RESULTS_DIR/defence_summary)")
    parser.add_argument("--top-types", type=int, default=40,
                        help="maximum system types shown in each PDF (TSVs always contain all types)")
    args = parser.parse_args()
    if args.top_types < 1:
        fail("--top-types must be at least 1")

    results = args.results_dir.resolve()
    out = (args.outdir or results / "defence_summary").resolve()
    out.mkdir(parents=True, exist_ok=True)
    system_rows = read_tsv(results / "systems.dereplicated.tsv")
    hotspot_rows = read_tsv(results / "hotspots.tsv")

    overall_calls = Counter()
    genomes_by_type: dict[str, set[str]] = defaultdict(set)
    caller_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in system_rows:
        system = row["system"]
        overall_calls[system] += 1
        genomes_by_type[system].add(row["genome"])
        callers = set(filter(None, row.get("callers", "").split(",")))
        if "padloc" in callers:
            caller_counts[system]["padloc"] += 1
        if "defensefinder" in callers:
            caller_counts[system]["defensefinder"] += 1
        if {"padloc", "defensefinder"}.issubset(callers):
            caller_counts[system]["both"] += 1

    overall_summary = [
        {
            "system_type": system,
            "system_calls": count,
            "genomes": len(genomes_by_type[system]),
            "padloc_supported": caller_counts[system]["padloc"],
            "defensefinder_supported": caller_counts[system]["defensefinder"],
            "both_callers": caller_counts[system]["both"],
        }
        for system, count in sorted(overall_calls.items(), key=lambda item: (-item[1], item[0]))
    ]
    write_tsv(out / "defence_counts_by_type.tsv",
              ["system_type", "system_calls", "genomes", "padloc_supported",
               "defensefinder_supported", "both_callers"],
              ([row[key] for key in ("system_type", "system_calls", "genomes", "padloc_supported",
                                     "defensefinder_supported", "both_callers")]
               for row in overall_summary))

    hotspot_total = Counter()
    hotspot_positive = Counter()
    hotspot_empty = Counter()
    hotspot_multi = Counter()
    hotspot_system_counts: Counter[tuple[str, str]] = Counter()
    hotspot_system_types: dict[str, set[str]] = defaultdict(set)
    for row in hotspot_rows:
        hotspot = row["hotspot_id"]
        systems = split_systems(row.get("defence_system", ""))
        hotspot_total[hotspot] += 1
        if systems:
            hotspot_positive[hotspot] += 1
        else:
            hotspot_empty[hotspot] += 1
        if len(systems) > 1:
            hotspot_multi[hotspot] += 1
        for system in systems:
            hotspot_system_counts[(hotspot, system)] += 1
            hotspot_system_types[hotspot].add(system)

    hotspots = sorted(hotspot_total)
    hotspot_summary = []
    for hotspot in hotspots:
        total = hotspot_total[hotspot]
        positive = hotspot_positive[hotspot]
        occurrences = sum(hotspot_system_counts[(hotspot, s)] for s in hotspot_system_types[hotspot])
        hotspot_summary.append((
            hotspot, total, positive, hotspot_empty[hotspot],
            f"{positive / total:.6f}" if total else "0.000000",
            len(hotspot_system_types[hotspot]), occurrences, hotspot_multi[hotspot],
        ))
    write_tsv(out / "hotspot_defence_summary.tsv",
              ["hotspot_id", "total_loci", "defence_positive_loci", "boundary_only_loci",
               "defence_positive_fraction", "distinct_system_types", "system_occurrences",
               "multi_system_loci"], hotspot_summary)

    distribution_rows = []
    all_hotspot_systems = sorted({system for _, system in hotspot_system_counts})
    for hotspot in hotspots:
        for system in sorted(hotspot_system_types[hotspot],
                             key=lambda s: (-hotspot_system_counts[(hotspot, s)], s)):
            count = hotspot_system_counts[(hotspot, system)]
            total = hotspot_total[hotspot]
            positive = hotspot_positive[hotspot]
            distribution_rows.append((
                hotspot, system, count, total, positive,
                f"{count / total:.6f}" if total else "0.000000",
                f"{count / positive:.6f}" if positive else "0.000000",
            ))
    write_tsv(out / "defence_distribution_by_hotspot.tsv",
              ["hotspot_id", "system_type", "loci_with_system", "total_hotspot_loci",
               "defence_positive_loci", "fraction_of_all_loci", "fraction_of_defence_loci"],
              distribution_rows)

    draw_type_counts(out / "defence_counts_by_type.pdf", overall_summary, args.top_types,
                     sum(overall_calls.values()), len(overall_calls))
    draw_hotspot_heatmap(out / "defence_distribution_by_hotspot.pdf", hotspots,
                         all_hotspot_systems, hotspot_system_counts, args.top_types)
    print(f"Wrote defence summaries and plots to {out}")


if __name__ == "__main__":
    main()
