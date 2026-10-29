#!/usr/bin/env python3
"""Downstream plasmid/prophage analysis for PhorificV2.

The program is deliberately run *after* the main PhorificV2 analysis.  It
filters geNomad calls, classifies existing hotspot instances by MGE context,
extracts high-confidence plasmids and prophages, translates the existing
PADLOC/DefenseFinder calls into element-local coordinates, and runs independent
MGE-specific hotspot analyses without rerunning either defence caller.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path

from hotspot_pipeline import Gene, pdf_escape, read_fasta, read_gff, write_fasta, write_vector_pdf


NA = {"", "na", "nan", "none", "."}


def fail(message: str) -> None:
    raise SystemExit(f"error: {message}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        fail(f"missing required table: {path}")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def number(value: str, default: float | None = None) -> float | None:
    if value.strip().lower() in NA:
        return default
    try:
        return float(value)
    except ValueError:
        return default


def integer(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def clean_contig(value: str) -> str:
    """Return the input-contig ID embedded in a geNomad sequence name."""
    value = value.split("|provirus_", 1)[0]
    value = re.sub(r"_\[gcode=\d+\]_\[topology=[^]]+\]$", "", value)
    accession = re.match(r"([A-Z]{2}_[A-Z0-9]+\.\d+)", value)
    if accession:
        return accession.group(1)
    return value.split()[0]


def clean_genome(value: str) -> str:
    """Match geNomad prefixes such as GCF_xxx_genomic to Phorific genome IDs."""
    value = value.strip()
    accession = re.match(r"(GC[AF]_\d+\.\d+)", value)
    return accession.group(1) if accession else re.sub(r"_genomic$", "", value)


def safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


@dataclass(frozen=True)
class Region:
    kind: str
    genome: str
    contig: str
    start: int
    end: int
    seq_name: str
    score: float
    topology: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def filter_predictions(rows: list[dict[str, str]], kind: str, args) -> tuple[list[dict[str, str]], list[list[object]]]:
    score_column = "plasmid_score" if kind == "plasmid" else "virus_score"
    kept, audit = [], []
    for row in rows:
        reasons = []
        score = number(row.get(score_column, ""), -1.0)
        fdr = number(row.get("fdr", ""), None)
        enrichment = number(row.get("marker_enrichment", ""), float("-inf"))
        if score is None or score < args.min_score:
            reasons.append(f"{score_column}<{args.min_score:g}")
        if fdr is not None and fdr > args.max_fdr:
            reasons.append(f"fdr>{args.max_fdr:g}")
        if integer(row.get("n_genes", "0")) < args.min_genes:
            reasons.append(f"n_genes<{args.min_genes}")
        if enrichment is None or enrichment < args.min_marker_enrichment:
            reasons.append(f"marker_enrichment<{args.min_marker_enrichment:g}")
        if integer(row.get("n_hallmarks", "0")) < args.min_hallmarks:
            reasons.append(f"n_hallmarks<{args.min_hallmarks}")
        retained = not reasons
        if retained:
            kept.append(row)
        audit.append([kind, row.get("genome", ""), row.get("seq_name", ""), int(retained), ";".join(reasons)])
    return kept, audit


def regions_from_rows(rows: list[dict[str, str]], kind: str, integrated_only: bool = False) -> list[Region]:
    regions = []
    score_column = "plasmid_score" if kind == "plasmid" else "virus_score"
    for row in rows:
        topology = row.get("topology", "")
        if integrated_only and topology.lower() != "provirus":
            continue
        contig = clean_contig(row["seq_name"])
        coordinates = row.get("coordinates", "")
        match = re.fullmatch(r"(\d+)-(\d+)", coordinates) if coordinates.lower() not in NA else None
        if match:
            start, end = map(int, match.groups())
        else:
            start, end = 1, integer(row.get("length", "0"))
        if start < 1 or end < start:
            continue
        regions.append(Region(kind, clean_genome(row["genome"]), contig, start, end, row["seq_name"],
                              float(row[score_column]), topology))
    return regions


def interval_overlap(a1: int, a2: int, b1: int, b2: int) -> int:
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def index_regions(regions: list[Region]) -> dict[tuple[str, str], list[Region]]:
    index: dict[tuple[str, str], list[Region]] = defaultdict(list)
    for region in regions:
        index[(clean_genome(region.genome), clean_contig(region.contig))].append(region)
    return index


def region_matches(regions: dict[tuple[str, str], list[Region]], genome: str, contig: str, start: int, end: int,
                   minimum_fraction: float) -> list[tuple[Region, float]]:
    result = []
    contig = clean_contig(contig)
    length = max(1, end - start + 1)
    for region in regions.get((clean_genome(genome), contig), []):
        fraction = interval_overlap(start, end, region.start, region.end) / length
        if fraction >= minimum_fraction:
            result.append((region, fraction))
    return sorted(result, key=lambda item: (-item[1], -item[0].score, item[0].start))


def classify_hotspots(hotspots: list[dict[str, str]], plasmids: list[Region], phages: list[Region],
                      minimum_fraction: float) -> tuple[list[list[object]], Counter[tuple[str, str]]]:
    output, counts = [], Counter()
    plasmid_index, phage_index = index_regions(plasmids), index_regions(phages)
    for row in hotspots:
        start, end = int(row["start"]), int(row["end"])
        pm = region_matches(plasmid_index, row["genome"], row["contig"], start, end, minimum_fraction)
        vm = region_matches(phage_index, row["genome"], row["contig"], start, end, minimum_fraction)
        # A provirus can be nested in a plasmid.  The plot needs exclusive
        # categories, so the more specific prophage interval takes precedence.
        context = "prophage" if vm else "plasmid" if pm else "chromosome"
        counts[(row["hotspot_id"], context)] += 1
        best = vm[0] if vm else pm[0] if pm else None
        output.append([
            row["hotspot_id"], row["locus_id"], row["genome"], row["contig"], start, end,
            row.get("defence_system", ""), row.get("defence_present", ""), context,
            int(bool(pm)), int(bool(vm)), int(bool(pm and vm)),
            f"{best[1]:.6f}" if best else "0.000000", best[0].seq_name if best else "",
        ])
    return output, counts


def pdf_text(commands, height, x, y, value, size=10, colour=(0.1, 0.1, 0.1), centred=False):
    value = str(value).encode("latin-1", errors="replace").decode("latin-1")
    if centred:
        x -= len(value) * size * 0.25
    commands.append(f"{colour[0]} {colour[1]} {colour[2]} rg BT /F1 {size} Tf "
                    f"{x:.1f} {height-y:.1f} Td ({pdf_escape(value)}) Tj ET")


def pdf_rect(commands, height, x, y, width, rectangle_height, colour):
    commands.append(f"{colour[0]} {colour[1]} {colour[2]} rg "
                    f"{x:.1f} {height-y-rectangle_height:.1f} {width:.1f} {rectangle_height:.1f} re f")


def draw_context_plot(path: Path, counts: Counter[tuple[str, str]]) -> None:
    hotspots = sorted({hotspot for hotspot, _ in counts})
    width, row_height = 1050, 42
    height = max(470, 205 + row_height * len(hotspots))
    left, plot_width = 185, 760
    # Publication palette without a red/green pairing: muted indigo,
    # peach-orange, and blue-grey, separated in hue and luminance.
    colours = {"prophage": (0.31, 0.39, 0.57), "plasmid": (0.95, 0.63, 0.40),
               "chromosome": (0.58, 0.65, 0.73)}
    commands = []
    pdf_text(commands, height, 32, 38, "Genomic context of Phorific hotspot instances", 20)
    pdf_text(commands, height, 32, 64,
             "High-confidence geNomad calls; prophage takes precedence over plasmid when predictions overlap",
             10, (0.35, 0.35, 0.35))
    legend_x = 510
    for index, category in enumerate(("prophage", "plasmid", "chromosome")):
        x = legend_x + index * 145
        pdf_rect(commands, height, x, 82, 18, 12, colours[category])
        pdf_text(commands, height, x + 25, 93, category, 10)
    for index, hotspot in enumerate(hotspots):
        y = 125 + index * row_height
        total = sum(counts[(hotspot, category)] for category in colours)
        pdf_text(commands, height, 32, y + 18, hotspot, 10)
        x = left
        for category in ("prophage", "plasmid", "chromosome"):
            value = counts[(hotspot, category)]
            bar_width = plot_width * value / total if total else 0
            if bar_width:
                pdf_rect(commands, height, x, y, bar_width, 25, colours[category])
                if bar_width >= 25:
                    label_colour = (0.18, 0.15, 0.13) if category == "plasmid" else (1, 1, 1)
                    pdf_text(commands, height, x + bar_width / 2, y + 17, value, 9, label_colour, centred=True)
            x += bar_width
        pdf_text(commands, height, left + plot_width + 12, y + 18, f"n={total}", 9, (0.35, 0.35, 0.35))
    axis_y = 133 + len(hotspots) * row_height
    top_axis_y = 125
    # Explicit categorical y-axis with one tick for every hotspot.
    commands.append(f"0.32 0.32 0.32 RG 0.7 w {left:.1f} {height-top_axis_y:.1f} m "
                    f"{left:.1f} {height-axis_y:.1f} l S")
    for index in range(len(hotspots)):
        tick_y = 125 + index * row_height + 12.5
        commands.append(f"0.32 0.32 0.32 RG 0.7 w {left-5:.1f} {height-tick_y:.1f} m "
                        f"{left:.1f} {height-tick_y:.1f} l S")
    y_title = "Hotspot architecture"
    title_midpoint = (top_axis_y + axis_y) / 2
    title_baseline = height - title_midpoint - len(y_title) * 10 * 0.25
    commands.append(f"0.20 0.20 0.20 rg BT /F1 10 Tf 0 1 -1 0 18 {title_baseline:.1f} Tm "
                    f"({pdf_escape(y_title)}) Tj ET")
    # Bottom percentage axis only: publication-ready and intentionally free of
    # the distracting full-height gridlines used by spreadsheet defaults.
    pdf_y = height - axis_y
    commands.append(f"0.32 0.32 0.32 RG 0.7 w {left:.1f} {pdf_y:.1f} m "
                    f"{left+plot_width:.1f} {pdf_y:.1f} l S")
    for value in (0, 25, 50, 75, 100):
        x = left + plot_width * value / 100
        commands.append(f"0.32 0.32 0.32 RG 0.7 w {x:.1f} {pdf_y:.1f} m {x:.1f} {pdf_y-5:.1f} l S")
        pdf_text(commands, height, x, axis_y + 18, str(value), 9, (0.28, 0.28, 0.28), centred=True)
    pdf_text(commands, height, left + plot_width/2, axis_y + 42,
             "Percentage of hotspot instances (%)", 10, (0.22, 0.22, 0.22), centred=True)
    write_vector_pdf(path, width, height, commands)


def accession_key(path: Path) -> str:
    match = re.search(r"(GC[AF]_\d+\.\d+)", str(path))
    if match:
        return match.group(1)
    stem = path.name
    for suffix in (".gff3", ".gff", ".faa", ".fna", ".fa", ".fasta"):
        if stem.endswith(suffix):
            stem = stem[:-len(suffix)]
    return re.sub(r"_(genomic|protein)$", "", stem)


def discover_sources(root: Path, genomes: set[str]) -> dict[str, dict[str, Path]]:
    result: dict[str, dict[str, Path]] = defaultdict(dict)
    suffix_kind = {".gff3": "gff", ".gff": "gff", ".faa": "faa",
                   ".fna": "fna", ".fa": "fna", ".fasta": "fna"}
    candidates: dict[tuple[str, str], list[Path]] = defaultdict(list)
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in suffix_kind:
            continue
        genome = accession_key(path)
        if genome in genomes:
            candidates[(genome, suffix_kind[path.suffix.lower()])].append(path)
    for (genome, kind), paths in candidates.items():
        def rank(path: Path):
            exact = path.stem in {genome, f"{genome}_genomic", f"{genome}_protein"}
            return (not exact, len(str(path)), str(path))
        result[genome][kind] = sorted(paths, key=rank)[0]
    return result


def existing_resolved_sources(path: Path, genomes: set[str]) -> dict[str, dict[str, Path]]:
    result: dict[str, dict[str, Path]] = defaultdict(dict)
    if not path.is_file():
        return result
    for row in read_tsv(path):
        genome = row.get("genome_id", "")
        if genome not in genomes:
            continue
        for kind in ("gff", "faa", "fna"):
            candidate = Path(row.get(kind, "")) if row.get(kind) else None
            if candidate and candidate.is_file():
                result[genome][kind] = candidate
    return result


def choose_contig(records: dict[str, str], wanted: str) -> str | None:
    wanted = clean_contig(wanted)
    if wanted in records:
        return wanted
    matches = [name for name in records if clean_contig(name) == wanted]
    return matches[0] if len(matches) == 1 else None


def write_extracted_gff(path: Path, mge_id: str, genes: list[Gene], length: int) -> None:
    with path.open("w") as handle:
        handle.write("##gff-version 3\n")
        handle.write(f"##sequence-region {mge_id} 1 {length}\n")
        for gene in genes:
            product = gene.product.replace(";", ",").replace("\t", " ")
            handle.write(f"{mge_id}\tPhorificMGE\tCDS\t{gene.start}\t{gene.end}\t.\t{gene.strand}\t.\t"
                         f"ID={gene.gene_id};protein_id={gene.gene_id};product={product}\n")


def materialize_parent_clusters(path: Path, element_roots: list[Path]) -> None:
    """Project the large parent MMseqs table onto only extracted MGE genes.

    Holding all 7+ million parent members in a Python dictionary is wasteful.
    The extraction stage writes a compact original/local UID request table;
    this function scans the parent table once and writes the two small MGE
    cluster tables directly.
    """
    if not path.is_file():
        fail(f"parent MMseqs cluster table not found: {path}")
    requests: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for index, root in enumerate(element_roots):
        request_path = root / "cluster_requests.tsv"
        for row in read_tsv(request_path):
            requests[row["original_uid"]].append((index, row["local_uid"]))
    handles = [(root / "precomputed_clusters.tsv").open("w") for root in element_roots]
    found = Counter()
    try:
        with path.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2 or fields[1] not in requests:
                    continue
                representative, original_uid = fields[0], fields[1]
                for index, local_uid in requests[original_uid]:
                    handles[index].write(f"{representative}\t{local_uid}\n")
                    found[index] += 1
    finally:
        for handle in handles:
            handle.close()
    for index, root in enumerate(element_roots):
        (root / "cluster_projection_summary.json").write_text(json.dumps({
            "requested_members": sum(1 for values in requests.values()
                                     for destination, _ in values if destination == index),
            "matched_members": found[index],
        }, indent=2) + "\n")


def extract_elements(kind: str, regions: list[Region], sources: dict[str, dict[str, Path]],
                     systems: list[dict[str, str]], root: Path):
    element_root = root / ("plasmids" if kind == "plasmid" else "prophages")
    element_root.mkdir(parents=True, exist_ok=True)
    manifest_rows, system_rows, cluster_requests, metadata_rows, skipped = [], [], [], [], []
    system_index: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in systems:
        system_index[(row["genome"], clean_contig(row["contig"]))].append(row)
    by_genome: dict[str, list[Region]] = defaultdict(list)
    for region in regions:
        by_genome[region.genome].append(region)
    counter = 0
    for genome, genome_regions in sorted(by_genome.items()):
        files = sources.get(genome, {})
        missing = [kind_ for kind_ in ("gff", "faa") if kind_ not in files]
        if missing:
            skipped.extend([r.seq_name, genome, ",".join(missing), "missing source file"] for r in genome_regions)
            continue
        proteins = read_fasta(files["faa"])
        genes, _ = read_gff(files["gff"], genome, proteins)
        sequences = read_fasta(files["fna"]) if "fna" in files else {}
        genes_by_contig: dict[str, list[Gene]] = defaultdict(list)
        for gene in genes:
            genes_by_contig[clean_contig(gene.contig)].append(gene)
        for region in sorted(genome_regions, key=lambda r: (r.contig, r.start, r.end)):
            annotation_contig = clean_contig(region.contig)
            if annotation_contig not in genes_by_contig:
                skipped.append([region.seq_name, genome, "gff contig", f"contig not found: {region.contig}"])
                continue
            source_contig = choose_contig(sequences, region.contig) if sequences else None
            sequence = sequences[source_contig] if source_contig else ""
            end = min(region.end, len(sequence)) if sequence else region.end
            selected_original = [g for g in genes_by_contig[annotation_contig]
                                 if g.start >= region.start and g.end <= end]
            selected_original = [g for g in selected_original if g.sequence]
            if not selected_original:
                skipped.append([region.seq_name, genome, "annotation", "no fully contained proteins"])
                continue
            counter += 1
            mge_id = f"{kind}_{counter:06d}__{safe_id(genome)}__{safe_id(region.contig)}"
            target = element_root / mge_id
            target.mkdir(parents=True, exist_ok=True)
            local_genes = [replace(g, genome=mge_id, contig=mge_id,
                                   start=g.start-region.start+1, end=g.end-region.start+1)
                           for g in selected_original]
            fna, faa, gff = target / f"{mge_id}.fna", target / f"{mge_id}.faa", target / f"{mge_id}.gff3"
            fna_value = ""
            if sequence:
                write_fasta([(mge_id, sequence[region.start-1:end])], fna)
                fna_value = str(fna.resolve())
            write_fasta(((g.gene_id, g.sequence) for g in local_genes), faa)
            write_extracted_gff(gff, mge_id, local_genes, end-region.start+1)
            manifest_rows.append([mge_id, fna_value, str(gff.resolve()), str(faa.resolve())])
            metadata_rows.append([mge_id, kind, genome, region.contig, region.start, end,
                                  region.seq_name, region.topology, f"{region.score:.6f}"])
            for original, local in zip(selected_original, local_genes):
                cluster_requests.append([original.uid, local.uid])
            for system in system_index[(genome, clean_contig(region.contig))]:
                s1, s2 = int(system["start"]), int(system["end"])
                overlap = interval_overlap(region.start, end, s1, s2)
                if not overlap or overlap / max(1, s2-s1+1) < 0.5:
                    continue
                local_start = max(s1, region.start) - region.start + 1
                local_end = min(s2, end) - region.start + 1
                system_rows.append([system.get("locus_id", ""), mge_id, mge_id, local_start, local_end,
                                    system["system"], system.get("callers", ""), system.get("source_ids", "")])
    write_tsv(element_root / "manifest.tsv", ["genome_id", "fna", "gff", "faa"], manifest_rows)
    write_tsv(element_root / "precomputed_systems.tsv",
              ["locus_id", "genome", "contig", "start", "end", "system", "callers", "source_ids"], system_rows)
    write_tsv(element_root / "cluster_requests.tsv", ["original_uid", "local_uid"], cluster_requests)
    write_tsv(element_root / "elements.tsv",
              ["mge_id", "mge_type", "source_genome", "source_contig", "source_start", "source_end",
               "genomad_seq_name", "topology", "score"], metadata_rows)
    write_tsv(element_root / "skipped_elements.tsv", ["seq_name", "genome", "problem", "detail"], skipped)
    return element_root, len(manifest_rows), len(system_rows)


def write_old_boundary_pairs(hotspots: list[dict[str, str]], path: Path) -> None:
    pairs = sorted({(row["left_boundary_cluster"], row["right_boundary_cluster"]) for row in hotspots})
    write_tsv(path, ["left_boundary_cluster", "right_boundary_cluster"], pairs)


def run_mge_hotspots(kind: str, element_root: Path, output: Path, old_pairs: Path, args) -> None:
    systems = read_tsv(element_root / "precomputed_systems.tsv")
    if not systems:
        print(f"No reused defence calls occur in retained {kind} elements; skipping hotspot inference.")
        return
    element_count = len(read_tsv(element_root / "manifest.tsv"))
    min_hotspot_members = (args.min_hotspot_members if args.min_hotspot_members is not None
                           else max(1, math.ceil(element_count * args.min_hotspot_prevalence)))
    print(f"Requiring each {kind} hotspot in at least {min_hotspot_members}/{element_count} "
          f"elements ({100 * min_hotspot_members / element_count:.2f}%).")
    command = [
        sys.executable, str(Path(__file__).with_name("phorificV2.py")),
        "--manifest", str(element_root / "manifest.tsv"),
        "--outdir", str(output),
        "--precomputed-systems", str(element_root / "precomputed_systems.tsv"),
        "--precomputed-clusters", str(element_root / "precomputed_clusters.tsv"),
        "--exclude-boundary-pairs", str(old_pairs),
        "--hotspot-prefix", f"{kind}_hotspot",
        "--plot-palette", args.plot_palette,
        "--window-kb", str(args.window_kb),
        "--threads", str(args.threads),
        "--min-boundary-frequency", str(args.min_boundary_frequency),
        "--min-comparable-loci", str(args.min_comparable_loci),
        "--min-hotspot-members", str(min_hotspot_members),
        "--boundary-scan-span-factor", str(args.boundary_scan_span_factor),
        "--representation-only",
    ]
    print("Running:", " ".join(command))
    subprocess.run(command, check=True)
    if args.include_edge_hotspots:
        edge_command = [
            sys.executable, str(Path(__file__).with_name("plot_edge_associated_architecture.py")),
            "--element-root", str(element_root),
            "--outdir", str(output / "loci"),
            "--hotspot-prefix", f"{kind}_edge_hotspot",
            "--window-kb", str(args.window_kb),
            "--min-boundary-frequency", str(args.min_boundary_frequency),
            "--min-comparable-loci", str(args.min_comparable_loci),
            "--min-hotspot-prevalence", str(args.min_hotspot_prevalence),
            "--min-hotspot-members", str(min_hotspot_members),
            "--boundary-scan-span-factor", str(args.boundary_scan_span_factor),
            "--plot-palette", args.plot_palette,
        ]
        print("Running experimental edge hotspots:", " ".join(edge_command))
        subprocess.run(edge_command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="PhorificMGE: geNomad-aware downstream analysis for PhorificV2")
    parser.add_argument("--phorific-results", type=Path, required=True,
                        help="completed main PhorificV2 output directory")
    parser.add_argument("--plasmids-tsv", type=Path, required=True,
                        help="combined geNomad plasmid summary")
    parser.add_argument("--phages-tsv", type=Path, required=True,
                        help="combined geNomad virus summary")
    parser.add_argument("--genome-files-dir", type=Path,
                        help="recursive source of FNA, GFF/GFF3 and FAA files (required unless --analysis-only)")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--genome", action="append", default=[],
                        help="optional genome ID to retain; repeatable (primarily for testing/subsets)")
    parser.add_argument("--analysis-only", action="store_true",
                        help="filter calls, classify existing hotspots and plot; do not extract or rerun Phorific")
    parser.add_argument("--prepare-only", action="store_true",
                        help="also extract MGEs, but do not run MGE hotspot inference")
    # geNomad conservative preset values available in the master summary tables.
    parser.add_argument("--min-score", type=float, default=0.80)
    parser.add_argument("--max-fdr", type=float, default=0.05)
    parser.add_argument("--min-genes", type=int, default=1)
    parser.add_argument("--min-marker-enrichment", type=float, default=1.50)
    parser.add_argument("--min-hallmarks", type=int, default=1)
    parser.add_argument("--min-hotspot-overlap", type=float, default=0.50,
                        help="minimum fraction of a hotspot interval inside an MGE prediction")
    parser.add_argument("--phage-extraction", choices=("provirus", "all"), default="provirus",
                        help="extract integrated proviruses only (default) or all viral predictions")
    parser.add_argument("--window-kb", type=float, default=20.0)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-boundary-frequency", type=float, default=0.5)
    parser.add_argument("--min-comparable-loci", type=int, default=3)
    parser.add_argument("--min-hotspot-prevalence", type=float, default=0.05,
                        help="minimum fraction of extracted elements carrying an architecture (default: 0.05)")
    parser.add_argument("--min-hotspot-members", type=int,
                        help="optional fixed member cutoff; overrides --min-hotspot-prevalence")
    parser.add_argument("--boundary-scan-span-factor", type=float, default=5.0)
    parser.add_argument("--include-edge-hotspots", action="store_true",
                        help="also run experimental one-edge hotspots (off by default)")
    parser.add_argument("--plot-palette", choices=("classic", "peach"), default="classic",
                        help="locus colours; classic matches original Phorific (blue boundaries, red defence)")
    args = parser.parse_args()
    if not (0 <= args.min_score <= 1 and 0 <= args.max_fdr <= 1
            and 0 < args.min_hotspot_overlap <= 1 and 0 < args.min_hotspot_prevalence <= 1
            and (args.min_hotspot_members is None or args.min_hotspot_members >= 1)):
        fail("filter thresholds are out of range")
    results, out = args.phorific_results.resolve(), args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    hotspots = read_tsv(results / "hotspots.tsv")
    systems = read_tsv(results / "systems.dereplicated.tsv")
    selected_genomes = {clean_genome(value) for value in args.genome}
    raw_plasmids, raw_phages = read_tsv(args.plasmids_tsv), read_tsv(args.phages_tsv)
    if selected_genomes:
        raw_plasmids = [row for row in raw_plasmids if clean_genome(row["genome"]) in selected_genomes]
        raw_phages = [row for row in raw_phages if clean_genome(row["genome"]) in selected_genomes]
        hotspots = [row for row in hotspots if clean_genome(row["genome"]) in selected_genomes]
        systems = [row for row in systems if clean_genome(row["genome"]) in selected_genomes]
    plasmid_rows, plasmid_audit = filter_predictions(raw_plasmids, "plasmid", args)
    phage_rows, phage_audit = filter_predictions(raw_phages, "phage", args)
    if not plasmid_rows or not phage_rows:
        fail("conservative filtering removed every prediction from one or both MGE classes")
    write_tsv(out / "filtered_plasmids.tsv", list(plasmid_rows[0]),
              ([row[key] for key in plasmid_rows[0]] for row in plasmid_rows))
    write_tsv(out / "filtered_phages.tsv", list(phage_rows[0]),
              ([row[key] for key in phage_rows[0]] for row in phage_rows))
    write_tsv(out / "filter_audit.tsv", ["mge_type", "genome", "seq_name", "retained", "rejection_reasons"],
              plasmid_audit + phage_audit)
    plasmid_regions = regions_from_rows(plasmid_rows, "plasmid")
    phage_regions = regions_from_rows(phage_rows, "phage")
    context_rows, context_counts = classify_hotspots(
        hotspots, plasmid_regions, phage_regions, args.min_hotspot_overlap)
    write_tsv(out / "hotspot_mge_context.tsv",
              ["hotspot_id", "locus_id", "genome", "contig", "start", "end", "defence_system",
               "defence_present", "mge_context", "plasmid_match", "phage_match", "overlap_conflict",
               "best_overlap_fraction", "matched_genomad_sequence"], context_rows)
    summary_rows = []
    for hotspot in sorted({row["hotspot_id"] for row in hotspots}):
        values = [context_counts[(hotspot, category)] for category in ("prophage", "plasmid", "chromosome")]
        summary_rows.append([hotspot, *values, sum(values)])
    write_tsv(out / "hotspot_mge_context_summary.tsv",
              ["hotspot_id", "prophage", "plasmid", "chromosome", "total"], summary_rows)
    draw_context_plot(out / "hotspot_mge_context_stacked.pdf", context_counts)
    write_tsv(out / "filter_summary.tsv", ["mge_type", "input_predictions", "retained_predictions", "removed_predictions"], [
        ["plasmid", len(plasmid_audit), len(plasmid_rows), len(plasmid_audit)-len(plasmid_rows)],
        ["phage", len(phage_audit), len(phage_rows), len(phage_audit)-len(phage_rows)],
    ])
    summary = {
        "plasmids_input": len(plasmid_audit), "plasmids_retained": len(plasmid_rows),
        "phages_input": len(phage_audit), "phages_retained": len(phage_rows),
        "hotspot_instances": len(hotspots),
        "hotspot_contexts": {category: sum(value for (hotspot, key), value in context_counts.items() if key == category)
                             for category in ("prophage", "plasmid", "chromosome")},
        "filter_note": "Matches geNomad conservative thresholds available in summary TSVs; max-USCG cannot be reapplied without gene tables.",
    }
    if args.analysis_only:
        (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary, indent=2))
        return
    if not args.genome_files_dir or not args.genome_files_dir.is_dir():
        fail("--genome-files-dir is required for MGE extraction")
    extraction_phages = regions_from_rows(phage_rows, "phage", integrated_only=args.phage_extraction == "provirus")
    genomes = {r.genome for r in plasmid_regions + extraction_phages}
    sources = existing_resolved_sources(results / "resolved_inputs.tsv", genomes)
    discovered = discover_sources(args.genome_files_dir.resolve(), genomes)
    for genome, files in discovered.items():
        sources[genome].update({kind: path for kind, path in files.items() if kind not in sources[genome]})
    extraction_root = out / "extracted"
    plasmid_root, plasmids_extracted, plasmid_systems = extract_elements(
        "plasmid", plasmid_regions, sources, systems, extraction_root)
    phage_root, phages_extracted, phage_systems = extract_elements(
        "prophage", extraction_phages, sources, systems, extraction_root)
    materialize_parent_clusters(results / "mmseqs" / "all_proteins_cluster.tsv", [plasmid_root, phage_root])
    old_pairs = out / "existing_chromosomal_boundary_pairs.tsv"
    write_old_boundary_pairs(hotspots, old_pairs)
    summary.update({"plasmids_extracted": plasmids_extracted, "prophages_extracted": phages_extracted,
                    "reused_plasmid_defence_calls": plasmid_systems,
                    "reused_prophage_defence_calls": phage_systems})
    if not args.prepare_only:
        run_mge_hotspots("plasmid", plasmid_root, out / "plasmid_hotspots", old_pairs, args)
        run_mge_hotspots("prophage", phage_root, out / "prophage_hotspots", old_pairs, args)
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
