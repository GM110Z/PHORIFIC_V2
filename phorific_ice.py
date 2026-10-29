#!/usr/bin/env python3
"""Analyse ICEfinder-opt calls downstream of a completed PhorificV2 run.

The workflow counts ICE/IME/rICE calls, maps ICEfinder's ordinal contig names
back to the source GFF contigs, tests existing Phorific hotspots for overlap,
extracts single-contig elements, projects existing defence calls and MMseqs
clusters, and discovers ICE-specific hotspot architectures without rerunning
PADLOC, DefenseFinder, or MMseqs2.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path

from hotspot_pipeline import Gene, read_fasta, read_gff, write_fasta
from phorific_mge import (
    Region, accession_key, clean_contig, clean_genome, discover_sources,
    existing_resolved_sources, index_regions, interval_overlap,
    materialize_parent_clusters, read_tsv, region_matches, run_mge_hotspots,
    safe_id, write_extracted_gff, write_old_boundary_pairs, write_tsv,
)


TYPE_RE = re.compile(r"_(rICE|AICE|ICE|IME)_\d+$", re.IGNORECASE)
SIMPLE_LOCATION_RE = re.compile(r"^([^:;]+):\(?(\d+)\.\.(\d+)\)?$")
RANGED_PART_RE = re.compile(r"^([^:;]+):\(?(-?\d+)\.\.(-?\d+)\)?$")
ORDINAL_CONTIG_RE = re.compile(r"^contig_(\d+)$", re.IGNORECASE)


def fail(message: str) -> None:
    raise SystemExit(f"error: {message}")


def element_type(mge_id: str) -> str:
    match = TYPE_RE.search(mge_id)
    if not match:
        return "Other"
    value = match.group(1)
    return "rICE" if value.lower() == "rice" else value.upper()


@dataclass(frozen=True)
class Prediction:
    genome: str
    mge_id: str
    mge_type: str
    location: str
    segments: tuple[Region, ...]


def gff_contigs(path: Path) -> tuple[list[str], dict[str, int]]:
    ordered: list[str] = []
    lengths: dict[str, int] = {}
    seen: set[str] = set()
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("##sequence-region"):
                fields = line.rstrip("\n").split()
                if len(fields) >= 4:
                    contig = clean_contig(fields[1])
                    if contig not in seen:
                        ordered.append(contig)
                        seen.add(contig)
                    try:
                        lengths[contig] = max(lengths.get(contig, 0), int(fields[3]))
                    except ValueError:
                        pass
            elif not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 5:
                    contig = clean_contig(fields[0])
                    if contig not in seen:
                        ordered.append(contig)
                        seen.add(contig)
                    try:
                        lengths[contig] = max(lengths.get(contig, 0), int(fields[4]))
                    except ValueError:
                        pass
    return ordered, lengths


def resolve_contig(token: str, ordered: list[str], lengths: dict[str, int]) -> str | None:
    token = clean_contig(token.strip())
    # ICEfinder-cleaned FASTA headers for the original strains retain metadata,
    # e.g. contig_1_length=2814993_circular=true, while their GFF uses contig_1.
    token = re.sub(r"_length=\d+(?:_circular=(?:true|false))?$", "", token,
                   flags=re.IGNORECASE)
    ordinal = ORDINAL_CONTIG_RE.fullmatch(token)
    if ordinal:
        index = int(ordinal.group(1)) - 1
        return ordered[index] if 0 <= index < len(ordered) else None
    if token in lengths:
        return token
    matches = [name for name in ordered if clean_contig(name) == token]
    return matches[0] if len(matches) == 1 else None


def parse_segments(genome: str, mge_id: str, mge_type: str, location: str,
                   ordered: list[str], lengths: dict[str, int]) -> tuple[list[Region], str]:
    parts = [part.strip() for part in location.split(";") if part.strip()]
    regions: list[Region] = []
    for part_index, part in enumerate(parts, 1):
        match = RANGED_PART_RE.fullmatch(part)
        if match:
            raw_contig, start_text, end_text = match.groups()
            start, end = int(start_text), int(end_text)
        else:
            raw_contig = part
            start = 1
            end = 0
        contig = resolve_contig(raw_contig, ordered, lengths)
        if contig is None:
            return [], f"unresolved contig: {raw_contig}"
        if end == 0:
            end = lengths.get(contig, 0)
        # A small number of fragmented rICE calls extend beyond an assembly
        # edge. Clamp those fragments to the actual contig rather than losing
        # the entire biological prediction.
        start = max(1, start)
        end = min(end, lengths.get(contig, end))
        if end < start:
            return [], f"invalid coordinates: {part}"
        regions.append(Region(
            kind=mge_type, genome=genome, contig=contig, start=start, end=end,
            seq_name=f"{mge_id}:part{part_index}", score=1.0, topology=mge_id,
        ))
    return regions, "" if regions else "no location segments"


def parse_predictions(rows: list[dict[str, str]], sources: dict[str, dict[str, Path]]) \
        -> tuple[list[Prediction], list[list[object]]]:
    cache: dict[str, tuple[list[str], dict[str, int]]] = {}
    predictions: list[Prediction] = []
    audit: list[list[object]] = []
    for row in rows:
        raw_genome = row.get("Isolate", "").strip()
        genome = clean_genome(raw_genome)
        mge_id = row.get("MGE", "").strip()
        location = row.get("Location", "").strip()
        mge_type = element_type(mge_id)
        gff = sources.get(genome, {}).get("gff")
        if not gff:
            audit.append([raw_genome, genome, mge_id, mge_type, location, 0, 0, "missing GFF"])
            continue
        if genome not in cache:
            cache[genome] = gff_contigs(gff)
        segments, problem = parse_segments(genome, mge_id, mge_type, location, *cache[genome])
        audit.append([raw_genome, genome, mge_id, mge_type, location, int(bool(segments)),
                      len(segments), problem])
        if segments:
            predictions.append(Prediction(genome, mge_id, mge_type, location, tuple(segments)))
    return predictions, audit


def write_counts(rows: list[dict[str, str]], genome_ids: list[str], output: Path) -> None:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        counts[clean_genome(row.get("Isolate", ""))][element_type(row.get("MGE", ""))] += 1
    types = ["ICE", "IME", "rICE", "AICE", "Other"]
    out_rows = []
    for genome in genome_ids:
        values = [counts[genome][kind] for kind in types]
        out_rows.append([genome, *values, sum(values), int(sum(values) > 0)])
    write_tsv(output, ["genome_id", *types, "total", "has_ICE_or_IME"], out_rows)


def classify_existing_hotspots(hotspots: list[dict[str, str]], predictions: list[Prediction],
                               minimum_fraction: float, output: Path) -> None:
    all_regions = [region for prediction in predictions for region in prediction.segments]
    lookup = index_regions(all_regions)
    prediction_by_id = {prediction.mge_id: prediction for prediction in predictions}
    detail = []
    counts: Counter[tuple[str, str]] = Counter()
    for row in hotspots:
        start, end = int(row["start"]), int(row["end"])
        matches = region_matches(lookup, row["genome"], row["contig"], start, end, minimum_fraction)
        best = matches[0] if matches else None
        if best:
            mge_id = best[0].topology
            context = prediction_by_id[mge_id].mge_type
            fraction = best[1]
        else:
            mge_id, context, fraction = "", "non-ICE", 0.0
        counts[(row["hotspot_id"], context)] += 1
        detail.append([
            row["hotspot_id"], row["locus_id"], row["genome"], row["contig"], start, end,
            row.get("defence_system", ""), row.get("defence_present", ""), context,
            mge_id, f"{fraction:.6f}",
        ])
    write_tsv(output / "existing_hotspot_ice_context.tsv",
              ["hotspot_id", "locus_id", "genome", "contig", "start", "end",
               "defence_system", "defence_present", "ice_context", "matched_mge", "overlap_fraction"], detail)
    contexts = ["ICE", "IME", "rICE", "AICE", "Other", "non-ICE"]
    summary = []
    for hotspot in sorted({row["hotspot_id"] for row in hotspots}):
        values = [counts[(hotspot, context)] for context in contexts]
        summary.append([hotspot, *values, sum(values)])
    write_tsv(output / "existing_hotspot_ice_context_summary.tsv",
              ["hotspot_id", *contexts, "total"], summary)


def extract_single_contig_elements(predictions: list[Prediction], sources: dict[str, dict[str, Path]],
                                    systems: list[dict[str, str]], root: Path) -> tuple[Path, int, int]:
    element_root = root / "ice"
    element_root.mkdir(parents=True, exist_ok=True)
    manifest_rows, system_rows, requests, metadata_rows, skipped = [], [], [], [], []
    system_index: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in systems:
        system_index[(clean_genome(row["genome"]), clean_contig(row["contig"]))].append(row)
    eligible = [prediction for prediction in predictions if len(prediction.segments) == 1]
    for prediction in predictions:
        if len(prediction.segments) != 1:
            skipped.append([prediction.mge_id, prediction.genome, "fragmented multi-contig call",
                            len(prediction.segments)])
    by_genome: dict[str, list[Prediction]] = defaultdict(list)
    for prediction in eligible:
        by_genome[prediction.genome].append(prediction)
    counter = 0
    for genome, genome_predictions in sorted(by_genome.items()):
        files = sources.get(genome, {})
        if "gff" not in files or "faa" not in files:
            for prediction in genome_predictions:
                skipped.append([prediction.mge_id, genome, "missing source", "GFF and/or FAA"])
            continue
        proteins = read_fasta(files["faa"])
        genes, _ = read_gff(files["gff"], genome, proteins)
        genes_by_contig: dict[str, list[Gene]] = defaultdict(list)
        for gene in genes:
            genes_by_contig[clean_contig(gene.contig)].append(gene)
        for prediction in genome_predictions:
            region = prediction.segments[0]
            selected = [gene for gene in genes_by_contig[clean_contig(region.contig)]
                        if gene.start >= region.start and gene.end <= region.end and gene.sequence]
            if not selected:
                skipped.append([prediction.mge_id, genome, "annotation", "no fully contained proteins"])
                continue
            counter += 1
            local_id = f"ice_{counter:06d}__{safe_id(genome)}__{safe_id(prediction.mge_type)}"
            target = element_root / local_id
            target.mkdir(parents=True, exist_ok=True)
            local_genes = [replace(gene, genome=local_id, contig=local_id,
                                   start=gene.start-region.start+1, end=gene.end-region.start+1)
                           for gene in selected]
            faa = target / f"{local_id}.faa"
            gff = target / f"{local_id}.gff3"
            write_fasta(((gene.gene_id, gene.sequence) for gene in local_genes), faa)
            write_extracted_gff(gff, local_id, local_genes, region.length)
            manifest_rows.append([local_id, "", str(gff.resolve()), str(faa.resolve())])
            metadata_rows.append([local_id, prediction.mge_id, prediction.mge_type, genome,
                                  region.contig, region.start, region.end, prediction.location])
            for original, local in zip(selected, local_genes):
                requests.append([original.uid, local.uid])
            for system in system_index[(genome, clean_contig(region.contig))]:
                s1, s2 = int(system["start"]), int(system["end"])
                overlap = interval_overlap(region.start, region.end, s1, s2)
                if overlap / max(1, s2-s1+1) < 0.5:
                    continue
                system_rows.append([
                    system.get("locus_id", ""), local_id, local_id,
                    max(s1, region.start)-region.start+1,
                    min(s2, region.end)-region.start+1,
                    system["system"], system.get("callers", ""), system.get("source_ids", ""),
                ])
    write_tsv(element_root / "manifest.tsv", ["genome_id", "fna", "gff", "faa"], manifest_rows)
    write_tsv(element_root / "precomputed_systems.tsv",
              ["locus_id", "genome", "contig", "start", "end", "system", "callers", "source_ids"],
              system_rows)
    write_tsv(element_root / "cluster_requests.tsv", ["original_uid", "local_uid"], requests)
    write_tsv(element_root / "elements.tsv",
              ["element_id", "icefinder_mge", "mge_type", "source_genome", "source_contig",
               "source_start", "source_end", "icefinder_location"], metadata_rows)
    write_tsv(element_root / "skipped_elements.tsv",
              ["icefinder_mge", "genome", "reason", "detail"], skipped)
    return element_root, len(manifest_rows), len(system_rows)


def run_plots(script: Path, counts: Path, contexts: Path, outdir: Path) -> None:
    subprocess.run(["Rscript", str(script), str(counts), str(contexts), str(outdir)], check=True)


def map_hotspots_to_sources(hotspot_dir: Path, element_root: Path, denominator: int) -> None:
    hotspot_path = hotspot_dir / "hotspots.tsv"
    if not hotspot_path.is_file():
        return
    elements = {row["element_id"]: row for row in read_tsv(element_root / "elements.tsv")}
    mapped_rows = []
    pair_by_hotspot: dict[str, tuple[str, str]] = {}
    summaries: dict[str, dict[str, object]] = defaultdict(
        lambda: {"instances": 0, "defence": 0, "systems": Counter(), "types": Counter()})
    for row in read_tsv(hotspot_path):
        source = elements.get(row["genome"])
        if not source:
            continue
        source_start = int(source["source_start"])
        original_start = source_start + int(row["start"]) - 1
        original_end = source_start + int(row["end"]) - 1
        mapped_rows.append([
            row["hotspot_id"], row["locus_id"], row["genome"], source["icefinder_mge"],
            source["mge_type"], source["source_genome"], source["source_contig"],
            original_start, original_end, row["left_boundary_id"], row["right_boundary_id"],
            row["left_boundary_cluster"], row["right_boundary_cluster"],
            row["defence_system"], row["defence_present"], row["instances_of_architecture"],
        ])
        summary = summaries[row["hotspot_id"]]
        pair_by_hotspot[row["hotspot_id"]] = tuple(sorted(
            (row["left_boundary_cluster"], row["right_boundary_cluster"])))
        summary["instances"] += 1
        summary["types"][source["mge_type"]] += 1
        if row["defence_present"] == "1":
            summary["defence"] += 1
            for system in (value.strip() for value in row["defence_system"].split(";")):
                if system:
                    summary["systems"][system] += 1
    write_tsv(hotspot_dir.parent / "ice_hotspots_source_coordinates.tsv",
              ["hotspot_id", "locus_id", "element_id", "icefinder_mge", "mge_type",
               "source_genome", "source_contig", "source_start", "source_end",
               "left_boundary_id", "right_boundary_id", "left_boundary_cluster",
               "right_boundary_cluster", "defence_system", "defence_present",
               "instances_of_architecture"], mapped_rows)
    summary_rows = []
    for hotspot, values in sorted(summaries.items()):
        instances = int(values["instances"])
        defence = int(values["defence"])
        systems = values["systems"]
        types = values["types"]
        system_text = "; ".join(f"{name} ({count})" for name, count in systems.most_common())
        summary_rows.append([
            hotspot, instances, f"{100*instances/max(1, denominator):.2f}", defence,
            instances-defence, f"{100*defence/max(1, instances):.2f}",
            types["ICE"], types["IME"], types["rICE"], system_text,
        ])
    write_tsv(hotspot_dir.parent / "ice_hotspot_summary.tsv",
              ["hotspot_id", "instances", "percent_of_analyzable_elements", "with_defence",
               "boundary_only", "percent_with_defence", "ICE", "IME", "rICE",
               "defence_system_counts"], summary_rows)
    families: dict[tuple[str, str], list[str]] = defaultdict(list)
    for hotspot, pair in pair_by_hotspot.items():
        families[pair].append(hotspot)
    family_rows = []
    for index, (pair, hotspot_ids) in enumerate(sorted(families.items()), 1):
        instances = sum(int(summaries[h]["instances"]) for h in hotspot_ids)
        defence = sum(int(summaries[h]["defence"]) for h in hotspot_ids)
        systems: Counter[str] = Counter()
        types: Counter[str] = Counter()
        for hotspot in hotspot_ids:
            systems.update(summaries[hotspot]["systems"])
            types.update(summaries[hotspot]["types"])
        family_rows.append([
            f"ice_hotspot_family_{index:04d}", ";".join(sorted(hotspot_ids)), pair[0], pair[1],
            instances, f"{100*instances/max(1, denominator):.2f}", defence, instances-defence,
            f"{100*defence/max(1, instances):.2f}", types["ICE"], types["IME"], types["rICE"],
            "; ".join(f"{name} ({count})" for name, count in systems.most_common()),
        ])
    write_tsv(hotspot_dir.parent / "ice_hotspot_unique_boundary_families.tsv",
              ["boundary_family", "directional_hotspot_ids", "boundary_cluster_a",
               "boundary_cluster_b", "instances", "percent_of_analyzable_elements",
               "with_defence", "boundary_only", "percent_with_defence", "ICE", "IME", "rICE",
               "defence_system_counts"], family_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="PhorificICE: ICEfinder-opt analysis after PhorificV2")
    parser.add_argument("--phorific-results", type=Path, required=True)
    parser.add_argument("--icefinder-tsv", type=Path, required=True)
    parser.add_argument("--genome-files-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--analysis-only", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--min-hotspot-overlap", type=float, default=0.50)
    parser.add_argument("--window-kb", type=float, default=20.0)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-boundary-frequency", type=float, default=0.5)
    parser.add_argument("--min-comparable-loci", type=int, default=3)
    parser.add_argument("--min-hotspot-prevalence", type=float, default=0.05)
    parser.add_argument("--min-hotspot-members", type=int)
    parser.add_argument("--boundary-scan-span-factor", type=float, default=5.0)
    parser.add_argument("--include-edge-hotspots", action="store_true")
    parser.add_argument("--plot-palette", choices=("classic", "peach"), default="classic")
    args = parser.parse_args()
    if not (0 < args.min_hotspot_overlap <= 1 and 0 < args.min_hotspot_prevalence <= 1):
        fail("prevalence and overlap thresholds must lie between zero and one")
    results = args.phorific_results.resolve()
    output = args.outdir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    raw = read_tsv(args.icefinder_tsv.resolve())
    required = {"Isolate", "MGE", "Location"}
    if raw and not required.issubset(raw[0]):
        fail("ICEfinder table must contain Isolate, MGE, and Location columns")
    resolved = read_tsv(results / "resolved_inputs.tsv")
    genome_ids = [row["genome_id"] for row in resolved]
    genomes = set(genome_ids) | {clean_genome(row.get("Isolate", "")) for row in raw}
    sources = existing_resolved_sources(results / "resolved_inputs.tsv", genomes)
    discovered = discover_sources(args.genome_files_dir.resolve(), genomes)
    for genome, files in discovered.items():
        sources[genome].update({kind: path for kind, path in files.items() if kind not in sources[genome]})
    predictions, audit = parse_predictions(raw, sources)
    write_tsv(output / "icefinder_mapping_audit.tsv",
              ["raw_isolate", "genome_id", "mge_id", "mge_type", "location", "mapped",
               "segments", "problem"], audit)
    counts_path = output / "icefinder_counts_per_genome.tsv"
    write_counts(raw, genome_ids, counts_path)
    hotspots = read_tsv(results / "hotspots.tsv")
    classify_existing_hotspots(hotspots, predictions, args.min_hotspot_overlap, output)
    run_plots(Path(__file__).with_name("plot_icefinder_summary.R"), counts_path,
              output / "existing_hotspot_ice_context_summary.tsv", output)
    summary = {
        "genomes": len(genome_ids),
        "icefinder_predictions": len(raw),
        "mapped_predictions": len(predictions),
        "unmapped_predictions": len(raw)-len(predictions),
        "single_contig_predictions": sum(len(prediction.segments) == 1 for prediction in predictions),
        "fragmented_multicontig_predictions": sum(len(prediction.segments) > 1 for prediction in predictions),
        "counts_by_type": dict(Counter(element_type(row.get("MGE", "")) for row in raw)),
    }
    if args.analysis_only:
        (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary, indent=2))
        return
    systems = read_tsv(results / "systems.dereplicated.tsv")
    element_root, extracted, reused_systems = extract_single_contig_elements(
        predictions, sources, systems, output / "extracted")
    materialize_parent_clusters(results / "mmseqs" / "all_proteins_cluster.tsv", [element_root])
    old_pairs = output / "existing_chromosomal_boundary_pairs.tsv"
    write_old_boundary_pairs(hotspots, old_pairs)
    summary.update({"elements_extracted_for_hotspot_inference": extracted,
                    "reused_defence_calls": reused_systems})
    if not args.prepare_only:
        run_mge_hotspots("ice", element_root, output / "ice_hotspots", old_pairs, args)
        map_hotspots_to_sources(output / "ice_hotspots", element_root, extracted)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
