#!/usr/bin/env python3
"""PhorificV2 annotation with Pfam, AMRFinderPlus, and protein VFDB.

The workflow is split into stages so each external program can be run after
activating its own environment.  Raw results are reused unless --force is set.
The final table contains one row per gene and hotspot-locus context; repeated
hits from a tool are collapsed into unique semicolon-separated values.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

from pfam_priority import excluded_reason as pfam_excluded_reason
from pfam_priority import interesting_category


@dataclass(frozen=True)
class FastaRecord:
    name: str
    description: str
    sequence: str


def read_fasta(path: Path) -> Iterator[FastaRecord]:
    name = ""
    description = ""
    chunks: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="latin-1") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    yield FastaRecord(name, description, "".join(chunks))
                description = line[1:].strip()
                name = description.split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name:
        yield FastaRecord(name, description, "".join(chunks))


def write_fasta_record(handle, name: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{name}\n")
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def truthy(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def norm_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lstrip("#").lower())


def get_column(row: dict[str, str], *names: str, default: str = "") -> str:
    normalized = {norm_header(key): value for key, value in row.items()}
    for name in names:
        key = norm_header(name)
        if key in normalized:
            return normalized[key]
    return default


def unique_join(values: Iterable[str]) -> str:
    cleaned = {str(value).strip() for value in values if str(value).strip() not in {"", ".", "NA"}}
    return ";".join(sorted(cleaned, key=str.casefold))


def run_command(command: list[str], stdout_path: Path | None = None) -> None:
    print("Running:", " ".join(command), file=sys.stderr)
    if stdout_path is None:
        subprocess.run(command, check=True)
        return
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("w") as stdout:
        subprocess.run(command, check=True, stdout=stdout)


def require_executable(name: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"Required executable is not available in this environment: {name}")


def hotspot_rows(results_dir: Path) -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    path = results_dir / "hotspots.tsv"
    if not path.is_file():
        raise SystemExit(f"Missing hotspot table: {path}")
    rows = read_tsv(path)
    by_locus: dict[str, dict[str, str]] = {}
    for row in rows:
        locus_id = row["locus_id"]
        if locus_id in by_locus:
            raise SystemExit(f"Duplicate locus_id in hotspots.tsv: {locus_id}")
        by_locus[locus_id] = row
    return rows, by_locus


def gene_key(genome: str, contig: str, gene_id: str) -> str:
    return f"{genome}|{contig}|{gene_id}"


def prepare_proteins(results_dir: Path, outdir: Path) -> tuple[int, int]:
    rows, by_locus = hotspot_rows(results_dir)
    loci_dir = results_dir / "loci"
    if not loci_dir.is_dir():
        raise SystemExit(f"Missing final loci directory: {loci_dir}")

    proteins: dict[str, str] = {}
    gene_rows: list[dict[str, object]] = []

    for hotspot in rows:
        locus_id = hotspot["locus_id"]
        hotspot_id = hotspot["hotspot_id"]
        prefix = loci_dir / hotspot_id / locus_id
        faa = prefix.with_suffix(".faa")
        genes_tsv = prefix.with_suffix(".genes.tsv")
        if not faa.is_file() or not genes_tsv.is_file():
            raise SystemExit(f"Missing FAA or genes table for final hotspot locus: {prefix}")

        records = list(read_fasta(faa))
        record_by_gene = {record.name.split("|")[-1]: record for record in records}
        for record in records:
            prior = proteins.get(record.name)
            if prior is not None and prior != record.sequence:
                raise SystemExit(f"Protein identifier has inconsistent sequences: {record.name}")
            proteins[record.name] = record.sequence

        for gene in read_tsv(genes_tsv):
            genome = hotspot["genome"]
            contig = gene.get("contig", hotspot["contig"])
            gene_id = gene["gene_id"]
            record = record_by_gene.get(gene_id)
            key = record.name if record else gene_key(genome, contig, gene_id)
            gene_rows.append(
                {
                    "hotspot_id": hotspot_id,
                    "locus_id": locus_id,
                    "genome": genome,
                    "contig": contig,
                    "locus_start": hotspot["start"],
                    "locus_end": hotspot["end"],
                    "gene_id": gene_id,
                    "gene_key": key,
                    "start": gene.get("start", ""),
                    "end": gene.get("end", ""),
                    "strand": gene.get("strand", ""),
                    "product": gene.get("product", ""),
                    "is_defence": "1" if truthy(gene.get("is_defence", "0")) else "0",
                    "defence_system": hotspot.get("defence_system", "")
                    if truthy(gene.get("is_defence", "0"))
                    else "",
                    "protein_present": "1" if record else "0",
                }
            )

    outdir.mkdir(parents=True, exist_ok=True)
    protein_path = outdir / "hotspot_proteins.faa"
    with protein_path.open("w") as handle:
        for key in sorted(proteins):
            write_fasta_record(handle, key, proteins[key])

    gene_fields = [
        "hotspot_id",
        "locus_id",
        "genome",
        "contig",
        "locus_start",
        "locus_end",
        "gene_id",
        "gene_key",
        "start",
        "end",
        "strand",
        "product",
        "is_defence",
        "defence_system",
        "protein_present",
    ]
    write_tsv(outdir / "hotspot_gene_map.tsv", gene_rows, gene_fields)
    return len(proteins), len(gene_rows)


def fna_genome_id(path: Path) -> str:
    match = re.search(r"(GC[AF]_\d+\.\d+)", path.name)
    if match:
        return match.group(1)
    name = re.sub(r"\.(fna|fa|fasta)(\.gz)?$", "", path.name, flags=re.I)
    return re.sub(r"_genomic$", "", name)


def find_genome_fnas(fna_dir: Path) -> dict[str, Path]:
    candidates: dict[str, list[Path]] = defaultdict(list)
    for path in fna_dir.rglob("*"):
        if path.is_file() and re.search(r"\.(fna|fa|fasta)$", path.name, flags=re.I):
            candidates[fna_genome_id(path)].append(path)
    resolved: dict[str, Path] = {}
    for genome, paths in candidates.items():
        resolved[genome] = sorted(paths, key=lambda path: (len(path.parts), len(str(path))))[0]
        if len(paths) > 1:
            print(
                f"Warning: {genome} has {len(paths)} FNA candidates; using {resolved[genome]}",
                file=sys.stderr,
            )
    return resolved


def prepare_hotspot_dna(results_dir: Path, outdir: Path, fna_dir: Path) -> int:
    rows, _ = hotspot_rows(results_dir)
    fna_by_genome = find_genome_fnas(fna_dir)
    by_genome: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_genome[row["genome"]].append(row)

    missing = sorted(set(by_genome) - set(fna_by_genome))
    if missing:
        examples = ", ".join(missing[:12])
        raise SystemExit(f"No FNA found for {len(missing)} hotspot genomes; examples: {examples}")

    output = outdir / "hotspot_regions.fna"
    with output.open("w") as handle:
        for genome in sorted(by_genome):
            contigs = {record.name: record.sequence for record in read_fasta(fna_by_genome[genome])}
            for row in by_genome[genome]:
                contig = row["contig"]
                if contig not in contigs:
                    raise SystemExit(f"Contig {contig} for {row['locus_id']} not found in {fna_by_genome[genome]}")
                start = int(row["start"])
                end = int(row["end"])
                if start < 1 or end > len(contigs[contig]) or start > end:
                    raise SystemExit(
                        f"Invalid coordinates for {row['locus_id']}: {contig}:{start}-{end} "
                        f"(length {len(contigs[contig])})"
                    )
                write_fasta_record(handle, row["locus_id"], contigs[contig][start - 1 : end])
    return len(rows)


def stage_prepare(args: argparse.Namespace) -> None:
    results_dir = Path(args.results_dir).resolve()
    outdir = Path(args.outdir).resolve()
    unique_proteins, gene_contexts = prepare_proteins(results_dir, outdir)
    print(f"Prepared {unique_proteins} unique hotspot proteins across {gene_contexts} gene contexts")


def stage_pfam(args: argparse.Namespace) -> None:
    require_executable(args.hmmsearch)
    outdir = Path(args.outdir).resolve()
    proteins = outdir / "hotspot_proteins.faa"
    output = outdir / "pfam.domtblout"
    if output.is_file() and output.stat().st_size and not args.force:
        print(f"Reusing existing Pfam output: {output}")
        return
    if not proteins.is_file():
        raise SystemExit(f"Run the prepare stage first; missing {proteins}")
    pfam_hmm = Path(args.pfam_hmm).resolve()
    if not pfam_hmm.is_file():
        raise SystemExit(f"Pfam HMM database does not exist: {pfam_hmm}")
    run_command(
        [
            args.hmmsearch,
            "--cpu",
            str(args.threads),
            "--cut_ga",
            "--noali",
            "--domtblout",
            str(output),
            "-o",
            "/dev/null",
            str(pfam_hmm),
            str(proteins),
        ]
    )


def stage_amr(args: argparse.Namespace) -> None:
    require_executable(args.amrfinder)
    if args.recovery_batch_size < 1:
        raise SystemExit("--recovery-batch-size must be at least 1")
    outdir = Path(args.outdir).resolve()
    proteins = outdir / "hotspot_proteins.faa"
    output = outdir / "amrfinder.tsv"
    if output.is_file() and output.stat().st_size and not args.force:
        print(f"Reusing existing AMRFinderPlus output: {output}")
        return
    if not proteins.is_file():
        raise SystemExit(f"Run the prepare stage first; missing {proteins}")
    # A missing/corrupt AMRFinder database is a global configuration failure,
    # not evidence of a bad protein. Detect it before the crash-recovery logic
    # can recursively split the entire FASTA into pointless single records.
    database_check = subprocess.run(
        [args.amrfinder, "--database_version"], text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if database_check.returncode:
        detail = (database_check.stderr or database_check.stdout).strip()
        raise SystemExit(
            "AMRFinder database is unavailable. Run 'amrfinder -u' in the "
            f"AMRFinder environment before retrying.\n{detail}"
        )
    def command(input_faa: Path, result_tsv: Path) -> list[str]:
        return [
            args.amrfinder,
            "--protein",
            str(input_faa),
            "--threads",
            str(args.threads),
            "--output",
            str(result_tsv),
        ]

    try:
        run_command(command(proteins, output))
        (outdir / "amrfinder_failed_proteins.tsv").unlink(missing_ok=True)
        return
    except subprocess.CalledProcessError:
        output.unlink(missing_ok=True)
        if not args.recover_crashes:
            raise
        print(
            "AMRFinder failed on the combined FASTA; recursively isolating problematic proteins",
            file=sys.stderr,
        )

    records = list(read_fasta(proteins))
    successful_outputs: list[Path] = []
    failed_records: list[FastaRecord] = []

    with tempfile.TemporaryDirectory(prefix="amrfinder_recovery_", dir=outdir) as tmp:
        tmpdir = Path(tmp)
        counter = 0

        def run_batch(batch: list[FastaRecord]) -> None:
            nonlocal counter
            counter += 1
            batch_faa = tmpdir / f"batch_{counter:05d}.faa"
            batch_tsv = tmpdir / f"batch_{counter:05d}.tsv"
            with batch_faa.open("w") as handle:
                for record in batch:
                    write_fasta_record(handle, record.name, record.sequence)
            try:
                run_command(command(batch_faa, batch_tsv))
                successful_outputs.append(batch_tsv)
            except subprocess.CalledProcessError:
                batch_tsv.unlink(missing_ok=True)
                if len(batch) == 1:
                    failed_records.append(batch[0])
                    print(f"Isolated AMRFinder-crashing protein: {batch[0].name}", file=sys.stderr)
                    return
                midpoint = len(batch) // 2
                run_batch(batch[:midpoint])
                run_batch(batch[midpoint:])

        # Start with bounded chunks so one report-crashing protein does not
        # force AMRFinder to repeat half of a very large proteome at every
        # recursion level. Only a failing chunk is subdivided further.
        for start in range(0, len(records), args.recovery_batch_size):
            run_batch(records[start:start + args.recovery_batch_size])

        header: str | None = None
        data_lines: list[str] = []
        for result in successful_outputs:
            lines = [line for line in result.read_text().splitlines() if line]
            if not lines:
                continue
            if header is None:
                header = lines[0]
            elif lines[0] != header:
                raise SystemExit(f"AMRFinder batch headers do not match: {result}")
            data_lines.extend(lines[1:])
        if header is None:
            raise SystemExit("Every AMRFinder recovery batch failed; no result table could be created")
        output.write_text("\n".join([header, *data_lines]) + "\n")

    write_tsv(
        outdir / "amrfinder_failed_proteins.tsv",
        ({"protein_id": record.name, "sequence_length": len(record.sequence)} for record in failed_records),
        ["protein_id", "sequence_length"],
    )
    print(
        f"AMRFinder recovery completed with {len(failed_records)} isolated crashing protein(s); "
        f"all other results were combined in {output}",
        file=sys.stderr,
    )


def stage_virulence(args: argparse.Namespace) -> None:
    """Search hotspot proteins against selected VFDB protein categories."""
    require_executable(args.diamond)
    outdir = Path(args.outdir).resolve()
    proteins = outdir / "hotspot_proteins.faa"
    vfdb = Path(args.vfdb_proteins).resolve()
    selected_faa = outdir / "vfdb_selected_proteins.faa"
    database = outdir / "vfdb_selected"
    output = outdir / "diamond_vfdb.tsv"
    if output.is_file() and output.stat().st_size and not args.force:
        print(f"Reusing existing DIAMOND/VFDB output: {output}")
        return
    if not proteins.is_file():
        raise SystemExit(f"Run the prepare stage first; missing {proteins}")
    if not vfdb.is_file():
        raise SystemExit(f"VFDB protein FASTA does not exist: {vfdb}")

    selected_count = 0
    with selected_faa.open("w") as handle:
        for record in read_fasta(vfdb):
            if not classify_virulence("", record.description, ""):
                continue
            handle.write(f">{record.description}\n")
            for start in range(0, len(record.sequence), 80):
                handle.write(record.sequence[start : start + 80] + "\n")
            selected_count += 1
    if not selected_count:
        raise SystemExit("No requested adhesin, toxin, or T7SS records were recognised in VFDB")

    run_command([args.diamond, "makedb", "--in", str(selected_faa), "--db", str(database)])
    run_command(
        [
            args.diamond,
            "blastp",
            "--query",
            str(proteins),
            "--db",
            str(database),
            "--out",
            str(output),
            "--threads",
            str(args.threads),
            "--id",
            str(args.min_identity),
            "--query-cover",
            str(args.min_coverage),
            "--subject-cover",
            str(args.min_coverage),
            "--evalue",
            str(args.evalue),
            "--max-target-seqs",
            str(args.max_targets),
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "qlen",
            "slen",
            "qcovhsp",
            "scovhsp",
            "evalue",
            "bitscore",
            "stitle",
        ],
    )


def stage_plots(args: argparse.Namespace) -> None:
    try:
        from plot_hotspot_amr import make_plots
        from plot_hotspot_pfam import make_pfam_plots
        from plot_hotspot_virulence import make_virulence_plot
    except ImportError as error:
        raise SystemExit(
            "The plots stage needs matplotlib, numpy, and seaborn. "
            "Run it from a Python/Conda environment containing those packages."
        ) from error
    results_dir = Path(args.results_dir).resolve()
    outdir = Path(args.outdir).resolve()
    output_prefix = Path(args.output_prefix).resolve() if args.output_prefix else outdir / "hotspot_amr_summary"
    make_plots(results_dir, outdir, output_prefix)
    pfam_output_prefix = outdir / "hotspot_pfam_summary"
    if (outdir / "pfam.domtblout").is_file():
        make_pfam_plots(results_dir, outdir, pfam_output_prefix, args.top_pfam_domains)
    else:
        print(f"Skipping Pfam plots because {outdir / 'pfam.domtblout'} is missing")
    if (outdir / "diamond_vfdb_selected.tsv").is_file():
        make_virulence_plot(results_dir, outdir)
    else:
        print(f"Skipping virulence plots because {outdir / 'diamond_vfdb_selected.tsv'} is missing")


VIRULENCE_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    (
        "T7SS",
        re.compile(
            r"type\s*(?:vii|7)\s*secretion|\bt7ss\b|\bwxg100\b|"
            r"\b(?:esx[a-z0-9]*|ess[abc]|esa[abc]|ecc[a-e][0-9]*|mycp[0-9]*)\b",
            re.I,
        ),
    ),
    (
        "endotoxin",
        re.compile(r"\bendotoxins?\b|lipopolysaccharide|\blipid\s*a\b|\blpx[a-z]\b", re.I),
    ),
    (
        "exotoxin",
        re.compile(
            r"\bexotoxins?\b|\bcytolysins?\b|\bhemolysins?\b|\bleukocidins?\b|"
            r"\benterotoxins?\b|pore[- ]forming toxin|secreted toxin|\bcyl[ablm]\b|\bhly[a-z]?\b",
            re.I,
        ),
    ),
    (
        "adhesin",
        re.compile(
            r"\badhes(?:in|ion|ive|ion factor|ion protein)\b|\badherence\b|"
            r"aggregation substance|\bmscramm\b|\bfimbr(?:ia|ial|iae)\b|"
            r"\bpilus\b|\bpili\b|\b(?:acm|scm|esp|ace|efaa|ebp[abc]|fss[123])\b",
            re.I,
        ),
    ),
]


def classify_virulence(gene: str, product: str, accession: str = "") -> str:
    text = " ".join((gene, product, accession))
    categories = [name for name, pattern in VIRULENCE_PATTERNS if pattern.search(text)]
    return unique_join(categories)


def parse_pfam(path: Path) -> dict[str, list[dict[str, str]]]:
    hits: dict[str, list[dict[str, str]]] = defaultdict(list)
    seen: set[tuple[str, ...]] = set()
    if not path.is_file():
        return hits
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split(maxsplit=22)
            if len(parts) < 22:
                continue
            protein = parts[0]
            pfam_name = parts[3]
            accession = parts[4]
            item = {
                "accession": accession,
                "name": pfam_name,
                "description": parts[22] if len(parts) > 22 else "",
                "i_evalue": parts[12],
                "score": parts[13],
                "ali_from": parts[17],
                "ali_to": parts[18],
            }
            signature = (protein, accession, parts[17], parts[18])
            if signature not in seen:
                hits[protein].append(item)
                seen.add(signature)
    return hits


def parse_amr(path: Path) -> dict[str, list[dict[str, str]]]:
    hits: dict[str, list[dict[str, str]]] = defaultdict(list)
    if not path.is_file() or not path.stat().st_size:
        return hits
    for row in read_tsv(path):
        protein = get_column(row, "Protein identifier", "Protein id")
        if not protein:
            continue
        hits[protein].append(
            {
                "gene": get_column(row, "Gene symbol", "Element symbol"),
                "sequence_name": get_column(
                    row, "Sequence name", "Name of closest sequence", "Element name"
                ),
                "element_type": get_column(row, "Element type", "Type"),
                "element_subtype": get_column(row, "Element subtype", "Subtype"),
                "class": get_column(row, "Class"),
                "subclass": get_column(row, "Subclass"),
                "method": get_column(row, "Method"),
                "identity": get_column(
                    row, "% Identity to reference sequence", "% Identity to reference"
                ),
                "coverage": get_column(
                    row, "% Coverage of reference sequence", "% Coverage of reference"
                ),
                "accession": get_column(
                    row, "Accession of closest sequence", "Closest reference accession"
                ),
            }
        )
    return hits


def overlap_size(start1: int, end1: int, start2: int, end2: int) -> int:
    return max(0, min(end1, end2) - max(start1, start2) + 1)


def parse_selected_virulence(
    path: Path, gene_rows: list[dict[str, str]]
) -> tuple[dict[tuple[str, str], list[dict[str, str]]], list[dict[str, str]], int]:
    selected: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    selected_rows: list[dict[str, str]] = []
    rejected = 0
    if not path.is_file() or not path.stat().st_size:
        return selected, selected_rows, rejected

    genes_by_locus: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in gene_rows:
        genes_by_locus[row["locus_id"]].append(row)

    for row in read_tsv(path):
        locus_id = get_column(row, "SEQUENCE")
        gene = get_column(row, "GENE")
        product = get_column(row, "PRODUCT")
        accession = get_column(row, "ACCESSION")
        category = classify_virulence(gene, product, accession)
        if not category:
            rejected += 1
            continue
        if locus_id not in genes_by_locus:
            continue

        local_start = int(float(get_column(row, "START")))
        local_end = int(float(get_column(row, "END")))
        locus_start = int(genes_by_locus[locus_id][0]["locus_start"])
        genomic_start = locus_start + min(local_start, local_end) - 1
        genomic_end = locus_start + max(local_start, local_end) - 1

        best_gene: dict[str, str] | None = None
        best_overlap = 0
        for candidate in genes_by_locus[locus_id]:
            if not candidate["start"] or not candidate["end"]:
                continue
            overlap = overlap_size(
                genomic_start,
                genomic_end,
                int(candidate["start"]),
                int(candidate["end"]),
            )
            if overlap > best_overlap:
                best_overlap = overlap
                best_gene = candidate
        if best_gene is None or best_overlap == 0:
            continue

        item = {
            "category": category,
            "gene": gene,
            "product": product,
            "accession": accession,
            "identity": get_column(row, "%IDENTITY"),
            "coverage": get_column(row, "%COVERAGE"),
            "database": get_column(row, "DATABASE"),
            "hit_start": str(genomic_start),
            "hit_end": str(genomic_end),
            "gene_overlap_bp": str(best_overlap),
        }
        selected[(locus_id, best_gene["gene_key"])].append(item)
        selected_rows.append(
            {
                "hotspot_id": best_gene["hotspot_id"],
                "locus_id": locus_id,
                "genome": best_gene["genome"],
                "contig": best_gene["contig"],
                "gene_id": best_gene["gene_id"],
                "gene_key": best_gene["gene_key"],
                **item,
            }
        )
    return selected, selected_rows, rejected


def parse_protein_virulence(
    path: Path, gene_rows: list[dict[str, str]]
) -> tuple[dict[tuple[str, str], list[dict[str, str]]], list[dict[str, str]], int]:
    """Parse DIAMOND VFDB matches and expand unique proteins to locus contexts."""
    selected: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    selected_rows: list[dict[str, str]] = []
    rejected = 0
    if not path.is_file() or not path.stat().st_size:
        return selected, selected_rows, rejected

    contexts: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in gene_rows:
        contexts[row["gene_key"]].append(row)

    best: dict[tuple[str, str], dict[str, str]] = {}
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t", 10)
            if len(parts) != 11:
                continue
            qseqid, sseqid, identity, length, qlen, slen, qcov, scov, evalue, bitscore, title = parts
            category_string = classify_virulence("", title, sseqid)
            if not category_string:
                rejected += 1
                continue
            accession_match = re.search(r"(?:gb|ref)\|([^)|\s]+)", title)
            accession = accession_match.group(1) if accession_match else ""
            gene_match = re.search(r"\s\(([^()]+)\)\s", title)
            gene = gene_match.group(1) if gene_match else sseqid.split("(", 1)[0]
            product_match = re.search(r"\)\s+([^[]+?)\s*\[", title)
            product = product_match.group(1).strip() if product_match else title
            for category in category_string.split(";"):
                item = {
                    "category": category,
                    "gene": gene,
                    "product": product,
                    "accession": accession,
                    "vfdb_id": sseqid.split("(", 1)[0],
                    "identity": identity,
                    "query_coverage": qcov,
                    "reference_coverage": scov,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "database": "VFDB_setA_pro",
                }
                signature = (qseqid, category)
                if signature not in best or float(bitscore) > float(best[signature]["bitscore"]):
                    best[signature] = item

    for (gene_key_value, _category), item in best.items():
        for context in contexts.get(gene_key_value, []):
            selected[(context["locus_id"], gene_key_value)].append(item)
            selected_rows.append(
                {
                    "hotspot_id": context["hotspot_id"],
                    "locus_id": context["locus_id"],
                    "genome": context["genome"],
                    "contig": context["contig"],
                    "gene_id": context["gene_id"],
                    "gene_key": gene_key_value,
                    **item,
                }
            )
    return selected, selected_rows, rejected


def stage_combine(args: argparse.Namespace) -> None:
    results_dir = Path(args.results_dir).resolve()
    outdir = Path(args.outdir).resolve()
    gene_map = outdir / "hotspot_gene_map.tsv"
    if not gene_map.is_file():
        prepare_proteins(results_dir, outdir)
    genes = read_tsv(gene_map)

    required_raw = [outdir / "pfam.domtblout", outdir / "amrfinder.tsv", outdir / "diamond_vfdb.tsv"]
    missing_raw = [str(path) for path in required_raw if not path.is_file()]
    if missing_raw:
        raise SystemExit(
            "Cannot combine before all three searches have completed. Missing: " + ", ".join(missing_raw)
        )

    pfam = parse_pfam(outdir / "pfam.domtblout")
    amr = parse_amr(outdir / "amrfinder.tsv")
    virulence, selected_virulence_rows, rejected_virulence = parse_protein_virulence(
        outdir / "diamond_vfdb.tsv", genes
    )

    combined: list[dict[str, object]] = []
    excluded_pfam_rows: list[dict[str, object]] = []
    for gene in genes:
        key = gene["gene_key"]
        pfam_hits = pfam.get(key, [])
        amr_hits = amr.get(key, [])
        is_defence = truthy(gene["is_defence"])
        exclusion_reasons = []
        if is_defence:
            exclusion_reasons.append("defence")
        if amr_hits:
            exclusion_reasons.append("AMR")
        if exclusion_reasons and pfam_hits:
            for hit in pfam_hits:
                excluded_pfam_rows.append(
                    {
                        "hotspot_id": gene["hotspot_id"],
                        "locus_id": gene["locus_id"],
                        "genome": gene["genome"],
                        "gene_id": gene["gene_id"],
                        "gene_key": key,
                        "defence_system": gene["defence_system"],
                        "amr_gene": unique_join(amr_hit["gene"] for amr_hit in amr_hits),
                        "exclusion_reason": unique_join(exclusion_reasons),
                        "pfam_accession": hit["accession"],
                        "pfam_name": hit["name"],
                        "pfam_description": hit["description"],
                    }
                )
            pfam_hits = []
        elif pfam_hits:
            prioritised_hits = []
            for hit in pfam_hits:
                reason = pfam_excluded_reason(hit)
                category = interesting_category(hit) if not reason else ""
                if reason or not category:
                    excluded_pfam_rows.append(
                        {
                            "hotspot_id": gene["hotspot_id"],
                            "locus_id": gene["locus_id"],
                            "genome": gene["genome"],
                            "gene_id": gene["gene_id"],
                            "gene_key": key,
                            "defence_system": gene["defence_system"],
                            "amr_gene": "",
                            "exclusion_reason": reason or "not_prioritised_accessory_function",
                            "pfam_accession": hit["accession"],
                            "pfam_name": hit["name"],
                            "pfam_description": hit["description"],
                        }
                    )
                    continue
                prioritised_hits.append({**hit, "candidate_category": category})
            pfam_hits = prioritised_hits

        virulence_hits = virulence.get((gene["locus_id"], key), [])
        combined.append(
            {
                **gene,
                "amr_gene": unique_join(hit["gene"] for hit in amr_hits),
                "amr_sequence_name": unique_join(hit["sequence_name"] for hit in amr_hits),
                "amr_element_type": unique_join(hit["element_type"] for hit in amr_hits),
                "amr_element_subtype": unique_join(hit["element_subtype"] for hit in amr_hits),
                "amr_class": unique_join(hit["class"] for hit in amr_hits),
                "amr_subclass": unique_join(hit["subclass"] for hit in amr_hits),
                "amr_method": unique_join(hit["method"] for hit in amr_hits),
                "amr_accession": unique_join(hit["accession"] for hit in amr_hits),
                "virulence_category": unique_join(hit["category"] for hit in virulence_hits),
                "virulence_gene": unique_join(hit["gene"] for hit in virulence_hits),
                "virulence_product": unique_join(hit["product"] for hit in virulence_hits),
                "virulence_accession": unique_join(hit["accession"] for hit in virulence_hits),
                "virulence_identity": unique_join(hit["identity"] for hit in virulence_hits),
                "virulence_query_coverage": unique_join(hit["query_coverage"] for hit in virulence_hits),
                "virulence_reference_coverage": unique_join(
                    hit["reference_coverage"] for hit in virulence_hits
                ),
                "pfam_accession": unique_join(hit["accession"] for hit in pfam_hits),
                "pfam_name": unique_join(hit["name"] for hit in pfam_hits),
                "pfam_description": unique_join(hit["description"] for hit in pfam_hits),
                "pfam_candidate_category": unique_join(
                    hit["candidate_category"] for hit in pfam_hits
                ),
            }
        )

    combined_fields = [
        "hotspot_id",
        "locus_id",
        "genome",
        "contig",
        "gene_id",
        "gene_key",
        "start",
        "end",
        "strand",
        "product",
        "is_defence",
        "defence_system",
        "amr_gene",
        "amr_sequence_name",
        "amr_element_type",
        "amr_element_subtype",
        "amr_class",
        "amr_subclass",
        "amr_method",
        "amr_accession",
        "virulence_category",
        "virulence_gene",
        "virulence_product",
        "virulence_accession",
        "virulence_identity",
        "virulence_query_coverage",
        "virulence_reference_coverage",
        "pfam_accession",
        "pfam_name",
        "pfam_description",
        "pfam_candidate_category",
    ]
    write_tsv(outdir / "hotspot_annotations.tsv", combined, combined_fields)

    selected_fields = [
        "hotspot_id",
        "locus_id",
        "genome",
        "contig",
        "gene_id",
        "gene_key",
        "category",
        "gene",
        "product",
        "accession",
        "identity",
        "vfdb_id",
        "query_coverage",
        "reference_coverage",
        "evalue",
        "bitscore",
        "database",
    ]
    write_tsv(outdir / "diamond_vfdb_selected.tsv", selected_virulence_rows, selected_fields)
    excluded_fields = [
        "hotspot_id",
        "locus_id",
        "genome",
        "gene_id",
        "gene_key",
        "defence_system",
        "amr_gene",
        "exclusion_reason",
        "pfam_accession",
        "pfam_name",
        "pfam_description",
    ]
    write_tsv(outdir / "pfam_hits_excluded.tsv", excluded_pfam_rows, excluded_fields)
    defence_excluded_rows = [
        row for row in excluded_pfam_rows if "defence" in str(row["exclusion_reason"]).split(";")
    ]
    write_tsv(
        outdir / "pfam_hits_excluded_for_defence_genes.tsv",
        defence_excluded_rows,
        excluded_fields,
    )

    summary_lines = [
        f"Hotspot gene contexts: {len(combined)}",
        f"Unique protein identifiers: {len({row['gene_key'] for row in genes if row['protein_present'] == '1'})}",
        f"Defence gene contexts: {sum(truthy(row['is_defence']) for row in genes)}",
        f"Gene contexts with AMRFinder hits: {sum(bool(row['amr_gene']) for row in combined)}",
        f"Gene contexts with selected VFDB hits: {sum(bool(row['virulence_gene']) for row in combined)}",
        f"Selected VFDB hit rows: {len(selected_virulence_rows)}",
        f"DIAMOND hits rejected by requested category filter: {rejected_virulence}",
        f"Gene contexts with retained Pfam hits: {sum(bool(row['pfam_accession']) for row in combined)}",
        f"Pfam hit rows excluded because the gene is defence: {len(defence_excluded_rows)}",
        f"Pfam hit rows excluded because the gene has an AMRFinder hit: "
        f"{sum('AMR' in str(row['exclusion_reason']).split(';') for row in excluded_pfam_rows)}",
        f"Pfam hit rows excluded as integration/transposition: "
        f"{sum(row['exclusion_reason'] == 'integration_or_transposition' for row in excluded_pfam_rows)}",
        f"Pfam hit rows excluded as generic regulation/DNA binding: "
        f"{sum(row['exclusion_reason'] == 'generic_regulation_or_dna_binding' for row in excluded_pfam_rows)}",
        f"Pfam hit rows excluded as phage structure/lysis: "
        f"{sum(row['exclusion_reason'] == 'phage_structure_or_lysis' for row in excluded_pfam_rows)}",
        f"Pfam hit rows excluded as unprioritised background function: "
        f"{sum(row['exclusion_reason'] == 'not_prioritised_accessory_function' for row in excluded_pfam_rows)}",
    ]
    (outdir / "annotation_summary.txt").write_text("\n".join(summary_lines) + "\n")
    print("\n".join(summary_lines))
    print(f"Wrote combined annotations to {outdir / 'hotspot_annotations.tsv'}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="PhorificV2 hotspot annotation: Pfam, AMRFinderPlus, and DIAMOND/VFDB"
    )
    subparsers = parser.add_subparsers(dest="stage", required=True)

    prepare = subparsers.add_parser("prepare", help="Build the unique hotspot-protein FASTA")
    prepare.add_argument("--results-dir", required=True)
    prepare.add_argument("--outdir", required=True)
    prepare.set_defaults(function=stage_prepare)

    pfam = subparsers.add_parser("pfam", help="Run hmmsearch against Pfam-A")
    pfam.add_argument("--outdir", required=True)
    pfam.add_argument("--pfam-hmm", required=True)
    pfam.add_argument("--threads", type=int, default=8)
    pfam.add_argument("--hmmsearch", default="hmmsearch")
    pfam.add_argument("--force", action="store_true")
    pfam.set_defaults(function=stage_pfam)

    amr = subparsers.add_parser("amr", help="Run AMRFinderPlus on hotspot proteins")
    amr.add_argument("--outdir", required=True)
    amr.add_argument("--threads", type=int, default=8)
    amr.add_argument("--amrfinder", default="amrfinder")
    amr.add_argument("--force", action="store_true")
    amr.add_argument(
        "--recovery-batch-size", type=int, default=1000,
        help="initial AMRFinder crash-recovery batch size (default: 1000 proteins)",
    )
    amr.add_argument(
        "--no-recover-crashes",
        dest="recover_crashes",
        action="store_false",
        help="Fail immediately instead of isolating proteins that trigger an AMRFinder report crash",
    )
    amr.set_defaults(recover_crashes=True)
    amr.set_defaults(function=stage_amr)

    virulence = subparsers.add_parser("virulence", help="Run DIAMOND BLASTP against VFDB proteins")
    virulence.add_argument("--outdir", required=True)
    virulence.add_argument("--vfdb-proteins", required=True)
    virulence.add_argument("--threads", type=int, default=8)
    virulence.add_argument("--diamond", default="diamond")
    virulence.add_argument("--min-identity", type=float, default=35.0)
    virulence.add_argument("--min-coverage", type=float, default=70.0)
    virulence.add_argument("--evalue", type=float, default=1e-10)
    virulence.add_argument("--max-targets", type=int, default=25)
    virulence.add_argument("--force", action="store_true")
    virulence.set_defaults(function=stage_virulence)

    combine = subparsers.add_parser("combine", help="Filter, map, deduplicate, and combine all annotations")
    combine.add_argument("--results-dir", required=True)
    combine.add_argument("--outdir", required=True)
    combine.set_defaults(function=stage_combine)

    plots = subparsers.add_parser("plots", help="Plot annotation frequencies by hotspot")
    plots.add_argument("--results-dir", required=True)
    plots.add_argument("--outdir", required=True)
    plots.add_argument("--output-prefix", help="Optional output prefix for PDF, PNG, and data TSV")
    plots.add_argument("--top-pfam-domains", type=int, default=15)
    plots.set_defaults(function=stage_plots)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if hasattr(args, "threads") and args.threads < 1:
        parser.error("--threads must be at least 1")
    if hasattr(args, "top_pfam_domains") and args.top_pfam_domains < 1:
        parser.error("--top-pfam-domains must be at least 1")
    args.function(args)


if __name__ == "__main__":
    main()
