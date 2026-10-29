#!/usr/bin/env python3
"""PhorificV2: find recurrent defence hotspots from paired genome annotations.

The pipeline deliberately uses only the Python standard library. PADLOC and
MMseqs2 are invoked as external programs, and all intermediate tables remain
available for auditing and threshold tuning.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path
from urllib.parse import unquote


@dataclass
class Gene:
    genome: str
    contig: str
    start: int
    end: int
    strand: str
    gene_id: str
    product: str = "hypothetical protein"
    sequence: str = ""
    cluster: str = ""
    defence: bool = False

    @property
    def uid(self) -> str:
        return f"{self.genome}|{self.contig}|{self.gene_id}"


@dataclass
class SystemCall:
    genome: str
    route: str
    source_id: str
    contig: str
    start: int
    end: int
    system: str
    proteins: set[str] = field(default_factory=set)
    routes: set[str] = field(default_factory=set)


def fail(message: str) -> None:
    raise SystemExit(f"error: {message}")


def run(cmd: list[str], log: Path, dry_run: bool = False) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a") as handle:
        handle.write("$ " + " ".join(cmd) + "\n")
        if dry_run:
            return
        proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT)
    if proc.returncode:
        fail(f"command failed ({proc.returncode}); see {log}")


def env_command(runner: str, environment: str | None, executable: str, arguments: list[str]) -> list[str]:
    """Build a command in the selected Conda environment, or use the current PATH."""
    if not environment:
        return [executable, *arguments]
    selector = "-p" if "/" in environment else "-n"
    return [runner, "run", "--no-capture-output", selector, environment, executable, *arguments]


def require_tool(runner: str, environment: str | None, executable: str, label: str) -> None:
    if environment:
        if shutil.which(runner) is None:
            fail(f"environment runner not found: {runner}")
    elif shutil.which(executable) is None:
        fail(f"{label} executable not found: {executable}")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = ""
    chunks: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name:
        records[name] = "".join(chunks)
    return records


def write_fasta(records, path: Path, width: int = 80) -> None:
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for i in range(0, len(sequence), width):
                handle.write(sequence[i:i + width] + "\n")


def attrs(text: str) -> dict[str, str]:
    result = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = unquote(value)
    return result


def read_gff(path: Path, genome: str, proteins: dict[str, str]) -> tuple[list[Gene], dict[str, int]]:
    genes, lengths = [], {}
    with path.open() as handle:
        for raw in handle:
            if raw.startswith("##sequence-region"):
                parts = raw.split()
                if len(parts) >= 4:
                    lengths[parts[1]] = int(parts[3])
            if raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2].lower() not in {"cds", "gene"}:
                continue
            a = attrs(cols[8])
            gene_id = a.get("protein_id") or a.get("locus_tag") or a.get("ID")
            if not gene_id:
                continue
            candidates = [gene_id, re.sub(r"^cds-", "", gene_id), a.get("ID", ""), a.get("locus_tag", "")]
            seq = next((proteins[x] for x in candidates if x in proteins), "")
            product = a.get("product") or a.get("Name") or "hypothetical protein"
            genes.append(Gene(genome, cols[0], int(cols[3]), int(cols[4]), cols[6], gene_id, product, seq))
            lengths[cols[0]] = max(lengths.get(cols[0], 0), int(cols[4]))
    # Prefer CDS if a GFF contains overlapping gene and CDS rows.
    unique = {(g.contig, g.start, g.end, g.strand): g for g in genes}
    return sorted(unique.values(), key=lambda g: (g.contig, g.start, g.end)), lengths


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"genome_id", "gff", "faa"}
    if not rows or not required.issubset(rows[0]):
        fail("manifest must be tab-separated with genome_id, gff, faa columns; fna is optional")
    base = path.parent
    for row in rows:
        for key in ("gff", "faa"):
            p = Path(row[key])
            row[key] = str(p if p.is_absolute() else (base / p).resolve())
            if not Path(row[key]).is_file():
                fail(f"missing {key} for {row['genome_id']}: {row[key]}")
        if row.get("fna"):
            p = Path(row["fna"])
            row["fna"] = str(p if p.is_absolute() else (base / p).resolve())
            if not Path(row["fna"]).is_file():
                fail(f"missing optional fna for {row['genome_id']}: {row['fna']}")
        else:
            row["fna"] = ""
    return rows


def discover_genomes(directory: Path) -> list[dict[str, str]]:
    """Recursively discover same-stem GFF/FAA pairs and optional FNA files."""
    directory = directory.resolve()
    if not directory.is_dir():
        fail(f"input directory does not exist: {directory}")
    annotations = sorted(set(directory.rglob("*.gff3")) | set(directory.rglob("*.gff")))
    rows, incomplete = [], []
    for gff in annotations:
        stem = gff.name[:-5] if gff.name.endswith(".gff3") else gff.name[:-4]
        faa, fna = gff.with_name(stem + ".faa"), gff.with_name(stem + ".fna")
        if not faa.is_file():
            incomplete.append(f"{gff}: missing .faa")
            continue
        rows.append({"genome_id": stem, "fna": str(fna) if fna.is_file() else "",
                     "gff": str(gff), "faa": str(faa)})
    if not rows:
        detail = "\n".join(incomplete[:10])
        fail(f"no complete same-stem GFF/FAA genome pairs found under {directory}"
             + (f"\n{detail}" if detail else ""))
    duplicates = [name for name, count in Counter(r["genome_id"] for r in rows).items() if count > 1]
    if duplicates:
        fail("duplicate genome basenames found in different folders: " + ", ".join(sorted(duplicates))
             + "; rename them so every genome basename is unique")
    if incomplete:
        print(f"warning: ignored {len(incomplete)} incomplete annotation set(s)", file=sys.stderr)
        for message in incomplete[:10]:
            print(f"warning: {message}", file=sys.stderr)
    return rows


def padloc_csv(directory: Path) -> Path | None:
    hits = sorted(directory.glob("*_padloc.csv"))
    if not hits:
        return None
    if len(hits) != 1:
        fail(f"expected one *_padloc.csv in {directory}, found {len(hits)}")
    return hits[0]


def padloc_run_complete(directory: Path, log: Path) -> bool:
    if padloc_csv(directory) is not None or (directory / ".hotspot_pipeline.complete").is_file():
        return True
    # PADLOC writes no CSV when a successful run finds no systems.
    return (any(directory.glob("*.domtblout")) and log.is_file()
            and "Nothing found" in log.read_text(errors="replace"))


def parse_padloc(path: Path, genome: str, route: str) -> list[SystemCall]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            number = row.get("system.number") or row.get("system_number") or "1"
            groups[(row.get("seqid", ""), number)].append(row)
    calls = []
    for (contig, number), rows in groups.items():
        starts = [int(float(r["start"])) for r in rows if r.get("start")]
        ends = [int(float(r["end"])) for r in rows if r.get("end")]
        if not contig or not starts or not ends:
            continue
        names = [r.get("system", "unknown") for r in rows]
        proteins = {r.get("target.name", "") for r in rows if r.get("target.name")}
        calls.append(SystemCall(genome, route, f"{route}:{number}", contig, min(starts), max(ends),
                                Counter(names).most_common(1)[0][0], proteins, {route}))
    return calls


def defensefinder_systems(directory: Path) -> Path | None:
    hits = sorted(directory.glob("*defense_finder_systems.tsv"))
    if not hits:
        return None
    if len(hits) != 1:
        fail(f"expected one *defense_finder_systems.tsv in {directory}, found {len(hits)}")
    return hits[0]


def make_defensefinder_faa(genes: list[Gene], path: Path) -> dict[str, Gene]:
    """Write ordered gembase protein IDs and retain an exact map to GFF genes."""
    mapping: dict[str, Gene] = {}
    records = []
    by_contig: dict[str, list[Gene]] = defaultdict(list)
    for gene in genes:
        if gene.sequence:
            by_contig[gene.contig].append(gene)
    for contig in sorted(by_contig):
        safe_contig = re.sub(r"[^A-Za-z0-9.-]", "-", contig)
        for position, gene in enumerate(sorted(by_contig[contig], key=lambda x: (x.start, x.end)), 1):
            df_id = f"{gene.genome}--{safe_contig}_{position}"
            mapping[df_id] = gene
            records.append((df_id, gene.sequence))
    write_fasta(records, path)
    return mapping


def split_ids(value: str) -> list[str]:
    return [x.strip().strip("[]'\"") for x in re.split(r"[,;\s]+", value or "")
            if x.strip().strip("[]'\"")]


def parse_defensefinder(path: Path, genome: str, mapping: dict[str, Gene]) -> list[SystemCall]:
    calls = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            ids = split_ids(row.get("protein_in_syst", ""))
            if not ids:
                ids = [x for x in (row.get("sys_beg", ""), row.get("sys_end", "")) if x]
            genes = [mapping[x] for x in ids if x in mapping]
            if not genes:
                continue
            if len({g.contig for g in genes}) != 1:
                fail(f"DefenseFinder system {row.get('sys_id', '?')} crosses contigs; check gembase IDs")
            system = row.get("subtype") or row.get("type") or "unknown"
            calls.append(SystemCall(genome, "defensefinder", row.get("sys_id", "defensefinder"),
                                    genes[0].contig, min(g.start for g in genes), max(g.end for g in genes),
                                    system, {g.gene_id for g in genes}, {"defensefinder"}))
    return calls


def normal_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def gene_id_aliases(gene_id: str) -> set[str]:
    """Return the common GFF/FASTA spellings of one protein identifier."""
    value = gene_id.strip()
    bare = re.sub(r"^cds-", "", value)
    return {value, bare, f"cds-{bare}"}


def genes_called_by_system(genes: list[Gene], system: SystemCall) -> list[Gene]:
    """Map a caller's component proteins back to genes, with a span fallback."""
    called = set().union(*(gene_id_aliases(x) for x in system.proteins)) if system.proteins else set()
    matched = [g for g in genes if gene_id_aliases(g.gene_id) & called]
    if matched:
        return matched
    # Retain visibility for unusual caller IDs that cannot be mapped directly.
    return [g for g in genes if g.end >= system.start and g.start <= system.end]


def combined_system_call(genome: str, contig: str, source_id: str,
                         systems: list[SystemCall], left: Gene, right: Gene) -> SystemCall:
    """Create the display-level call for all systems inside one boundary pair."""
    if not systems:
        return SystemCall(genome, "boundary_scan", source_id, contig,
                          left.end + 1, right.start - 1, "no defence detected", set(), {"boundary_scan"})
    return SystemCall(
        genome, "boundary_scan", source_id, contig,
        min(s.start for s in systems), max(s.end for s in systems),
        "; ".join(sorted({s.system for s in systems})),
        set().union(*(s.proteins for s in systems)),
        set().union(*(s.routes for s in systems)) | {"boundary_scan"},
    )


def overlap_fraction(a: SystemCall, b: SystemCall) -> float:
    overlap = max(0, min(a.end, b.end) - max(a.start, b.start) + 1)
    return overlap / max(1, min(a.end - a.start + 1, b.end - b.start + 1))


def dereplicate(calls: list[SystemCall], min_overlap: float) -> list[SystemCall]:
    # First collapse genuine duplicates within each caller. Cross-caller matching
    # is handled separately so one broad call cannot absorb multiple neighbours.
    single_caller: list[SystemCall] = []
    for call in sorted(calls, key=lambda x: (x.genome, x.contig, x.start, x.end)):
        match = next((x for x in single_caller if x.genome == call.genome and x.contig == call.contig
                      and x.route == call.route
                      and overlap_fraction(x, call) >= min_overlap
                      and normal_name(x.system) == normal_name(call.system)), None)
        if match:
            match.start, match.end = min(match.start, call.start), max(match.end, call.end)
            match.proteins |= call.proteins
            match.routes |= call.routes
            match.source_id = ",".join(sorted(set(match.source_id.split(",")) | {call.source_id}))
        else:
            single_caller.append(call)

    candidates = []
    for i, left in enumerate(single_caller):
        for j in range(i + 1, len(single_caller)):
            right = single_caller[j]
            if left.route == right.route or left.genome != right.genome or left.contig != right.contig:
                continue
            overlap = overlap_fraction(left, right)
            if overlap < min_overlap:
                continue
            union = left.proteins | right.proteins
            protein_jaccard = len(left.proteins & right.proteins) / len(union) if union else 0
            candidates.append((overlap + protein_jaccard, overlap, i, j))
    used, merged = set(), []
    for _, _, i, j in sorted(candidates, reverse=True):
        if i in used or j in used:
            continue
        left, right = single_caller[i], single_caller[j]
        left.start, left.end = min(left.start, right.start), max(left.end, right.end)
        left.proteins |= right.proteins
        left.routes |= right.routes
        left.source_id = ",".join(sorted(set(left.source_id.split(",")) | set(right.source_id.split(","))))
        if normal_name(left.system) != normal_name(right.system):
            left.system = "/".join(sorted(set(left.system.split("/")) | set(right.system.split("/"))))
        merged.append(left)
        used |= {i, j}
    merged.extend(call for i, call in enumerate(single_caller) if i not in used)
    merged.sort(key=lambda x: (x.genome, x.contig, x.start, x.end))
    return merged


def tsv(path: Path, header: list[str], rows) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def load_precomputed_systems(path: Path, allowed_genomes: set[str]) -> list[SystemCall]:
    """Load already-dereplicated PADLOC/DefenseFinder calls.

    This is used by downstream MGE analyses, where rerunning the callers on an
    excised element would both waste work and lose the original genomic
    coordinate context.
    """
    if not path.is_file():
        fail(f"precomputed systems table not found: {path}")
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"genome", "contig", "start", "end", "system"}
    if rows and not required.issubset(rows[0]):
        fail(f"precomputed systems table lacks columns: {', '.join(sorted(required - set(rows[0])))}")
    systems = []
    for row in rows:
        if row["genome"] not in allowed_genomes:
            continue
        callers = {x for x in row.get("callers", "precomputed").split(",") if x}
        source = row.get("source_ids") or row.get("source_id") or row.get("locus_id") or "precomputed"
        systems.append(SystemCall(
            genome=row["genome"], route="precomputed", source_id=source,
            contig=row["contig"], start=int(row["start"]), end=int(row["end"]),
            system=row["system"], routes=callers or {"precomputed"},
        ))
    return systems


def load_excluded_boundary_pairs(path: Path | None) -> set[tuple[str, str]]:
    if path is None:
        return set()
    if not path.is_file():
        fail(f"excluded boundary-pair table not found: {path}")
    pairs = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"left_boundary_cluster", "right_boundary_cluster"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            fail("excluded boundary-pair table must contain left_boundary_cluster and right_boundary_cluster")
        for row in reader:
            pairs.add((row["left_boundary_cluster"], row["right_boundary_cluster"]))
    return pairs


def write_genes_gff(path: Path, genes: list[Gene], source: str = "PhorificV2",
                    seqid_prefix: str = "") -> None:
    with path.open("w") as handle:
        handle.write("##gff-version 3\n")
        for g in genes:
            product = g.product.replace(";", ",")
            seqid = f"{seqid_prefix}{g.contig}"
            handle.write(f"{seqid}\t{source}\tCDS\t{g.start}\t{g.end}\t.\t{g.strand}\t.\t"
                         f"ID={g.gene_id};product={product};cluster={g.cluster};"
                         f"defence={str(g.defence).lower()}\n")


def svg_locus(path: Path, locus_id: str, genes: list[Gene], left: Gene, right: Gene) -> None:
    width, row_y, margin = 1400, 90, 60
    lo, hi = left.start, right.end
    scale = (width - 2 * margin) / max(1, hi - lo + 1)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="190" viewBox="0 0 {width} 190">',
             '<style>text{font-family:Arial,sans-serif;font-size:11px}.title{font-size:16px;font-weight:bold}</style>',
             f'<text x="{margin}" y="25" class="title">{html.escape(locus_id)}</text>',
             f'<line x1="{margin}" y1="{row_y}" x2="{width-margin}" y2="{row_y}" stroke="#555"/>']
    for gene in genes:
        x1 = margin + (gene.start - lo) * scale
        x2 = margin + (gene.end - lo + 1) * scale
        tip = min(10, max(2, (x2 - x1) / 3))
        if gene.strand == "-":
            points = [(x2, row_y-12), (x1+tip, row_y-12), (x1, row_y), (x1+tip, row_y+12), (x2, row_y+12)]
        else:
            points = [(x1, row_y-12), (x2-tip, row_y-12), (x2, row_y), (x2-tip, row_y+12), (x1, row_y+12)]
        colour = "#d73027" if gene.defence else ("#4575b4" if gene in (left, right) else "#bdbdbd")
        ptext = " ".join(f"{x:.1f},{y}" for x, y in points)
        title = html.escape(f"{gene.gene_id} | {gene.product} | {gene.start}-{gene.end} | {gene.cluster}")
        parts.append(f'<polygon points="{ptext}" fill="{colour}" stroke="#333"><title>{title}</title></polygon>')
        if gene in (left, right) or gene.defence:
            parts.append(f'<text x="{(x1+x2)/2:.1f}" y="{row_y+30}" text-anchor="middle" transform="rotate(30 {(x1+x2)/2:.1f},{row_y+30})">{html.escape(gene.gene_id)}</text>')
    parts += [f'<text x="{margin}" y="175">boundary</text><rect x="115" y="165" width="14" height="10" fill="#4575b4"/>',
              '<text x="150" y="175">defence</text><rect x="205" y="165" width="14" height="10" fill="#d73027"/>',
              '<text x="240" y="175">other annotated / hypothetical genes</text><rect x="435" y="165" width="14" height="10" fill="#bdbdbd"/>', '</svg>']
    path.write_text("\n".join(parts) + "\n")


def pdf_escape(value: str) -> str:
    return value.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def write_vector_pdf(path: Path, width: int, height: int, commands: list[str]) -> None:
    """Write a one-page vector PDF using only built-in Helvetica."""
    stream = ("\n".join(commands) + "\n").encode("latin-1", errors="replace")
    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        (f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {width} {height}] "
         "/Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>").encode(),
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
        f"<< /Length {len(stream)} >>\nstream\n".encode() + stream + b"endstream",
    ]
    output = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offsets = [0]
    for number, obj in enumerate(objects, 1):
        offsets.append(len(output))
        output.extend(f"{number} 0 obj\n".encode() + obj + b"\nendobj\n")
    xref = len(output)
    output.extend(f"xref\n0 {len(objects)+1}\n".encode())
    output.extend(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        output.extend(f"{offset:010d} 00000 n \n".encode())
    output.extend((f"trailer\n<< /Size {len(objects)+1} /Root 1 0 R >>\n"
                   f"startxref\n{xref}\n%%EOF\n").encode())
    path.write_bytes(output)


def pdf_hotspot_group(path: Path, hotspot: str, instances, locus_genes,
                      system_for_locus, systems_for_locus, palette: str = "classic",
                      edge_anchor_side: str | None = None) -> None:
    """Draw every instance of one boundary-pair architecture in one vector PDF."""
    width, left_margin, right_margin = 1700, 440, 60
    spans = [right.end - left.start + 1 for _, left, right, _ in instances]
    max_span = max(spans)
    scale = (width - left_margin - right_margin) / max_span
    row_heights = [78 + 14 * max(0, len(systems_for_locus.get(lid, [])) - 1)
                   for lid, _, _, _ in instances]
    height = 75 + sum(row_heights) + 55
    commands = ["0.8 w"]
    if palette == "peach":
        boundary_colour = (0.95, 0.63, 0.40)
        defence_colour = (0.31, 0.39, 0.57)
        other_colour = (0.76, 0.80, 0.85)
        bracket_colour = (0.22, 0.28, 0.43)
    else:
        boundary_colour = (0.27, 0.46, 0.71)
        defence_colour = (0.84, 0.19, 0.15)
        other_colour = (0.74, 0.74, 0.74)
        bracket_colour = (0.65, 0, 0.15)

    def py(y_from_top: float) -> float:
        return height - y_from_top

    def text(x: float, y: float, value: str, size: int = 10, colour=(0, 0, 0), centred=False) -> None:
        safe = value.encode("latin-1", errors="replace").decode("latin-1")
        if centred:
            x -= len(safe) * size * 0.25
        commands.append(f"{colour[0]} {colour[1]} {colour[2]} rg BT /F1 {size} Tf {x:.1f} {py(y):.1f} Td ({pdf_escape(safe)}) Tj ET")

    def line(x1: float, y1: float, x2: float, y2: float, colour=(0, 0, 0), line_width=0.8) -> None:
        commands.append(f"{colour[0]} {colour[1]} {colour[2]} RG {line_width} w {x1:.1f} {py(y1):.1f} m {x2:.1f} {py(y2):.1f} l S")

    def polygon(points, fill) -> None:
        commands.append(f"{fill[0]} {fill[1]} {fill[2]} rg 0.2 0.2 0.2 RG")
        path_cmd = [f"{points[0][0]:.1f} {py(points[0][1]):.1f} m"]
        path_cmd += [f"{x:.1f} {py(y):.1f} l" for x, y in points[1:]]
        commands.append(" ".join(path_cmd) + " h B")

    if edge_anchor_side == "left":
        title = f"{hotspot} - {len(instances)} instances aligned at conserved right boundary; open left edge"
    elif edge_anchor_side == "right":
        title = f"{hotspot} - {len(instances)} instances aligned at conserved left boundary; open right edge"
    else:
        title = f"{hotspot} - {len(instances)} instances aligned at left boundary"
    text(30, 28, title, 17)
    row_top = 75
    for row_index, (lid, left, right, _) in enumerate(instances):
        system = system_for_locus[lid]
        locus_systems = systems_for_locus.get(lid, [])
        y = row_top + 14 * max(0, len(locus_systems) - 1)
        row_top += row_heights[row_index]
        between = [g for g in locus_genes[lid] if g.start >= left.start and g.end <= right.end]
        if system.genome == system.contig and "__" in system.genome:
            parts = system.genome.split("__", 2)
            display_location = f"{parts[1]} | {parts[2]}:{left.start}-{right.end}" if len(parts) == 3 else system.genome
        else:
            display_location = f"{system.genome} | {system.contig}:{left.start}-{right.end}"
        # Keep the row identifier concise and clear of the gene arrows. Defence
        # systems are already labelled once above their corresponding blocks.
        text(20, y + 3, display_location, 9)
        locus_span = right.end - left.start + 1
        row_x0 = (left_margin + (max_span - locus_span) * scale
                  if edge_anchor_side == "left" else left_margin)
        line(row_x0, y, row_x0 + locus_span*scale, y, (0.33, 0.33, 0.33))
        for gene in between:
            x1 = row_x0 + (gene.start - left.start) * scale
            x2 = row_x0 + (gene.end - left.start + 1) * scale
            tip = min(9, max(2, (x2 - x1) / 3))
            if gene.strand == "-":
                points = [(x2, y-11), (x1+tip, y-11), (x1, y), (x1+tip, y+11), (x2, y+11)]
            else:
                points = [(x1, y-11), (x2-tip, y-11), (x2, y), (x2-tip, y+11), (x1, y+11)]
            # Boundary identity is the primary visual encoding: keep both ends
            # orange even when a boundary gene is itself part of a defence call.
            # The system bracket above the gene still records its defence role.
            is_edge_anchor = ((edge_anchor_side == "left" and gene is left)
                              or (edge_anchor_side == "right" and gene is right))
            is_boundary = gene in (left, right) and not is_edge_anchor
            colour = (boundary_colour if is_boundary
                      else defence_colour if gene.defence else other_colour)
            polygon(points, colour)
            if gene in (left, right):
                text((x1+x2)/2, y+27, gene.gene_id, 8, centred=True)
        if locus_systems:
            for system_index, internal_system in enumerate(locus_systems):
                defence_genes = genes_called_by_system(between, internal_system)
                if not defence_genes:
                    continue
                block_x1 = row_x0 + (min(g.start for g in defence_genes) - left.start) * scale
                block_x2 = row_x0 + (max(g.end for g in defence_genes) - left.start + 1) * scale
                bracket_y = y - 20 - 14 * system_index
                line(block_x1, bracket_y, block_x2, bracket_y, bracket_colour, 1.2)
                line(block_x1, bracket_y, block_x1, bracket_y+5, bracket_colour, 1.2)
                line(block_x2, bracket_y, block_x2, bracket_y+5, bracket_colour, 1.2)
                text((block_x1+block_x2)/2, bracket_y-5, internal_system.system,
                     10, bracket_colour, centred=True)
    legend_y = height - 25
    polygon([(30, legend_y-10), (44, legend_y-10), (44, legend_y), (30, legend_y)], boundary_colour)
    text(50, legend_y, "boundary", 10)
    polygon([(125, legend_y-10), (139, legend_y-10), (139, legend_y), (125, legend_y)], defence_colour)
    text(145, legend_y, "defence", 10)
    polygon([(215, legend_y-10), (229, legend_y-10), (229, legend_y), (215, legend_y)], other_colour)
    text(235, legend_y, "other/hypothetical gene", 10)
    write_vector_pdf(path, width, height, commands)


def main() -> None:
    p = argparse.ArgumentParser(
        description="PhorificV2: PADLOC + DefenseFinder comparative defence-hotspot finder"
    )
    inputs = p.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--input-dir", type=Path,
                        help="folder recursively containing same-stem .gff/.gff3 and .faa files; .fna is optional")
    inputs.add_argument("--manifest", type=Path, help="optional explicit TSV input table")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--window-kb", type=float, default=20.0)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--min-seq-id", type=float, default=0.5)
    p.add_argument("--coverage", type=float, default=0.8)
    p.add_argument("--min-call-overlap", type=float, default=0.5)
    p.add_argument("--exclude-system", action="append", default=["VSPR"],
                   help="caller system name to exclude before hotspot inference; repeatable (default: VSPR)")
    p.add_argument("--min-boundary-frequency", type=float, default=0.5)
    p.add_argument("--min-comparable-loci", type=int, default=3,
                   help="minimum loci in a defence-system group before boundary inference")
    p.add_argument("--min-hotspot-members", type=int, default=2,
                   help="minimum instances required to emit a hotspot locus representation")
    p.add_argument("--boundary-scan-span-factor", type=float, default=2.0,
                   help="maximum boundary-match span as a multiple of the longest detected instance")
    p.add_argument("--min-side-balance", type=float, default=0.5,
                   help="min(left frequency,right frequency)/max(...)")
    p.add_argument("--max-boundary-distance-kb", type=float, default=None)
    p.add_argument("--padloc", default="padloc")
    p.add_argument("--padloc-env", default=None, help="Conda environment name or path containing PADLOC")
    p.add_argument("--defense-finder", default="defense-finder")
    p.add_argument("--defense-finder-env", default=None,
                   help="Conda environment name or path containing DefenseFinder")
    p.add_argument("--mmseqs", default="mmseqs")
    p.add_argument("--mmseqs-env", default=None, help="optional Conda environment name or path containing MMseqs2")
    p.add_argument("--env-runner", default="conda", help="Conda-compatible executable used to run environments")
    p.add_argument("--skip-padloc", action="store_true", help="reuse PADLOC files already in outdir/padloc")
    p.add_argument("--skip-defensefinder", action="store_true",
                   help="reuse DefenseFinder files already in outdir/defensefinder")
    p.add_argument("--skip-mmseqs", action="store_true", help="reuse outdir/mmseqs/all_proteins_cluster.tsv")
    p.add_argument("--representation-only", action="store_true",
                   help="never run external tools; require cached PADLOC, DefenseFinder and MMseqs results")
    p.add_argument("--precomputed-systems", type=Path,
                   help="already-dereplicated system calls; bypass PADLOC and DefenseFinder")
    p.add_argument("--precomputed-clusters", type=Path,
                   help="MMseqs-style representative/member table; bypass protein clustering")
    p.add_argument("--exclude-boundary-pairs", type=Path,
                   help="TSV of existing boundary-cluster pairs that must not become hotspots")
    p.add_argument("--hotspot-prefix", default="hotspot",
                   help="output hotspot ID prefix (default: hotspot)")
    p.add_argument("--plot-palette", choices=("classic", "peach"), default="classic",
                   help="comparative-locus colour palette (default: classic)")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()
    if (args.window_kb <= 0 or not (0 <= args.min_side_balance <= 1)
            or args.min_hotspot_members < 1 or args.boundary_scan_span_factor < 1):
        fail("window and threshold arguments are out of range")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.hotspot_prefix):
        fail("--hotspot-prefix may contain only letters, numbers, dot, underscore and hyphen")
    rows = discover_genomes(args.input_dir) if args.input_dir else load_manifest(args.manifest)
    out = args.outdir.resolve()
    for sub in ("padloc", "defensefinder", "neighbourhoods", "mmseqs", "loci", "logs"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    tsv(out / "resolved_inputs.tsv", ["genome_id", "fna", "gff", "faa"],
        ((r["genome_id"], r["fna"], r["gff"], r["faa"]) for r in rows))
    print(f"Discovered {len(rows)} complete genome(s).")
    calls: list[SystemCall] = []
    genes_by_genome: dict[str, list[Gene]] = {}
    contigs_by_genome: dict[str, dict[str, str]] = {}
    lengths_by_genome: dict[str, dict[str, int]] = {}
    for row in rows:
        genome = row["genome_id"]
        proteins = read_fasta(Path(row["faa"]))
        genes, lengths = read_gff(Path(row["gff"]), genome, proteins)
        genes_by_genome[genome], lengths_by_genome[genome] = genes, lengths
        contigs_by_genome[genome] = read_fasta(Path(row["fna"])) if row.get("fna") else {}
        if not args.precomputed_systems:
            padloc_dest = out / "padloc" / genome
            padloc_dest.mkdir(parents=True, exist_ok=True)
            cached_padloc = padloc_csv(padloc_dest)
            padloc_complete = padloc_run_complete(padloc_dest, out / "logs" / f"padloc.{genome}.log")
            if padloc_complete:
                print(f"Reusing PADLOC results for {genome}.")
            elif args.representation_only:
                fail(f"representation-only mode: cached PADLOC results missing for {genome} in {padloc_dest}")
            elif args.skip_padloc:
                fail(f"--skip-padloc was set but cached PADLOC results are missing for {genome} in {padloc_dest}")
            else:
                if not args.dry_run:
                    require_tool(args.env_runner, args.padloc_env, args.padloc, "PADLOC")
                run(env_command(args.env_runner, args.padloc_env, args.padloc,
                                ["--outdir", str(padloc_dest), "--cpu", str(args.threads),
                                 "--faa", row["faa"], "--gff", row["gff"]]),
                    out / "logs" / f"padloc.{genome}.log", args.dry_run)
                if not args.dry_run:
                    (padloc_dest / ".hotspot_pipeline.complete").write_text("complete\n")

            df_dest = out / "defensefinder" / genome
            df_dest.mkdir(parents=True, exist_ok=True)
            df_faa = df_dest / f"{genome}.gembase.faa"
            df_mapping = make_defensefinder_faa(genes, df_faa)
            tsv(df_dest / "id_map.tsv", ["defensefinder_id", "original_gene_id", "contig", "start", "end"],
                ((df_id, gene.gene_id, gene.contig, gene.start, gene.end) for df_id, gene in df_mapping.items()))
            cached_defensefinder = defensefinder_systems(df_dest)
            defensefinder_complete = (cached_defensefinder is not None
                                      or (df_dest / ".hotspot_pipeline.complete").is_file())
            if defensefinder_complete:
                print(f"Reusing DefenseFinder results for {genome}.")
            elif args.representation_only:
                fail(f"representation-only mode: cached DefenseFinder results missing for {genome} in {df_dest}")
            elif args.skip_defensefinder:
                fail(f"--skip-defensefinder was set but cached results are missing for {genome} in {df_dest}")
            else:
                if not args.dry_run:
                    require_tool(args.env_runner, args.defense_finder_env, args.defense_finder, "DefenseFinder")
                run(env_command(args.env_runner, args.defense_finder_env, args.defense_finder,
                                ["run", "--out-dir", str(df_dest), "--workers", str(args.threads),
                                 "--db-type", "gembase", str(df_faa)]),
                    out / "logs" / f"defensefinder.{genome}.log", args.dry_run)
                if not args.dry_run:
                    (df_dest / ".hotspot_pipeline.complete").write_text("complete\n")

            if not args.dry_run:
                result_csv = padloc_csv(padloc_dest)
                if result_csv:
                    calls += parse_padloc(result_csv, genome, "padloc")
                result_tsv = defensefinder_systems(df_dest)
                if result_tsv:
                    calls += parse_defensefinder(result_tsv, genome, df_mapping)
    if args.dry_run:
        print(f"Dry run complete; commands are in {out / 'logs'}")
        return

    if args.precomputed_systems:
        calls = load_precomputed_systems(args.precomputed_systems.resolve(), set(genes_by_genome))
    raw_call_count = len(calls)
    raw_padloc_call_count = sum("padloc" in c.routes or c.route == "padloc" for c in calls)
    raw_defensefinder_call_count = sum("defensefinder" in c.routes or c.route == "defensefinder" for c in calls)
    excluded_systems = {normal_name(x) for x in args.exclude_system}
    excluded_calls = [c for c in calls if normal_name(c.system) in excluded_systems]
    calls = [c for c in calls if normal_name(c.system) not in excluded_systems]
    if excluded_calls:
        print(f"Excluded {len(excluded_calls)} call(s) from hotspot inference: "
              + ", ".join(sorted({c.system for c in excluded_calls})))
    systems = calls if args.precomputed_systems else dereplicate(calls, args.min_call_overlap)
    tsv(out / "systems.dereplicated.tsv",
        ["locus_id", "genome", "contig", "start", "end", "system", "callers", "source_ids"],
        ((f"locus_{i:06d}", s.genome, s.contig, s.start, s.end, s.system, ",".join(sorted(s.routes)), s.source_id)
         for i, s in enumerate(systems, 1)))
    if not systems:
        fail("PADLOC and DefenseFinder produced no defence-system calls")

    systems_by_location = defaultdict(list)
    for system in systems:
        systems_by_location[(system.genome, system.contig)].append(system)

    def overlaps_any_defence(gene: Gene) -> bool:
        """True when a gene overlaps any retained defence call at its location."""
        return any(system.end >= gene.start and system.start <= gene.end
                   for system in systems_by_location[(gene.genome, gene.contig)])

    window = round(args.window_kb * 1000)
    locus_genes: dict[str, list[Gene]] = {}
    system_for_locus: dict[str, SystemCall] = {}
    systems_for_locus: dict[str, list[SystemCall]] = {}
    for i, system in enumerate(systems, 1):
        lid = f"locus_{i:06d}"
        selected = [replace(g, defence=overlaps_any_defence(g))
                    for g in genes_by_genome[system.genome] if g.contig == system.contig
                    and g.end >= system.start - window and g.start <= system.end + window]
        locus_genes[lid], system_for_locus[lid] = selected, system
        systems_for_locus[lid] = [system]
        write_fasta(((g.uid, g.sequence) for g in selected if g.sequence), out / "neighbourhoods" / f"{lid}.faa")
        tsv(out / "neighbourhoods" / f"{lid}.genes.tsv",
            ["gene_id", "genome", "contig", "start", "end", "strand", "product", "is_defence"],
            ((g.gene_id, g.genome, g.contig, g.start, g.end, g.strand, g.product, int(g.defence)) for g in selected))
        write_genes_gff(out / "neighbourhoods" / f"{lid}.gff3", selected)
    # Cluster every annotated protein, not only defence neighbourhood proteins.
    # This is required to find boundary pairs in genomes with no detected defence.
    unique_records = {g.uid: g.sequence for genes in genes_by_genome.values() for g in genes if g.sequence}
    all_genome_faa = out / "neighbourhoods" / "all_genome_proteins.faa"
    write_fasta(unique_records.items(), all_genome_faa)
    cluster_tsv = out / "mmseqs" / "all_proteins_cluster.tsv"
    if args.precomputed_clusters:
        source_clusters = args.precomputed_clusters.resolve()
        if not source_clusters.is_file():
            fail(f"precomputed cluster table not found: {source_clusters}")
        if source_clusters != cluster_tsv.resolve():
            shutil.copyfile(source_clusters, cluster_tsv)
        print(f"Reusing precomputed protein clusters from {source_clusters}.")
    elif cluster_tsv.is_file():
        print(f"Reusing MMseqs2 clusters from {cluster_tsv}.")
    elif args.representation_only:
        fail(f"representation-only mode: cached MMseqs2 clusters missing: {cluster_tsv}")
    elif args.skip_mmseqs:
        fail(f"--skip-mmseqs was set but cached cluster table is missing: {cluster_tsv}")
    else:
        require_tool(args.env_runner, args.mmseqs_env, args.mmseqs, "MMseqs2")
        run(env_command(args.env_runner, args.mmseqs_env, args.mmseqs,
                        ["easy-cluster", str(all_genome_faa),
                         str(out / "mmseqs" / "all_proteins"), str(out / "mmseqs" / "tmp_all_proteins"),
                         "--min-seq-id", str(args.min_seq_id), "-c", str(args.coverage), "--cov-mode", "0",
                         "--threads", str(args.threads)]), out / "logs" / "mmseqs.log")
    if not cluster_tsv.is_file():
        fail(f"MMseqs cluster table not found: {cluster_tsv}")
    member_to_cluster = {}
    with cluster_tsv.open() as handle:
        for line in handle:
            rep, member = line.rstrip("\n").split("\t")[:2]
            member_to_cluster[member] = rep
    for genes in genes_by_genome.values():
        for gene in genes:
            gene.cluster = member_to_cluster.get(gene.uid, "")
    for genes in locus_genes.values():
        for gene in genes:
            gene.cluster = member_to_cluster.get(gene.uid, "")

    # Count each cluster once per locus and side; overlapping windows do not inflate within-locus counts.
    left_count, right_count, both_count, group_size = Counter(), Counter(), Counter(), Counter()
    per_locus_sides = {}
    for lid, genes in locus_genes.items():
        s = system_for_locus[lid]
        group = normal_name(s.system)
        group_size[group] += 1
        left = {g.cluster for g in genes if g.cluster and not g.defence and g.end < s.start}
        right = {g.cluster for g in genes if g.cluster and not g.defence and g.start > s.end}
        per_locus_sides[lid] = (group, left, right)
        left_count.update((group, x) for x in left)
        right_count.update((group, x) for x in right)
        both_count.update((group, x) for x in left | right)
    stats = []
    qualifying_left, qualifying_right = set(), set()
    for group, cluster in sorted(set(left_count) | set(right_count)):
        n = group_size[group]
        lf, rf, co = left_count[(group, cluster)] / n, right_count[(group, cluster)] / n, both_count[(group, cluster)] / n
        balance = min(lf, rf) / max(lf, rf) if max(lf, rf) else 0
        enough_data = n >= args.min_comparable_loci
        left_ok = enough_data and lf >= args.min_boundary_frequency
        right_ok = enough_data and rf >= args.min_boundary_frequency
        if left_ok:
            qualifying_left.add((group, cluster))
        if right_ok:
            qualifying_right.add((group, cluster))
        stats.append((group, n, int(enough_data), cluster, left_count[(group, cluster)], right_count[(group, cluster)], both_count[(group, cluster)], f"{lf:.6f}",
                      f"{rf:.6f}", f"{co:.6f}", f"{balance:.6f}", int(left_ok), int(right_ok)))
    tsv(out / "cluster_frequencies.tsv",
        ["system_group", "group_loci", "enough_comparable_loci", "cluster_id", "left_loci", "right_loci", "any_side_loci", "left_frequency", "right_frequency",
         "cooccurrence_frequency", "same_cluster_side_balance", "qualifies_left", "qualifies_right"], stats)

    max_dist = round((args.max_boundary_distance_kb or args.window_kb) * 1000)
    hotspot_rows, architecture_count = [], Counter()
    signature_pairs = {}
    selected_loci = []
    excluded_boundary_pairs = load_excluded_boundary_pairs(
        args.exclude_boundary_pairs.resolve() if args.exclude_boundary_pairs else None)
    for lid, genes in locus_genes.items():
        s = system_for_locus[lid]
        group = normal_name(s.system)
        lefts = [g for g in genes if (group, g.cluster) in qualifying_left and not g.defence and g.end < s.start and s.start-g.end <= max_dist]
        rights = [g for g in genes if (group, g.cluster) in qualifying_right and not g.defence and g.start > s.end and g.start-s.end <= max_dist]
        if not lefts or not rights:
            continue
        left, right = max(lefts, key=lambda g: g.end), min(rights, key=lambda g: g.start)
        if (left.cluster, right.cluster) in excluded_boundary_pairs:
            continue
        n = group_size[group]
        lf, rf = left_count[(group, left.cluster)] / n, right_count[(group, right.cluster)] / n
        pair_balance = min(lf, rf) / max(lf, rf)
        if pair_balance < args.min_side_balance:
            continue
        signature = f"{left.cluster}__{right.cluster}"
        signature_pairs[signature] = (left.cluster, right.cluster)
        architecture_count[signature] += 1
        selected_loci.append((lid, left, right, signature))

    # Multiple defence calls in the same interval can nominate the identical
    # boundary pair. Count that genomic interval once, not once per caller/system.
    unique_selected = {}
    for instance in selected_loci:
        lid, left, right, signature = instance
        anchor = system_for_locus[lid]
        key = (anchor.genome, anchor.contig, left.gene_id, right.gene_id, signature)
        unique_selected.setdefault(key, instance)
    selected_loci = list(unique_selected.values())

    # Search all genomes for the accepted boundary-cluster pairs, including
    # loci where neither caller detected any defence system between them.
    existing_pairs = {(system_for_locus[lid].genome, system_for_locus[lid].contig,
                       left.gene_id, right.gene_id, signature)
                      for lid, left, right, signature in selected_loci}
    max_span_by_signature = defaultdict(int)
    for _, left, right, signature in selected_loci:
        max_span_by_signature[signature] = max(max_span_by_signature[signature], right.end-left.start+1)
    boundary_scan_number = 0
    for signature, (left_cluster, right_cluster) in signature_pairs.items():
        scan_span = round(max_span_by_signature[signature] * args.boundary_scan_span_factor)
        for genome, genome_genes in genes_by_genome.items():
            by_contig = defaultdict(list)
            for gene in genome_genes:
                by_contig[gene.contig].append(gene)
            for contig, contig_genes in by_contig.items():
                ordered = sorted(contig_genes, key=lambda g: (g.start, g.end))
                rights = [g for g in ordered if g.cluster == right_cluster]
                used_right = set()
                for left in (g for g in ordered if g.cluster == left_cluster
                             and not overlaps_any_defence(g)):
                    possible = [g for g in rights if g.gene_id not in used_right and g.start > left.end
                                and g.end-left.start+1 <= scan_span and not overlaps_any_defence(g)]
                    if not possible:
                        continue
                    right = min(possible, key=lambda g: g.start)
                    key = (genome, contig, left.gene_id, right.gene_id, signature)
                    if key in existing_pairs:
                        used_right.add(right.gene_id)
                        continue
                    internal_systems = [s for s in systems_by_location[(genome, contig)]
                                        if s.end >= left.start and s.start <= right.end]
                    boundary_scan_number += 1
                    prefix = "boundary_defence" if internal_systems else "boundary_only"
                    lid = f"{prefix}_{boundary_scan_number:06d}"
                    internal_gene_ids = {
                        g.gene_id for internal in internal_systems
                        for g in genes_called_by_system(ordered, internal)
                    }
                    between = [replace(g, defence=g.gene_id in internal_gene_ids)
                               for g in ordered
                               if g.start >= left.start and g.end <= right.end]
                    # Use the copied boundary objects so identity-based colouring works.
                    copied_left = next(g for g in between if g.gene_id == left.gene_id)
                    copied_right = next(g for g in between if g.gene_id == right.gene_id)
                    locus_genes[lid] = between
                    systems_for_locus[lid] = internal_systems
                    system_for_locus[lid] = combined_system_call(
                        genome, contig, lid, internal_systems, copied_left, copied_right)
                    selected_loci.append((lid, copied_left, copied_right, signature))
                    existing_pairs.add(key)
                    used_right.add(right.gene_id)

    total_architecture_count = Counter(instance[3] for instance in selected_loci)
    eligible_signatures = [sig for sig in sorted(total_architecture_count)
                           if total_architecture_count[sig] >= args.min_hotspot_members]
    signature_ids = {sig: f"{args.hotspot_prefix}_{i:04d}" for i, sig in enumerate(eligible_signatures, 1)}
    selected_loci = [instance for instance in selected_loci if instance[3] in signature_ids]
    # Re-annotate every retained locus using all calls inside its final boundaries.
    # This captures alternate and multiple systems, including calls recovered by
    # the boundary scan, rather than displaying only the system that seeded it.
    updated_loci = []
    for lid, left, right, signature in selected_loci:
        anchor = system_for_locus[lid]
        internal_systems = [s for s in systems_by_location[(anchor.genome, anchor.contig)]
                            if s.end >= left.start and s.start <= right.end]
        systems_for_locus[lid] = internal_systems
        genes = locus_genes[lid]
        called_gene_ids = {
            g.gene_id for internal in internal_systems for g in genes_called_by_system(genes, internal)
        }
        genes = [replace(g, defence=g.gene_id in called_gene_ids) for g in genes]
        locus_genes[lid] = genes
        copied_left = next(g for g in genes if g.gene_id == left.gene_id)
        copied_right = next(g for g in genes if g.gene_id == right.gene_id)
        # Boundaries flank defence content; they cannot themselves be part of a
        # retained PADLOC/DefenseFinder call in that instance.
        if copied_left.defence or copied_right.defence:
            continue
        system_for_locus[lid] = combined_system_call(
            anchor.genome, anchor.contig, lid, internal_systems, copied_left, copied_right)
        updated_loci.append((lid, copied_left, copied_right, signature))
    selected_loci = updated_loci
    # Recheck prevalence after strict boundary validation and assign stable IDs
    # only to architectures that still meet the requested member threshold.
    total_architecture_count = Counter(instance[3] for instance in selected_loci)
    eligible_signatures = [sig for sig in sorted(total_architecture_count)
                           if total_architecture_count[sig] >= args.min_hotspot_members]
    signature_ids = {sig: f"{args.hotspot_prefix}_{i:04d}" for i, sig in enumerate(eligible_signatures, 1)}
    selected_loci = [instance for instance in selected_loci if instance[3] in signature_ids]
    # These directories are fully generated artifacts. Rebuild them exactly so
    # representations from no-longer-eligible singleton architectures cannot linger.
    loci_root = out / "loci"
    for old_dir in loci_root.glob(f"{args.hotspot_prefix}_[0-9][0-9][0-9][0-9]"):
        if old_dir.is_dir():
            shutil.rmtree(old_dir)
    for lid, left, right, signature in selected_loci:
        s = system_for_locus[lid]
        hotspot = signature_ids[signature]
        defence_present = s.system != "no defence detected"
        hotspot_rows.append((hotspot, lid, s.genome, s.contig, left.start, right.end, left.gene_id, right.gene_id,
                             left.cluster, right.cluster, s.system, int(defence_present),
                             total_architecture_count[signature]))
        locus_dir = out / "loci" / hotspot
        locus_dir.mkdir(parents=True, exist_ok=True)
        between = [g for g in locus_genes[lid] if g.start >= left.start and g.end <= right.end]
        write_fasta(((g.uid, g.sequence) for g in between if g.sequence), locus_dir / f"{lid}.faa")
        seq = contigs_by_genome[s.genome].get(s.contig, "")
        if seq:
            write_fasta([(f"{s.contig}:{left.start}-{right.end}", seq[left.start-1:right.end])], locus_dir / f"{lid}.fna")
        tsv(locus_dir / f"{lid}.genes.tsv",
            ["gene_id", "contig", "start", "end", "strand", "product", "cluster_id", "is_defence"],
            ((g.gene_id, g.contig, g.start, g.end, g.strand, g.product, g.cluster, int(g.defence)) for g in between))
        write_genes_gff(locus_dir / f"{lid}.gff3", between)
    # One comparative figure and one combined data bundle per recurrent hotspot.
    instances_by_hotspot = defaultdict(list)
    for instance in selected_loci:
        instances_by_hotspot[signature_ids[instance[3]]].append(instance)
    for hotspot, instances in sorted(instances_by_hotspot.items()):
        group_dir = out / "loci" / hotspot
        pdf_hotspot_group(group_dir / f"{hotspot}.comparative.pdf", hotspot, instances,
                          locus_genes, system_for_locus, systems_for_locus, args.plot_palette)
        combined_genes = []
        combined_faa = []
        combined_rows = []
        for lid, left, right, _ in instances:
            system = system_for_locus[lid]
            between = [g for g in locus_genes[lid] if g.start >= left.start and g.end <= right.end]
            combined_genes.extend((lid, g) for g in between)
            combined_faa.extend((f"{lid}|{g.uid}", g.sequence) for g in between if g.sequence)
            combined_rows.extend((lid, system.genome, g.gene_id, g.contig, g.start, g.end, g.strand,
                                  g.product, g.cluster, int(g.defence), int(g in (left, right))) for g in between)
        write_fasta(combined_faa, group_dir / f"{hotspot}.all_instances.faa")
        tsv(group_dir / f"{hotspot}.all_instances.genes.tsv",
            ["locus_id", "genome", "gene_id", "contig", "start", "end", "strand", "product",
             "cluster_id", "is_defence", "is_boundary"], combined_rows)
        with (group_dir / f"{hotspot}.all_instances.gff3").open("w") as handle:
            handle.write("##gff-version 3\n")
            for lid, g in combined_genes:
                product = g.product.replace(";", ",")
                handle.write(f"{lid}|{g.genome}|{g.contig}\tPhorificV2\tCDS\t{g.start}\t{g.end}\t.\t{g.strand}\t.\t"
                             f"ID={lid}|{g.gene_id};product={product};cluster={g.cluster};"
                             f"defence={str(g.defence).lower()}\n")
    tsv(out / "hotspots.tsv",
        ["hotspot_id", "locus_id", "genome", "contig", "start", "end", "left_boundary_id", "right_boundary_id",
         "left_boundary_cluster", "right_boundary_cluster", "defence_system", "defence_present",
         "instances_of_architecture"], hotspot_rows)
    summary = {"genomes": len(rows), "raw_caller_calls": raw_call_count,
               "padloc_calls": raw_padloc_call_count,
               "defensefinder_calls": raw_defensefinder_call_count,
               "dereplicated_systems": len(systems),
               "protein_clusters": len(set(member_to_cluster.values())),
               "qualifying_left_clusters": len(qualifying_left), "qualifying_right_clusters": len(qualifying_right),
               "hotspot_architectures": len(signature_ids), "hotspot_instances": len(hotspot_rows),
               "boundary_only_instances": sum(row[-2] == 0 for row in hotspot_rows),
               "boundary_recovered_defence_instances": sum(row[1].startswith("boundary_defence_") for row in hotspot_rows),
               "excluded_system_calls": len(excluded_calls),
               "excluded_existing_boundary_pairs": len(excluded_boundary_pairs),
               "parameters": vars(args) | {"manifest": str(args.manifest) if args.manifest else None,
                                              "input_dir": str(args.input_dir) if args.input_dir else None,
                                              "outdir": str(out)}}
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "parameters"}, indent=2))


if __name__ == "__main__":
    main()
