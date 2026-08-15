# PHORIFIC: A Defence hotspot pipeline

It identifies recurrent defence-system hotspots across a
collection of annotated prokaryotic genomes.

It combines PADLOC and DefenseFinder calls, dereplicates overlapping predictions,
clusters proteins with MMseqs2, learns recurrent flanking boundary proteins, and
searches every genome for both defence-containing and boundary-only versions of
each hotspot. Final architectures are presented as grouped vector PDFs.

## Pipeline overview

For every genome, the pipeline:

1. discovers a matching GFF3 and protein FASTA;
2. runs PADLOC on the `FAA + GFF3` pair;
3. runs DefenseFinder on an ordered, multi-replicon-safe copy of the FAA;
4. maps DefenseFinder identifiers back to the original GFF3 genes;
5. merges overlapping PADLOC and DefenseFinder system calls;
6. extracts a  neighbourhood around each dereplicated system;
7. clusters every annotated protein from every input genome with MMseqs2;
8. measures how often non-defence protein clusters occur on each side of defence systems;
9. selects recurrent left and right boundary clusters;
10. groups loci sharing the same boundary-cluster pair into hotspot architectures;
11. searches all contigs for the same boundaries with no detected defence system between them; and
12. writes coordinates, sequences, annotations, and comparative PDF diagrams.

## Requirements

- Python 3.9 or newer. No third-party Python package is required.
- PADLOC 2.x and a compatible PADLOC database.
- DefenseFinder and its compatible models.
- MMseqs2.

PADLOC and DefenseFinder may be installed in different Conda environments. The
pipeline does not update their databases or models automatically.


## Input files

The recommended input is a folder containing matching `.gff`/`.gff3` and `.faa`
files. The pipeline searches recursively, so files may be together or in
per-genome subdirectories.

```text
genomes/
├── strain_A/
│   ├── strain_A.gff3
│   ├── strain_A.faa
│   └── strain_A.fna       # optional
└── strain_B/
    ├── strain_B.gff3
    └── strain_B.faa
```

Requirements:

- GFF and FAA files belonging to one genome must have the same basename.
- Genome basenames must be unique across the input collection.
- FAA identifiers must match a GFF CDS `ID`, `protein_id`, or `locus_tag`.
- An FNA with the same basename is optional. When available, nucleotide locus
  FASTA files are also written.
- If an FNA is supplied, its contig identifiers must match the GFF.

Incomplete sets lacking an FAA are reported and ignored. The complete resolved
input list is saved automatically as `resolved_inputs.tsv`.

### Optional manifest

For an unusual directory layout, use a tab-separated manifest instead of
`--input-dir`:

```text
genome_id	gff	faa	fna
strain_A	/path/strain_A.gff3	/path/strain_A.faa	/path/strain_A.fna
strain_B	/path/strain_B.gff3	/path/strain_B.faa	
```

The required columns are `genome_id`, `gff`, and `faa`; `fna` is optional.

## Packages environments

1. MMseqs2 is taken from the current `PATH`. If it is also in a separate Conda
environment, add `--mmseqs-env ENVIRONMENT_NAME`.

Environment names use `conda run -n NAME`; environment paths use
`conda run -p PATH`. A different compatible runner can be supplied with
`--env-runner`.

2. For PADLOC, provide environment name under ` --padloc-env `
3. For DefenseFinder, provide environment name under ` --defense-finder-env`

## Resume and representation-only mode

The output directory is resumable. Subsequent runs check for and reuse completed
PADLOC, DefenseFinder, and all-protein MMseqs results. PADLOC runs that completed
with no systems are recognized from their domain table and `Nothing found` log.
Missing stages are run normally.

To guarantee that no external analysis program is launched, use:

```bash
python3 hotspot_pipeline.py \
  --input-dir collected_gff3 \
  --outdir hotspot_results \
  --window-kb 20 \
  --threads 8 \
  --representation-only \
  --min-hotspot-members 2 \
  --boundary-scan-span-factor 2
```

Representation-only mode reparses cached results and rebuilds hotspot tables,
locus bundles, and PDFs. It stops if a required cached result is missing.

Changing boundary-frequency, span, membership, or drawing parameters is safe in
representation-only mode. Changing `--min-seq-id` or `--coverage` requires new
MMseqs clustering; use a new output directory or remove the complete old
`outdir/mmseqs/` cache before rerunning.

## PADLOC and DefenseFinder dereplication

The callers often use different system names, protein names, and system extents.
Predictions are therefore not matched by annotation text.

Within one caller, overlapping calls are merged only when their normalized
system names agree. PADLOC-to-DefenseFinder candidates are matched using the same
genome and contig, interval overlap over the shorter call, and shared original
protein identifiers when available.

Cross-caller matching is one-to-one, preventing one broad prediction from
absorbing multiple neighbouring systems. The union coordinate interval and
protein identifiers are retained. `systems.dereplicated.tsv` records caller
support and original call IDs.

## Boundary inference

For each defence system, genes within `--window-kb` on both sides are extracted.
All proteins across all genomes are clustered with MMseqs2.

Within each normalized defence-system group, the pipeline counts whether each
non-defence protein cluster occurs to the left, to the right, or on either side
of defence.

A side-specific boundary candidate must meet `--min-boundary-frequency`. At
least `--min-comparable-loci` systems must exist before frequency-based inference
is permitted. The closest qualifying cluster on each side becomes the
provisional boundary, and their frequencies must satisfy `--min-side-balance`.

Each distinct ordered pair `(left boundary cluster, right boundary cluster)`
defines a hotspot architecture. By default, only architectures with at least two
total instances are represented.

This is an operational comparative-genomics boundary definition. Supporting
counts are retained in `cluster_frequencies.tsv` for inspection and tuning.

## Genomes with boundaries but no detected defence

After learning boundaries from defence-containing loci, every annotated contig
is searched for the same ordered boundary-cluster pairs.

A boundary-only locus is included when:

- both boundary clusters occur on the same contig;
- the left boundary precedes the right boundary;
- no dereplicated PADLOC or DefenseFinder system overlaps the interval; and
- its span is within `--boundary-scan-span-factor` times the longest detected
  instance of that architecture.

These loci are added to the same hotspot PDF and data bundle with
`defence_system = no defence detected` and `defence_present = 0`.

Increasing the span factor recovers more distant pairs but also increases the
risk of joining unrelated conserved proteins. The default `2.0` is moderately
permissive; `1.0` is strict.

## Principal parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `--window-kb` | `20` | Kilobases extracted on each side of a defence system. |
| `--threads` | `4` | Worker threads supplied to external tools. |
| `--min-call-overlap` | `0.5` | Minimum overlap over the shorter caller interval. |
| `--min-seq-id` | `0.5` | MMseqs minimum protein sequence identity. |
| `--coverage` | `0.8` | MMseqs bidirectional alignment coverage threshold. |
| `--min-boundary-frequency` | `0.5` | Minimum fraction of comparable loci containing a side-specific boundary cluster. |
| `--min-comparable-loci` | `3` | Minimum systems in a group before boundary inference. |
| `--min-side-balance` | `0.5` | Minimum ratio between selected left- and right-boundary frequencies. |
| `--max-boundary-distance-kb` | window size | Maximum defence-to-boundary distance. |
| `--min-hotspot-members` | `2` | Minimum total instances required for locus output and a PDF. |
| `--boundary-scan-span-factor` | `2.0` | Maximum boundary-only span relative to the longest detected instance. |

Run `python3 PHORIFICV2.py --help` for executable-path, environment,
resume, skip, and dry-run options.

## Output structure

```text
hotspot_results/
├── resolved_inputs.tsv
├── systems.dereplicated.tsv
├── cluster_frequencies.tsv
├── hotspots.tsv
├── summary.json
├── padloc/<genome>/
├── defensefinder/<genome>/
│   └── id_map.tsv
├── mmseqs/
│   └── all_proteins_cluster.tsv
├── neighbourhoods/
│   ├── all_genome_proteins.faa
│   ├── locus_XXXXXX.faa
│   ├── locus_XXXXXX.gff3
│   └── locus_XXXXXX.genes.tsv
└── loci/
    └── hotspot_XXXX/
        ├── hotspot_XXXX.comparative.pdf
        ├── hotspot_XXXX.all_instances.faa
        ├── hotspot_XXXX.all_instances.gff3
        ├── hotspot_XXXX.all_instances.genes.tsv
        ├── {locus_XXXXXX,boundary_only_XXXXXX}.faa
        ├── {locus_XXXXXX,boundary_only_XXXXXX}.gff3
        ├── {locus_XXXXXX,boundary_only_XXXXXX}.genes.tsv
        └── {locus_XXXXXX,boundary_only_XXXXXX}.fna  # only with input FNA
```

### Key tables

`systems.dereplicated.tsv`

- coordinates and annotation of every merged caller system;
- `callers`: PADLOC, DefenseFinder, or both;
- `source_ids`: original caller system identifiers.

`cluster_frequencies.tsv`

- side-specific protein-cluster counts and frequencies;
- comparable group size;
- whether each cluster qualifies on the left or right.

`hotspots.tsv`

- hotspot and locus IDs;
- genome, contig, start, and end;
- original boundary gene IDs and MMseqs clusters;
- defence-system annotation;
- `defence_present`: `1` for caller-detected defence and `0` for boundary-only loci;
- total instances of the architecture.

`summary.json`

- genome and caller counts;
- dereplicated-system and protein-cluster counts;
- hotspot architecture and instance counts;
- number of boundary-only instances;
- exact command-line parameters used.



## Interpretation notes

- PADLOC and DefenseFinder disagreement is preserved in the merged system name
  and caller-support columns.
- A boundary-only locus means neither caller identified a system within the
  selected interval; it does not prove the locus has no defensive function.
- Very large span factors can join unrelated occurrences of common clusters.
- A genome may be absent even when it contains defence systems if its flanking
  proteins do not meet the recurrent-boundary criteria.
- A genome with no detected defence may appear when it contains a complete
  accepted boundary pair within the configured span.
