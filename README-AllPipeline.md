# PhorificV2

PhorificV2 identifies recurrent defence-system hotspots across annotated
prokaryotic genomes. It combines PADLOC and DefenseFinder predictions, learns
conserved protein boundaries around defence loci, and recovers every locus
between those boundaries, including loci with no detected defence system.

## What it does

1. Finds matching GFF3 and FAA files for each genome.
2. Runs PADLOC and DefenseFinder.
3. Merges overlapping calls from the two callers while retaining both names
   when their annotations differ.
4. Extracts a configurable neighbourhood around each defence locus.
5. Clusters proteins across the dataset with MMseqs2.
6. Identifies recurrent non-defence protein clusters on both sides of defence.
7. Uses recurrent left/right cluster pairs as hotspot boundaries.
8. Searches every genome for those boundary pairs.
9. Retains defence-positive, alternative-defence, and boundary-only loci.
10. Groups loci with the same boundaries and produces coordinates, FAA/GFF3
    bundles, summary tables, and comparative PDF diagrams.

Multiple independent hotspot families can be recovered in one analysis.

## Main scripts

| Script | Purpose |
|---|---|
| `phorificV2.py` | Main hotspot-discovery command. |
| `hotspot_pipeline.py` | Implementation imported by `phorificV2.py`. |
| `summarise_defences.py` | Optional defence and hotspot summaries. |
| `phorificV2_annotation.py` | Optional Pfam, AMRFinderPlus, and VFDB annotation. |
| `hotspot_annotation_pipeline.py` | Annotation implementation. |
| `phorific_mge.py` | Optional plasmid/prophage analysis after geNomad. |
| `phorific_ice.py` | Optional ICE/IME/rICE analysis after ICEfinder-opt. |

`pfam_priority.py` and the `plot_hotspot_*.py` files are internal annotation
and plotting helpers. They are not separate stages to run.

## Requirements

Main hotspot discovery requires:

- Python 3.9+
- PADLOC 2.x and its database
- DefenseFinder and its models
- MMseqs2

Optional annotation requires HMMER/Pfam, AMRFinderPlus, DIAMOND, a VFDB Set A
protein FASTA, and matplotlib/NumPy/seaborn for plotting.

PADLOC, DefenseFinder, and MMseqs2 may be installed in different Conda
environments.

## Input

The input folder is searched recursively. Each genome requires matching GFF3
and FAA basenames:

```text
genomes/
├── strain_A.gff3
├── strain_A.faa
├── strain_A.fna       # optional
├── strain_B.gff3
└── strain_B.faa
```

FNA files are optional and are used only for nucleotide export. FAA identifiers
must match a GFF3 CDS `ID`, `protein_id`, or `locus_tag`.

## Run PhorificV2

```bash
python3 phorificV2.py \
  --input-dir /path/to/genomes \
  --outdir /path/to/phorific_results \
  --window-kb 20 \
  --threads 8 \
  --padloc-env padloc \
  --defense-finder-env defensefinder \
  --min-hotspot-members 2 \
  --boundary-scan-span-factor 2
```

MMseqs2 is taken from the current `PATH`. Add `--mmseqs-env ENV_NAME` if it is
installed in another Conda environment. See `python3 phorificV2.py --help` for
all thresholds and executable options.

Completed PADLOC, DefenseFinder, and MMseqs2 results are reused automatically.
To rebuild tables and PDFs without launching external programs, rerun with:

```text
--representation-only
```

## Summarise defence systems

```bash
python3 summarise_defences.py \
  --results-dir /path/to/phorific_results \
  --top-types 40
```

## Optional hotspot annotation

Run these stages in order because their external programs may be installed in
different environments:

```bash
python3 phorificV2_annotation.py prepare \
  --results-dir /path/to/phorific_results \
  --outdir /path/to/annotation_results

python3 phorificV2_annotation.py pfam \
  --outdir /path/to/annotation_results \
  --pfam-hmm /path/to/Pfam-A.hmm --threads 8

python3 phorificV2_annotation.py amr \
  --outdir /path/to/annotation_results --threads 8

python3 phorificV2_annotation.py virulence \
  --outdir /path/to/annotation_results \
  --vfdb-proteins /path/to/VFDB_setA_pro.fas.gz \
  --min-identity 35 --min-coverage 70 --threads 8

python3 phorificV2_annotation.py combine \
  --results-dir /path/to/phorific_results \
  --outdir /path/to/annotation_results

python3 phorificV2_annotation.py plots \
  --results-dir /path/to/phorific_results \
  --outdir /path/to/annotation_results
```

The Pfam summary removes proteins already assigned as defence or AMR and
deprioritises generic integration, transposition, regulation, and phage
structural functions. Raw hits remain available for auditing.

## Optional plasmid and prophage analysis

After running geNomad and combining its plasmid and virus summary tables:

```bash
python3 phorific_mge.py \
  --phorific-results /path/to/phorific_results \
  --plasmids-tsv /path/to/all_plasmids.tsv \
  --phages-tsv /path/to/all_phages.tsv \
  --genome-files-dir /path/to/genome_files \
  --outdir /path/to/phorific_mge_results \
  --min-hotspot-prevalence 0.05 \
  --boundary-scan-span-factor 5 \
  --threads 8
```

This classifies existing hotspots as chromosome, plasmid, or prophage and then
searches plasmids and prophages independently for new two-boundary hotspots. It
reuses the original defence calls and MMseqs clusters. Detailed filtering and
rerun options are in `PHORIFIC_MGE_WORKFLOW.md`.

## Optional ICEfinder analysis

After combining ICEfinder-opt `ICE_details.tsv` files:

```bash
python3 phorific_ice.py \
  --phorific-results /path/to/phorific_results \
  --icefinder-tsv /path/to/icefinder_master.tsv \
  --genome-files-dir /path/to/genome_files \
  --outdir /path/to/phorific_ice_results \
  --min-hotspot-prevalence 0.05 \
  --boundary-scan-span-factor 5 \
  --threads 8
```

Keep `plot_icefinder_summary.R` beside `phorific_ice.py` for the ICE summary
figures.

## Principal outputs

```text
phorific_results/
├── resolved_inputs.tsv
├── systems.dereplicated.tsv
├── boundary_cluster_stats.tsv
├── hotspots.tsv
├── mmseqs/
├── neighbourhoods/
└── loci/hotspot_####/
    ├── hotspot_####.comparative.pdf
    ├── hotspot_####.all_instances.faa
    ├── hotspot_####.all_instances.gff3
    └── hotspot_####.all_instances.genes.tsv
```

`hotspots.tsv` contains the genome, contig, coordinates, boundary clusters,
defence content, and hotspot membership of every recovered locus.


