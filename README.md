# phylogenetics_project

Downloads fungal genomes from [JGI MycoCosm](https://mycocosm.jgi.doe.gov), checks assembly and
annotation quality, and builds one clean proteome per genome (one protein per gene) for phylogenomic
analyses. The cohort is listed in `00_inputs/mycocosm_datasets.tsv` (portal IDs in column 1).

## Requirements

- A Slurm cluster with BUSCO 6.1.0 (`fungi_odb12`), SeqKit 2.13.0, and Python 3.12 with pandas,
  numpy and networkx. Module names and Slurm settings are set in the scripts.
- A JGI API key in `~/.jgi_api_key` for downloads.
- Local settings (Slurm account, BUSCO datasets folder): `cp config.example.sh config.sh`, edit it,
  and run `source config.sh` before `sbatch`.

Run all scripts from the project folder. Each script's header describes its arguments and outputs.

## Pipeline

| Step | Command | Output |
|---|---|---|
| 1. List published genomes | `scripts/list_mycocosm_genomes.sh` | `00_inputs/mycocosm_published.tsv` |
| 2. Download and decompress proteomes | `sbatch scripts/download_mycocosm_genomes.sbatch` | `01_genomes/`, `03_BUSCO/proteomes/` |
| 3. Assembly statistics | `scripts/assembly_stats.sh` | `02_assembly_stats/assembly_stats.tsv` |
| 4. BUSCO | `sbatch scripts/run_busco.sbatch` | `03_BUSCO/results/proteomes/` |
| 5. Collapse isoforms | `scripts/collapse_isoforms.sh` | `04_isoforms/` |
| 6. BUSCO on collapsed proteomes | `sbatch scripts/run_busco.sbatch 04_isoforms/collapsed_proteomes` | `03_BUSCO/results/collapsed_proteomes/` |
| 7. Combine BUSCO tables | `scripts/combine_busco_tables.sh <batch,batch,...>` | `03_BUSCO/busco_full_table_<batches>.tsv` |

Scripts skip finished work, so new genomes can be added to the cohort and the steps re-run. BUSCO
needs a new batch folder for new proteomes.

## Isoform collapse

MycoCosm GFFs store each isoform as a separate gene. Models sharing ≥ 50% of the shorter model's CDS
on the same strand are merged into one locus, and the longest protein is kept (ties: alphabetical
gene name, then lowest proteinId). Proteins are taken from the JGI protein FASTA, never translated.
On *Neurospora crassa* (`Neucr2`) this reproduces the 9,730 Broad loci exactly.

`04_isoforms/` contains `summary.tsv` (counts per genome), `removed_models.tsv` (removed → kept
model) and `unmerged_overlaps.tsv` (overlaps below 50%).

## MycoCosm notes

- FASTA files start with `;` comment lines: removed before SeqKit and BUSCO.
- `/` in protein headers is replaced by `_` (BUSCO rejects it).
- Gene names can contain `#`, so only lines starting with `#` are treated as GFF comments.
- Repeat-masked assemblies are soft-masked.
- Some GFF models have no protein (≤ 2% per genome, mostly in-frame stop codons); they are excluded.
