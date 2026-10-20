# phylogenetics_project

Downloads fungal genomes from [JGI MycoCosm](https://mycocosm.jgi.doe.gov), checks assembly and
annotation quality, builds one clean proteome per genome (one protein per gene) and prepares the
BUSCO markers for the species tree.

Genome lists (TSV, header row, portal IDs in column 1) in `00_inputs/`:
- `mycocosm_datasets.tsv`: all 186 candidates, used by steps 2–7;
- `cohort_datasets.tsv`: the 150 selected genomes, which enter the species tree and OrthoFinder;
- `shadows_datasets.tsv`: 16 second strains, haplotypes or older annotations, only searched for presence.

How the cohort was chosen is described in `00_inputs/cohort_selection.md`.

## Requirements

- A Slurm cluster with BUSCO 6.1.0 (`fungi_odb12`), SeqKit 2.13.0, MAFFT 7.525, ClipKIT 2.12.1,
  IQ-TREE 3.1.4, ASTER 1.25 (ASTRAL-IV) and Python 3.12 with pandas, numpy and networkx. Module names (CSC Roihu) and Slurm settings are
  set in the scripts.
- A JGI API key in `~/.jgi_api_key` for downloads.
- Local settings (Slurm account, BUSCO datasets folder): `cp config.example.sh config.sh`, edit it,
  and run `source config.sh` before `sbatch`.

Run all scripts from the project folder. Each script's header describes its arguments and outputs.

## Pipeline

| Step | Command | Output |
|---|---|---|
| 1. List published genomes | `scripts/list_mycocosm_genomes.sh` | `00_inputs/mycocosm_published.tsv` |
| 2. Download and decompress proteomes | `sbatch scripts/download_mycocosm_genomes.sbatch` | `01_genomes/`, `04_BUSCO/proteomes/` |
| 3. Assembly statistics | `scripts/assembly_stats.sh` | `02_assembly_stats/assembly_stats.tsv` |
| 4. BUSCO | `sbatch scripts/run_busco.sbatch` | `04_BUSCO/results/proteomes/` |
| 5. Collapse isoforms | `scripts/collapse_isoforms.sh` | `03_isoforms/` |
| 6. BUSCO on collapsed proteomes | `sbatch scripts/run_busco.sbatch 03_isoforms/collapsed_proteomes` | `04_BUSCO/results/collapsed_proteomes/` |
| 7. Combine BUSCO tables | `scripts/combine_busco_tables.sh collapsed_proteomes` (any `<batch,batch,...>`) | `04_BUSCO/busco_full_table_collapsed_proteomes.tsv` |
| 8. Clean protein IDs (cohort) | `python3 scripts/clean_proteomes.py` | `05_clean_proteomes/cohort/` |
| 9. Clean protein IDs (shadows) | `python3 scripts/clean_proteomes.py 00_inputs/shadows_datasets.tsv 03_isoforms/collapsed_proteomes 05_clean_proteomes/shadows` | `05_clean_proteomes/shadows/` |
| 10. Species-tree markers | `python3 scripts/select_busco_markers.py` | `06_species_tree/` |
| 11. Marker FASTAs | `python3 scripts/extract_marker_sequences.py` | `06_species_tree/markers/raw/` |
| 12. Align and trim markers | `sbatch scripts/align_markers.sbatch` | `06_species_tree/markers/aligned/`, `trimmed/` |
| 13. Supermatrix | `python3 scripts/concatenate_markers.py` | `06_species_tree/concatenated/` |
| 14. Marker trees | `sbatch scripts/marker_trees.sbatch` | `06_species_tree/marker_trees/` |
| 15. Collect trees and models | `python3 scripts/collect_marker_trees.py` | `06_species_tree/marker_trees.nwk`, `marker_models.tsv`, `concatenated/partitions_models.nex` |
| 16. Species trees | `sbatch scripts/species_trees.sbatch` | `06_species_tree/astral/`, `concatenated/iqtree/` |

Steps 8–11, 13 and 15 use Python from `module load python-data`. Step 10 reads the combined table of
the collapsed proteomes from step 7.

The download, extraction, collapse, BUSCO and alignment steps skip finished work, so new genomes can
be added and the steps re-run; BUSCO needs a new batch folder for new proteomes. Steps 8–11, 13 and 15
are fast and rewrite their outputs on every run. The tree jobs (14, 16) resume when resubmitted.

## Isoform collapse

MycoCosm GFFs store each isoform as a separate gene. Models sharing ≥ 50% of the shorter model's CDS
on the same strand are merged into one locus, and the longest protein is kept (ties: alphabetical
gene name, then lowest proteinId). Proteins are taken from the JGI protein FASTA, never translated.
On *Neurospora crassa* (`Neucr2`) this reproduces the 9,730 Broad loci exactly.

`03_isoforms/` contains `summary.tsv` (counts per genome), `removed_models.tsv` (removed → kept
model) and `unmerged_overlaps.tsv` (overlaps below 50%).

## Protein IDs

From step 8 on, proteins are named `<portal>-<proteinId>` (e.g. `Dicsqu464_2-3083339`): MycoCosm
headers contain `|`, and some gene names `#`, which break Newick trees. Portal IDs contain only
letters, digits and `_`, so the last `-` separates portal and proteinId. `05_clean_proteomes/*/id_map.tsv`
links each clean ID to its gene name and original header, for anchors and literature IDs. The
terminal `*` is removed from every sequence. BUSCO results keep the original headers;
`select_busco_markers.py` adds the clean ID to each selected sequence.

## Species-tree markers

BUSCO markers complete and single-copy in ≥ 90% of the cohort, and duplicated in ≤ 5% of the genomes
with ≤ 6% duplicated markers overall (so whole-genome duplications do not reject markers). In each
genome only a complete single-copy hit is used. `06_species_tree/genome_coverage.tsv` lists the
selected markers found per genome.

Each marker is aligned with MAFFT L-INS-i and trimmed with ClipKIT in gappy mode
(columns with > 50% gaps removed; the smart-gap mode kept most columns of these kingdom-wide alignments). Sequences are named
by portal. Markers with ≥ 100 columns after trimming and sequences from ≥ 90% of the cohort are
concatenated into `06_species_tree/concatenated/supermatrix.faa`, with one partition per marker
(`partitions.nex`); `marker_summary.tsv` and `taxon_occupancy.tsv` report what was kept.

## Species trees

Each kept marker gets a maximum-likelihood tree (IQ-TREE 3, ModelFinder among nuclear models,
1000 ultrafast bootstraps with NNI refinement). Two species trees are inferred:
- **coalescent:** ASTRAL-IV on the marker trees (`06_species_tree/astral/species_tree_astral.nwk`);
- **concatenated:** IQ-TREE 3 on the supermatrix, one partition per marker with that marker's best-fit
  model and edge-linked proportional branch lengths, with SH-aLRT and ultrafast bootstrap support
  (`06_species_tree/concatenated/iqtree/supermatrix.treefile`).

Both are rooted on *Spizellomyces punctatus* (`Spipu1`); ASTRAL needs the root to estimate branch
lengths accurately.

## MycoCosm notes

- FASTA files start with `;` comment lines: removed before SeqKit and BUSCO.
- `/` in protein headers is replaced by `_` (BUSCO rejects it).
- Gene names can contain `#`, so only lines starting with `#` are treated as GFF comments.
- Repeat-masked assemblies are soft-masked.
- Some GFF models have no protein (≤ 2% per genome, mostly in-frame stop codons); they are excluded.
