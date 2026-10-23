# phylogenetics_project

Downloads fungal genomes from [JGI MycoCosm](https://mycocosm.jgi.doe.gov), checks assembly and
annotation quality, builds one clean proteome per genome (one protein per gene) and prepares the
BUSCO markers for the species tree.

Genome lists (TSV, header row, portal IDs in column 1) in `00_inputs/`:
- `mycocosm_datasets.tsv`: all 186 candidates, used by steps 2–7;
- `cohort_datasets.tsv`: the 150 selected genomes, which enter the species tree and OrthoFinder;
- `cohort_taxonomy.tsv`: phylum, subphylum, class and order of each cohort genome (from `cohort_selection.md`);
- `shadows_datasets.tsv`: 16 second strains, haplotypes or older annotations, only searched for presence.

How the cohort was chosen is described in a separate report (not part of the repository).

## Requirements

- A Slurm cluster with BUSCO 6.1.0 (`fungi_odb12`), SeqKit 2.13.0, MAFFT 7.525, ClipKIT 2.12.1,
  IQ-TREE 3.1.4, ASTER 1.25 (ASTRAL-IV), DIAMOND 2.2.5, HMMER 3.4, FastTree (a module, or the copy in
  the OrthoFinder container) and Python 3.12 with pandas, numpy and networkx. Module names (CSC Roihu)
  and Slurm settings are set in the scripts.
- OrthoFinder 3.1.5, installed from `envs/of3_env.yaml` (see [OrthoFinder](#orthofinder)).
- A JGI API key in `~/.jgi_api_key` for downloads.
- Local settings (Slurm account, BUSCO datasets folder, OrthoFinder container): `cp config.example.sh config.sh`, edit it,
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
| 16b. iTOL taxonomy annotation | `python3 scripts/make_itol_taxonomy.py [tree]` | `06_species_tree/itol/` |
| 17. Species tree for OrthoFinder | `python3 scripts/prepare_species_tree.py --outgroup Catan2` | `07_orthofinder/species_tree_rooted.nwk` |
| 18. OrthoFinder | `sbatch scripts/run_orthofinder.sbatch` | `07_orthofinder/run/Results_cohort/` |
| 19. AlphaFold3 inputs | `python3 scripts/make_af3_inputs.py` | `08_structures/af3_inputs/`, `structure_ids.tsv` |
| 20. AlphaFold3 structures | `sbatch [--partition=… --time=…] scripts/run_af3.sbatch` | `08_structures/af3_outputs/` |
| 21. AlphaFold DB models (concept test) | `python3 scripts/fetch_afdb_models.py` | `08_structures/models_trimmed/`, `models.tsv` |
| 22. Structure trees (Foldtree-style) | `python3 scripts/structure_trees.py` | `08_structures/trees/structure_*.nwk` |
| 23. Sequence tree (baseline) | `scripts/sequence_tree.sh` | `08_structures/trees/sequence_ml.treefile` |
| 24. Substrate clustering scores | `python3 scripts/evaluate_substrate_clustering.py` | `08_structures/evaluation/` |
| 25. Pocket mapping (concept test) | `python3 scripts/map_pockets.py` | `08_structures/pockets/` |
| 26. Pocket scores | `python3 scripts/evaluate_pockets.py` | `08_structures/pockets/evaluation/`, `trees/` |

Gene families (generic; `<family>` is a folder in `families/`, see [Gene families](#gene-families)):

| Step | Command | Output (in the family's `OUT_DIR`, e.g. `09_families/<family>/`) |
|---|---|---|
| F1. Searches (Roihu) | `sbatch scripts/search_family.sbatch <family>` | `search/` |
| F2. Anchors, orphans, members, orthogroups | `python3 scripts/summarize_family.py <family>` | `anchors.tsv`, `members.tsv`, `orthogroup_summary.tsv`, … |
| F3. Outgroup sequences | `python3 scripts/fetch_outgroups.py <family>` | `outgroups.faa` |
| F4. Tree sequences (Roihu) | `python3 scripts/build_family_fasta.py <family>` | `tree/sequences.faa`, `tree/tips.tsv` |
| F5. Alignment (Roihu) | `sbatch scripts/align_family.sbatch <family>` | `tree/aligned.faa`, `tree/trimmed.faa` |
| F6. Quick tree (Roihu) | `sbatch --time=04:00:00 scripts/family_tree.sbatch <family> fasttree` | `tree/fasttree.nwk` |
| F7. iTOL files | `python3 scripts/make_family_itol.py <family>` | `itol/` |
| ⏸ | Review of the quick tree against `families/<family>/expectations.md` | |
| F8. Model selection (Roihu) | `sbatch scripts/family_tree.sbatch <family> modelfinder` | `tree/modelfinder.*` |
| F9. Tree (Roihu) | `sbatch scripts/family_tree.sbatch <family> tree` | `tree/ml.*` |
| ⏸ | Review of the tree (F7 again for the iTOL files) | |

Family-specific step for the sugar transporters: `python3 families/sugar_transporters/make_references.py`
converts the two curated tables in `00_inputs` into `families/sugar_transporters/references.tsv`.

Steps 8–11, 13, 15 and 17 use Python from `module load python-data`. Step 10 reads the combined table of
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

Both are rooted on *Catenaria anguillulae* (`Catan2`; `--outgroup` to change it), so that
Blastocladiomycota diverges first, as in the concatenated tree of Li et al. (2021); ASTRAL needs the
root to estimate branch lengths accurately. The first runs used `Spipu1`, which falls inside
Chytridiomyceta. A finished ASTRAL tree and IQ-TREE run are not repeated, so to rerun only ASTRAL,
move `06_species_tree/astral/` away and resubmit with a short time limit
(`sbatch --time=01:00:00 scripts/species_trees.sbatch`).

## OrthoFinder

OrthoFinder 3 runs on the 150 cohort proteomes in one standard run (all-vs-all DIAMOND search, MSA
gene trees with FAMSA and FastTree), with the ASTRAL species tree (`-s`) so that the
hierarchical orthogroups follow our species tree. For OrthoFinder, the tree is rerooted on
*Catenaria anguillulae* (`Catan2`), so that Blastocladiomycota diverges first, as in the
concatenated tree of Li et al. (2021); the `Spipu1` root falls inside Chytridiomyceta. `-X` keeps the `<portal>-<proteinId>` gene IDs.
`Orthogroups/Orthogroups.tsv` holds the hierarchical orthogroups (the former `N0.tsv`), and
`Orthogroups.GeneCount.tsv` the copy numbers. After the run, the folders with one file per
orthogroup or species pair are packed into `.tar.gz` archives (CSC scratch has a file-count quota).

Installation on Roihu as a Tykky container, once:
```bash
module load tykky
conda-containerize new --prefix /projappl/<project>/orthofinder envs/of3_env.yaml
```
Then set `ORTHOFINDER_DIR=/projappl/<project>/orthofinder` in `config.sh`. DIAMOND is pinned below
2.2, as OrthoFinder 3.1.5 requires.
The bioconda launcher loads its code from a `src/` folder next to itself, which the Tykky wrapper
hides; `run_orthofinder.sbatch` therefore sets `PYTHONPATH` to that folder inside the container
(mount point read from the container's `common.sh`).

OrthoFinder refuses an existing output folder and cannot resume. If a run stops after the DIAMOND
search, restart from its search results with `orthofinder -b <Results>/WorkingDirectory` (to rerun
a finished run with another species tree, unpack `WorkingDirectory.tar.gz` first).

## Sugar transporter structures

`00_inputs/characterized_sugar_transporters.tsv` lists characterized fungal sugar transporters (from
Xu et al. 2024, Table S1, plus later literature and TCDB family 2.A.1.1), with substrates, references,
a `role` column (transporter or sensor) and `use_in_structure_analysis` (FALSE for non-MFS proteins and
partial sequences). `00_inputs/dicsqu_sugar_transporters.tsv` adds the *D. squalens* transporters with
predicted specificities. AlphaFold3 runs in CSC's container on a GH200 GPU node, with your own model
parameters (`AF3_MODEL_DIR` in `config.sh`; the licence forbids sharing them). `run_af3.sbatch` skips
proteins that already have a model, so it can be resubmitted until all are done.

Concept test (steps 21–24, environment `envs/structure_test.yaml`, runs on a laptop): AlphaFold DB
models (identical sequence, or ≥ 98% identity in the same species) are trimmed of termini with
pLDDT < 70, compared all-vs-all with Foldseek, and turned into neighbour-joining trees (fident,
alntmscore, lddt), next to a sequence ML tree of the same proteins. Characterized transporters are
scored by (a) whether the nearest neighbour shares a substrate (p-value from shuffled labels) and
(b) the purity of the smallest clade holding all transporters of a substrate.

Pocket test (steps 25–26, same environment): the sugar-binding pocket of three sugar-porter crystal
structures with a bound sugar (XylE, PDB 4GBY and 4GBZ; GLUT3, 4ZW9; STP10, 6H7D), i.e. residues within
5 Å (core) or 8 Å (extended) of the sugar, is mapped onto every protein through a MAFFT alignment with
the reference sequences and, as a check, through Foldseek alignments of the models to the
references; columns are kept where both mappings agree for ≥ 80% of the models. Pocket signatures
are compared by BLOSUM62 distance and scored by (a), by (a') the nearest neighbour among transporters
< 40% identical, and by (c) the most common substrate of the 3 nearest, against the sequence tree
(paired sign-flip test) and random sets of well-aligned columns of the pocket's size.

## Gene families

Steps F1–F9 work on any gene family described by a folder `families/<family>/`. The scripts never
transfer a function to an uncharacterized protein; they stop at alignments, trees and annotation
files, and the reviews between them are done by a person.

**Files of a family**
- `family.conf`: `NAME`, `HMMS` (Pfam accessions defining the family, space-separated),
  `MIN_HMM_COVERAGE` and `OUT_DIR` (results folder).
- `references.tsv`: the reference proteins, one row each: `reference_id` (lowercase letters, digits,
  `_`), `evidence` (`characterized` or `predicted`), `species`, `short_name` (species abbreviation,
  e.g. `Dekbr`; used as the pseudo-portal of orphans), `protein_name`, `role`, `function`
  (`;`-separated), `cohort_portal` (the protein's own cohort genome, if known), `cohort_id` (its clean
  ID, if known), `use` (`TRUE`/`FALSE`), `references`, `notes`, `sequence`.
- `orphans.tsv`: one row per orphan pseudo-portal: `pseudo_portal`, `species`, `phylum`, `subphylum`,
  `class`, `order`, `nearest_cohort_portal`, `notes`. F2 lists pseudo-portals still missing here.
- `outgroups.tsv`: `outgroup_id` (`<Abbrev>-<protein>`), `group`, `organism`, `protein`, `uniprot`,
  `root` (`TRUE` for the outgroups the tree is rooted on), `notes`.
- `expectations.md`: what the outgroup and the characterized subfamilies should look like in the tree,
  written before any tree is built and checked by eye at the reviews.
- `function_columns.tsv` (optional): `label`, `colour`, `terms`: the iTOL function columns, each ticked
  by any of its terms (lowercase, comma-separated) in a characterized reference's `function`.

**Anchors and orphans (F2).** The references are searched against all cohort proteins with DIAMOND
(`--ultra-sensitive`). The first rule that applies gives each reference its status:
- `anchor_listed_id`: its `cohort_id`, at ≥ 98% identity over ≥ 90% of the reference;
- `anchor`: in its own genome (`cohort_portal` if that is a cohort genome, else the cohort genome of its
  species), ≥ 98% identity over ≥ 90% of both proteins;
- `anchor_gene_model_differs`: the same over ≥ 90% of the shorter protein only;
- `probable_anchor`: in its own genome, 95–98% identity over ≥ 90% of the shorter protein, and the
  cohort protein's closest reference is this one;
- orphan (`orphan_species_not_in_cohort`, `orphan_no_copy_in_genome`): enters the tree with its own
  sequence as `<short_name>-<protein_name>`. The species-only pseudo-portal (`Dekbr`) cannot clash with
  MycoCosm portals, which carry a version or strain (`Dekbr2`), and `^Dekbr` finds both. Orphans are
  labels only: they have no tip in the species tree and are left out of counts and reconciliation.

**Members (F2).** Cohort proteins whose domains of one family HMM cover ≥ `MIN_HMM_COVERAGE` of the HMM,
at Pfam's gathering threshold (`hmmsearch --cut_ga`). There is deliberately no score filter beyond
that (e.g. the lowest score of the characterized proteins), so that unusual members are kept; hits
below the coverage are listed in `hmm_hits.tsv`. Orthogroups are OrthoFinder's `Orthogroups.tsv`.

**Trees (F4–F9).** Members (full sequences from the clean proteomes; shadows never enter), orphans
that pass the membership rule themselves, and outgroups are aligned with MAFFT `--auto` and trimmed
with ClipKIT gappy 0.5, as for the species tree. Extra domains (e.g. from fusions or misannotations)
are not removed beforehand; their columns are mostly gaps and are trimmed. FastTree gives a quick
tree for the first review; the final tree is IQ-TREE 3 with ModelFinder on the full alignment, then
1000 ultrafast bootstraps and 1000 SH-aLRT replicates, rooted on the outgroups marked `root`.
`make_family_itol.py` shows function only for characterized references.

## MycoCosm notes

- FASTA files start with `;` comment lines: removed before SeqKit and BUSCO.
- `/` in protein headers is replaced by `_` (BUSCO rejects it).
- Gene names can contain `#`, so only lines starting with `#` are treated as GFF comments.
- Repeat-masked assemblies are soft-masked.
- Some GFF models have no protein (≤ 2% per genome, mostly in-frame stop codons); they are excluded.
