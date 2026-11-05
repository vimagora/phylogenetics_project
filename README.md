# phylogenetics_project

A pipeline for the phylogenetic study of fungal gene families across a cohort of
[JGI MycoCosm](https://mycocosm.jgi.doe.gov) genomes. It has two halves:
1. **Genome assessment** (steps 1–7): downloads candidate genomes, computes assembly statistics,
   reduces each proteome to one protein per gene and runs BUSCO before and after, as the basis for
   choosing a cohort.
2. **Analysis** (steps 8–18 and F1–F11): for the chosen genomes, clean protein IDs, a species tree from
   BUSCO markers (coalescent and concatenated), OrthoFinder orthogroups, and per gene family:
   anchoring of characterized proteins, family members and orthogroups, and family trees.

The repository holds the inputs and results of our study (150 genomes, sugar transporters), which can
be reproduced or replaced by your own. How our cohort was chosen is described in a separate report.

## Inputs

**Environment**, once per machine:
- `config.sh`: `cp config.example.sh config.sh`, fill in the Slurm account, the BUSCO datasets folder
  and the OrthoFinder container, and run `source config.sh` before `sbatch`;
- a JGI API key in `~/.jgi_api_key`, for downloads;
- the BUSCO `fungi_odb12` dataset, downloaded once for offline runs;
- the software under [Requirements](#requirements).

**Genomes**, in `00_inputs/`:
- `candidates.tsv`: the genomes to assess in steps 2–7: `portal`, `ncbi_taxid`, `mycocosm_name`,
  `notes`. Choose them from `mycocosm_published.tsv` (step 1), which gives the NCBI taxon ID and name.
- `genomes.tsv`: the genomes chosen for the analysis, each a candidate already processed by steps 2–7:
  `portal`, `ncbi_taxid`, `species`, `phylum`, `subphylum`, `class`, `order`, `role`, `notes`. `role` is
  `cohort` (species tree, OrthoFinder, family members) or `shadow` (second strains, haplotypes or older
  annotations, cleaned for presence checks only). Species and ranks are required for the cohort.

Keeping two tables lets you assess and filter candidates freely, and pass on only finished genomes:
`check_genomes.py` (also run by step 8) stops within seconds if a genome in `genomes.tsv` is not a
candidate, has no collapsed proteome or BUSCO results, or a cohort genome lacks its taxonomy.

**Gene families**, one folder each in `families/<family>/` (see [Gene families](#gene-families)):
`family.conf`, `references.tsv`, `outgroups.tsv`, and `orphans.tsv` if the family has orphans;
optionally `function_columns.tsv`.

## Requirements

- A Slurm cluster with BUSCO 6.1.0 (`fungi_odb12`), SeqKit 2.13.0, MAFFT 7.525, ClipKIT 2.12.1,
  IQ-TREE 3.1.4, ASTER 1.25 (ASTRAL-IV), DIAMOND 2.2.5, HMMER 3.4, FastTree (a module, or the copy in
  the OrthoFinder container) and Python 3.12 with pandas, numpy and networkx. Module names (CSC Roihu)
  and Slurm settings are set in the scripts.
- OrthoFinder 3.1.5, installed from `envs/of3_env.yaml` (see [OrthoFinder](#orthofinder)).

Run all scripts from the project folder. Each script's header describes its arguments and outputs.

## Pipeline

Genome assessment, on `candidates.tsv`:

| Step | Command | Output |
|---|---|---|
| 1. List published genomes | `scripts/list_mycocosm_genomes.sh` | `00_inputs/mycocosm_published.tsv` |
| 2. Download and decompress proteomes | `sbatch scripts/download_mycocosm_genomes.sbatch` | `01_genomes/`, `04_BUSCO/proteomes/` |
| 3. Assembly statistics | `scripts/assembly_stats.sh` | `02_assembly_stats/assembly_stats.tsv` |
| 4. BUSCO on the original proteomes | `sbatch scripts/run_busco.sbatch` | `04_BUSCO/results/proteomes/` |
| 5. Collapse isoforms | `scripts/collapse_isoforms.sh` | `03_isoforms/` |
| 6. BUSCO on collapsed proteomes | `sbatch scripts/run_busco.sbatch 03_isoforms/collapsed_proteomes` | `04_BUSCO/results/collapsed_proteomes/` |
| 7. Combine BUSCO tables | `scripts/combine_busco_tables.sh collapsed_proteomes` (any `<batch,batch,...>`) | `04_BUSCO/busco_full_table_collapsed_proteomes.tsv` |

Step 4 is not needed by later steps, but comparing it with step 6 shows how much duplication the
collapse removed, which matters when choosing proteomes.

Analysis, on `genomes.tsv`:

| Step | Command | Output |
|---|---|---|
| 8. Clean protein IDs (cohort) | `python3 scripts/clean_proteomes.py` | `05_clean_proteomes/cohort/` |
| 9. Clean protein IDs (shadows) | `python3 scripts/clean_proteomes.py --role shadow` | `05_clean_proteomes/shadows/` |
| 10. Species-tree markers | `python3 scripts/select_busco_markers.py` | `06_species_tree/` |
| 11. Marker FASTAs | `python3 scripts/extract_marker_sequences.py` | `06_species_tree/markers/raw/` |
| 12. Align and trim markers | `sbatch scripts/align_markers.sbatch` | `06_species_tree/markers/aligned/`, `trimmed/` |
| 13. Supermatrix | `python3 scripts/concatenate_markers.py` | `06_species_tree/concatenated/` |
| 14. Marker trees | `sbatch scripts/marker_trees.sbatch` | `06_species_tree/marker_trees/` |
| 15. Collect trees and models | `python3 scripts/collect_marker_trees.py` | `06_species_tree/marker_trees.nwk`, `marker_models.tsv`, `concatenated/partitions_models.nex` |
| 16. Species trees | `sbatch scripts/species_trees.sbatch [--outgroup PORTAL]` | `06_species_tree/astral/`, `concatenated/iqtree/` |
| 16b. iTOL taxonomy annotation | `python3 scripts/make_itol_taxonomy.py [tree]` | `06_species_tree/itol/` |
| 17. Species tree for OrthoFinder | `python3 scripts/prepare_species_tree.py --outgroup Catan2` | `07_orthofinder/species_tree_rooted.nwk` |
| 18. OrthoFinder | `sbatch scripts/run_orthofinder.sbatch` | `07_orthofinder/run/Results_cohort/` |

`python3 scripts/check_genomes.py` runs the checks of step 8 on their own, e.g. while choosing the
cohort.

Gene families (`<family>` is a folder in `families/`):

| Step | Command | Output (in the family's `OUT_DIR`, e.g. `09_families/<family>/`) |
|---|---|---|
| F1. Searches | `sbatch scripts/search_family.sbatch <family>` | `search/` |
| F2. Anchors, orphans, members, orthogroups | `python3 scripts/summarize_family.py <family>` | `anchors.tsv`, `members.tsv`, `orthogroup_summary.tsv`, … |
| F3. Outgroup sequences | `python3 scripts/fetch_outgroups.py <family>` | `outgroups.faa` |
| F4. Tree sequences | `python3 scripts/build_family_fasta.py <family>` | `tree/sequences.faa`, `tree/tips.tsv` |
| F5. Alignment | `sbatch scripts/align_family.sbatch <family>` | `tree/aligned.faa`, `tree/trimmed.faa` |
| F6. Quick tree | `sbatch --time=04:00:00 scripts/family_tree.sbatch <family> fasttree` | `tree/fasttree.nwk` |
| F7. iTOL files and figure | `python3 scripts/make_family_itol.py <family>`; `python3 scripts/make_family_figure.py <family>` | `itol/`; `figures/<family>_fasttree.svg` |
| ⏸ | Review of the quick tree against the criteria you set for the family | |
| F8. Model selection | `sbatch scripts/family_tree.sbatch <family> modelfinder` | `tree/modelfinder.*` |
| F9. Tree | `sbatch scripts/family_tree.sbatch <family> tree` | `tree/ml.*` |
| ⏸ | Review of the tree (F7 again for the iTOL files) | |
| F10. Subfamilies | Choose them on the trees and record them in `families/<family>/subfamilies.tsv` (by hand) | |
| F11. Subfamily tip lists | `python3 scripts/build_subfamilies.py <family>` | `subfamilies/report.md`, `subfamilies/<subfamily>/tips.tsv`, `sequences.faa` |
| ⏸ | Review of the report before aligning | |

Steps 8–11, 13, 15 and 17 use Python from `module load python-data`. Step 10 reads the combined table of
the collapsed proteomes from step 7.

The download, extraction, collapse, BUSCO and alignment steps skip finished work, so new candidates can
be added and the steps re-run; BUSCO needs a new batch folder for new proteomes. Steps 8–11, 13 and 15
are fast and rewrite their outputs on every run. The tree jobs (14, 16, F8, F9) resume when resubmitted.

`70_structure/` holds a concept test on protein structures of sugar transporters, which is not part of
the pipeline (see `70_structure/README.md`).

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

Both are rooted on the outgroup given with `--outgroup` (default *Catenaria anguillulae*, `Catan2`, so
that Blastocladiomycota diverges first, as in the concatenated tree of Li et al. 2021); ASTRAL needs the
root to estimate branch lengths accurately. A finished ASTRAL tree and IQ-TREE run are not repeated,
so to rerun only ASTRAL, move `06_species_tree/astral/` away and resubmit with a short time limit
(`sbatch --time=01:00:00 scripts/species_trees.sbatch`).

## OrthoFinder

OrthoFinder 3 runs on the cohort proteomes in one standard run (all-vs-all DIAMOND search, MSA gene
trees with FAMSA and FastTree), with the rooted ASTRAL species tree (`-s`, rerooted on `Catan2` by
step 17) so that the hierarchical orthogroups follow the species tree. `-X` keeps the
`<portal>-<proteinId>` gene IDs. `Orthogroups/Orthogroups.tsv` holds the hierarchical orthogroups (the
former `N0.tsv`), and `Orthogroups.GeneCount.tsv` the copy numbers. After the run, the folders with one
file per orthogroup or species pair are packed into `.tar.gz` archives (CSC scratch has a file-count
quota; a run on 150 genomes peaks at about 700,000 files).

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

## Gene families

Steps F1–F11 work on any gene family described by a folder `families/<family>/`. The scripts never
transfer a function to an uncharacterized protein; they stop at alignments, trees and annotation
files, and the reviews between them are done by a person.

**Files of a family**
- `family.conf`: `NAME`, `HMMS` (Pfam accessions defining the family, space-separated),
  `MIN_HMM_COVERAGE`, `OUT_DIR` (results folder) and, optionally, `MARK_PORTALS` (portals whose proteins
  get a tick in the figure, space-separated).
- `references.tsv`: the reference proteins, one row each: `reference_id` (lowercase letters, digits,
  `_`), `evidence` (`characterized`, `genetic` for a deletion phenotype without an assay, or `predicted`), `species`, `short_name` (species abbreviation,
  e.g. `Dekbr`; used as the pseudo-portal of orphans), `protein_name`, `role`, `function`
  (`;`-separated), `cohort_portal` (the protein's own cohort genome, if known), `cohort_id` (its clean
  ID, if known), `use` (`TRUE`/`FALSE`), `references`, `notes`, `sequence`.
- `outgroups.tsv`: `outgroup_id` (`<Abbrev>-<protein>`), `group`, `organism`, `protein`, `uniprot`,
  `root` (`TRUE` for the outgroups the tree is rooted on), `notes`.
- `orphans.tsv`: one row per orphan pseudo-portal: `pseudo_portal`, `species`, `phylum`, `subphylum`,
  `class`, `order`, `nearest_cohort_portal`, `notes`. F2 lists pseudo-portals still missing here.
- `function_columns.tsv` (optional): `label`, `colour`, `terms`: the iTOL function columns, each ticked
  by any of its terms (lowercase, comma-separated) in a characterized reference's `function`.
- `subfamilies.tsv` (F10): one row per set of tips: `subfamily`, `set` (`ingroup`, whose rows add up, or
  the name of an outgroup lineage), `kind` (`clade`, `orthogroups`, `tips` or `uniprot`), `tree` (`ml` or
  `fasttree`, for `clade`), `definition`, `sample` (`all`, or the number of tips to take from a clade) and
  `notes`. A `clade` is the most recent common ancestor of the tips in `definition`, in that tree rooted
  on the family's root outgroups; a `uniprot` item is `<tip_id>=<accession>`.

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
1000 ultrafast bootstraps and 1000 SH-aLRT replicates, rooted on the outgroups marked `root`. IQ-TREE
chooses its number of threads (`-T AUTO`), since trimmed family alignments are short.
`make_family_itol.py` shows function only for characterized references.

**Figure (F7).** `make_family_figure.py` draws the tree as a circular SVG with editable text, using only
the Python standard library. From the centre outward: the tree; a ring coloured by orthogroup (those
with at least 1% of the tips, in the same colours as the iTOL files); for each family HMM, a ring of the
full-sequence bit score (0 to the 99th percentile) and a ring of the HMM's coverage (with the membership
threshold on its colour bar); ticks for `MARK_PORTALS`; and the names of `characterized` and `genetic`
references (black) and outgroups (red), moved apart along the circle as little as possible. The
full-sequence score is used rather than the best or mean domain score because it adds up the pieces of
a split domain; coverage counts each HMM position once. Read together, low score with low coverage is a
fragment, low score with high coverage a poor match to the model, and a score above the usual range
with full coverage a possible tandem duplicate. The tree is rooted on the outgroups with `root` = TRUE,
or at its midpoint if there are none or they do not form a clade, as the caption states. `--tree ml`
draws the IQ-TREE tree and `--out` writes elsewhere.

**Subfamilies (F10–F11).** Large families are split into subfamily trees: the whole-family tree shows
where to cut, and the subfamily trees, each with its own alignment, are the ones interpreted. A lineage
outside a supported ingroup roots it whatever its exact position, so `build_subfamilies.py` samples
each outgroup lineage: it cuts the clade into as many parts as tips wanted and takes from each part
the characterized protein, else the member with HMM coverage ≥ 0.8, with the shortest terminal branch.
A UniProt protein identical to a family tip is replaced by that tip. The report lists each subfamily's
sets, orthogroups, anchors and outgroups, and the tips in no subfamily. `--lists-only` writes the lists
and report without sequences.

## MycoCosm notes

- FASTA files start with `;` comment lines: removed before SeqKit and BUSCO.
- `/` in protein headers is replaced by `_` (BUSCO rejects it).
- Gene names can contain `#`, so only lines starting with `#` are treated as GFF comments.
- Repeat-masked assemblies are soft-masked.
- Some GFF models have no protein (≤ 2% per genome, mostly in-frame stop codons); they are excluded.
- The published-genomes catalog comes in Windows-1252; step 1 stores it as UTF-8.
