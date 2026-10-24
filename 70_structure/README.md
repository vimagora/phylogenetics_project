# Structure concept test (sugar transporters)

A test of whether protein structures group fungal sugar transporters by substrate better than
sequence does. It is not part of the pipeline: the result was negative, and the sugar transporter trees
of the pipeline (`families/sugar_transporters/`) are built from sequence. The scripts are specific to
the sugar transporters and read the family's reference table,
`families/sugar_transporters/references.tsv`.

Environment: `70_structure/envs/structure_test.yaml` (runs on a laptop), except AlphaFold3 (step S2),
which runs in CSC's container on a GH200 GPU node with your own model parameters (`AF3_MODEL_DIR`,
`AF3_SIF` and `AF3_DB_DIR` in `config.sh`; the licence forbids sharing the parameters).

| Step | Command | Output |
|---|---|---|
| S1. AlphaFold3 inputs | `python3 70_structure/scripts/make_af3_inputs.py` | `70_structure/af3_inputs/`, `structure_ids.tsv` |
| S2. AlphaFold3 structures | `sbatch [--partition=… --time=…] 70_structure/scripts/run_af3.sbatch` | `70_structure/af3_outputs/` |
| S3. AlphaFold DB models | `python3 70_structure/scripts/fetch_afdb_models.py` | `70_structure/models_trimmed/`, `models.tsv` |
| S4. Structure trees (Foldtree-style) | `python3 70_structure/scripts/structure_trees.py` | `70_structure/trees/structure_*.nwk` |
| S5. Sequence tree (baseline) | `70_structure/scripts/sequence_tree.sh` | `70_structure/trees/sequence_ml.treefile` |
| S6. Substrate clustering scores | `python3 70_structure/scripts/evaluate_substrate_clustering.py` | `70_structure/evaluation/` |
| S7. Pocket mapping | `python3 70_structure/scripts/map_pockets.py` | `70_structure/pockets/` |
| S8. Pocket scores | `python3 70_structure/scripts/evaluate_pockets.py` | `70_structure/pockets/evaluation/`, `trees/` |

S1 also writes `structure_ids.tsv`, which the later steps read. `run_af3.sbatch` skips proteins that
already have a model, so it can be resubmitted until all are done; only a single-sequence test was run,
since the AlphaFold DB covered 199 of the 210 proteins.

**Structure trees (S3–S6).** AlphaFold DB models (identical sequence, or ≥ 98% identity in the same
species) are trimmed of termini with pLDDT < 70, compared all-vs-all with Foldseek, and turned into
neighbour-joining trees (fident, alntmscore, lddt), next to a sequence ML tree of the same proteins.
Characterized transporters are scored by (a) whether the nearest neighbour shares a substrate (p-value
from shuffled labels) and (b) the purity of the smallest clade holding all transporters of a substrate.

**Pocket test (S7–S8).** The sugar-binding pocket of three sugar-porter crystal structures with a bound
sugar (XylE, PDB 4GBY and 4GBZ; GLUT3, 4ZW9; STP10, 6H7D), i.e. residues within 5 Å (core) or 8 Å
(extended) of the sugar, is mapped onto every protein through a MAFFT alignment with the reference
sequences and, as a check, through Foldseek alignments of the models to the references; columns are
kept where both mappings agree for ≥ 80% of the models. Pocket signatures are compared by BLOSUM62
distance and scored by (a), by (a') the nearest neighbour among transporters < 40% identical, and by
(c) the most common substrate of the 3 nearest, against the sequence tree (paired sign-flip test) and
random sets of well-aligned columns of the pocket's size.
