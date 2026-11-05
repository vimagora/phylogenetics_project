# Subfamily plan

Plan for splitting the sugar transporter family into subfamily trees, written on 2026-10-29 after
comparing the FastTree tree (`review_fasttree.md`) with the IQ-TREE tree of the whole family
(`09_families/sugar_transporters/tree/ml.*`). Steps F10–F11 are written (§5); the rest is still a plan.
The generic steps are to serve every family, starting with the aromatic families.

## 1. Why subfamily trees

The whole-family trees (8,681 sequences, 455 columns after trimming) resolve clades of moderate
depth but not the backbone:

- The two trees share 79% of their splits (Robinson–Foulds 3,674 of 17,356). 43.5% of the IQ-TREE
  nodes have SH-aLRT ≥ 80 and UFBoot ≥ 95; the splits nearest the root have UFBoot 14–40.
- 455 columns are the core that ClipKIT (gappy 0.5) keeps across a family spanning bacteria to fungi.
  A subfamily aligns over much more of its length, so its own alignment keeps more columns for fewer
  sequences.
- The plant STPs and human GLUTs fall inside the fungal ITR + OG0003609 clade in both trees, and XylE
  groups with anaerobic gut fungi. The main lineages of the sugar porters are probably older than the
  split of fungi, animals and plants, so no eukaryotic protein roots the whole family. Each
  subfamily is rooted on its own (§3).

The whole-family trees are therefore a map for choosing subfamilies, and the subfamily trees are the
ones interpreted. The IQ-TREE run on the whole family (101 h) was useful once, to show which
boundaries are stable, and does not become a standard step (§5).

## 2. Subfamilies

Selection is based on the FastTree tree, because future families will have no whole-family IQ-TREE
tree. A subfamily is a FastTree clade that is mostly one orthogroup or a coherent set of orthogroups;
looks consistent in the score and coverage rings of the figure, or has a pattern we can explain; and
holds 100–2,500 tips, with priority for anchors and *D. squalens* proteins. Local support is reported
as a flag, not used as a cutoff, since it does not predict the IQ-TREE result in either direction (MAL's
group is 0.58 here and 99.6/100 there; CDT is 0.89 here and 35.6/50 there). Low support means looking
harder at that clade's own tree. The IQ-TREE run is kept as a calibration of these boundaries for this
family, not as an input.

Anchors are tips with characterized references; *D. squalens* counts are all its proteins in the clade.

| # | Subfamily | Ingroup | FastTree tips (support) | Agreement with IQ-TREE (SH-aLRT / UFBoot) | Anchors | Priority |
|---|---|---|---|---|---|---|
| 1 | `hexose` | OG0000008 piece B: hexose transporters and sensors | 690 (1.00) | same clade (99% shared), 85.4/69 | 77, incl. *A. niger* MstA, MstC, MstH, HxtB–E, sensors Snf3, Rgt2, RAG4, Hxs1, RCO-3; 4 *D. squalens* proteins | 1 |
| 2 | `mal` | MAL (OG0000174) plus the 186 OG0000008 members grouped with it | 903 (0.58) | same clade: IQ-TREE's 898 tips are all inside, 99.6/100; only the arrangement inside it varies | 10 MAL | 2 |
| 3 | `cdt` | CDT/lactose (OG0000100) | 1,080 (0.89) | same size, one tip differs, 35.6/50 | 18, incl. CDT-1, CDT-2, CltB/LacpB, LAC12 | 2 |
| 4 | `fsy` | Fsy (OG0000345) | 462 (0.88) | two clades in IQ-TREE (402, 77.7/86, and 60, 100/100); its own tree tests this | 9, incl. Fsy1 orthologs and Frt1 | 2 |
| 5 | `itr` | ITR (OG0001120) + OG0003609, with the STPs and GLUTs inside as landmarks | 390 (0.84) | same clade, 88.9/46; ITR core 100/99 | 4, incl. ITR1, ITR2, Itr2; STPs, GLUTs | 2 |
| 6 | `galacturonate_quinate` | OG0000008 piece C | 575 (0.98) | all 572 of its OG0000008 members are in IQ-TREE's 730-tip clade, which joins it with tree 8 without support (13.1/12) | 8: GatA, Gat-1, Gat1, Gat2, GalA, Qut-D, qa-y, AHA12287.1; 2 *D. squalens* proteins | 3 |
| 7 | `pentose` | OG0000008 piece A | 2,160 (0.92) | core 96% shared, 61.2/43; a 76-tip clade inside it stands apart in IQ-TREE (86.9/81) | 43, incl. XltA–C, Xut1, Xut3, Str1–2, LAT-1, AraT, HGT1; 11 *D. squalens* proteins | 3 |
| 8 | `gut_fungi` | OG0000008 clade of early-diverging fungi | 155 (1.00) | IQ-TREE joins it with tree 6 without support (13.1/12) | none; no *D. squalens* proteins | 4 |
| 9 | `weak` | all 335 members of OG0000772 + OG0000102 | by orthogroup; OG0000772 is a clade of 227 (1.00), OG0000102 of 114 (0.16), together 342 (0.33) | OG0000772 stays together; only 7 of OG0000102 join it (236 tips, 92.8/94) | 1 (Tr Tr44175, in OG0000772) | on hold (§6) |

Tree 8 comes almost entirely from anaerobic gut fungi: *Neocallimastix* (106 tips), *Anaeromyces* (21),
*Caecomyces* (12) and *Piromyces* (12). *Neocallimastix* has many duplicated genes (36.3% duplicated BUSCO
markers), so this is probably a lineage-specific expansion.

The two small OG0000008 clades first listed as a tree of their own (71 and 54 members) fall inside the
ingroup of tree 2, so they have no separate tree. The earlier `galacturonate_quinate` ingroup (730 tips)
is split into trees 6 and 8, which are not sister clades in FastTree.

**Not in any subfamily: 1,934 of 8,681 tips (22%), 12 of the 42 *D. squalens* proteins, no
characterized protein.** They are mostly OG0000018 (476), OG0000297 (464), OG0000544 (327), OG0000296
(310) and OG0000015 (115), which were missed when this plan was first written. With OG0000772 and
OG0000102 they form one lineage (2,135 tips, 94.6/45; FastTree 2,088, 0.97) whose only characterized
protein is Tr44175. See §6. The rest are small groups, the four proteins branching first in the
IQ-TREE tree (`Klula1-2178`, `Sacce1-2186`, `Zygrou1-124717`, `Lepor2-574981`) and XylE with its four
gut-fungal proteins; they stay in the whole-family trees.

The exact ingroup of each subfamily is set in a table (§5, `subfamilies.tsv`) by most recent common
ancestors of named tips in the FastTree tree, so that it can be rebuilt; the member lists it gives are
reviewed before alignment.

**Tree 9** asks whether OG0000772 + OG0000102 belong in the family: they are the weakest PF00083
members (mean HMM score 127 and 38, against 264 for all members), and OG0000102 may belong to another
MFS family. The one characterized protein in the clade, *T. reesei* Tr44175, is annotated as a member
of the monocarboxylate transporter family. Unlike the other trees, this one needs candidate families
as well as a root, so its comparison set replaces the rule of §3:

| Group | Role | Proteins (UniProt; TCDB family) |
|---|---|---|
| Sugar porters | are the weak proteins with them? | two characterized proteins from each of `hexose`, `pentose`, `mal`, `cdt` and `itr`, chosen in F10 |
| Other PF00083 hits | the nearest non-sugar families | Pho84 (P25297) and Git1 (P25346), phosphate:H+ symporters 2.A.1.9; Jen1 (P36035), 2.A.1.12 |
| Monocarboxylate porters | the hypothesis from Tr44175 | human MCT1 (P53985) and yeast Mch5 (Q08777), 2.A.1.13 |
| More distant MFS | root | Flr1 (P38124) and Tpo1 (Q07824), 2.A.1.2; Atr1 (P13090), 2.A.1.3; Dal5 (P15365), 2.A.1.14; *A. nidulans* CrnA (P22152), 2.A.1.8 |

The tree is rooted on the distant MFS proteins. Identity across MFS families is low, so a placement
counts only with strong support. As an independent test, the 335 proteins are scored with the MFS_1
HMM (PF07690) against PF00083 on Roihu, which needs no tree. Pho84, Git1 and Jen1 hit PF00083, so
their *S. cerevisiae* copies may already be family members; F11 uses the member instead of the UniProt
sequence if so, and that check also shows whether non-sugar porters are in the family.

## 3. Outgroups and rooting

If the ingroup is a supported clade, any sequence outside it roots it correctly, whatever its exact
sister. Rule per subfamily:

- **Outgroup sample:** 10–20 sequences from two or three neighbouring lineages in the whole-family
  trees (neighbouring, not necessarily sister, since the backbone is uncertain). Within each lineage,
  sequences spread over its diversity, preferring short terminal branches against long-branch
  attraction. Characterized proteins are preferred when a lineage has them.
- **One run, two outgroup lineages.** Each lineage must form a clade outside the ingroup, and rooting
  on either lineage must give the same ingroup root. This tests the root without a second run.
- **Checks on each subfamily tree:**
  1. the ingroup is monophyletic (if not: the boundary in `subfamilies.tsv` is revised and the tree
     rerun, with the decision recorded);
  2. each outgroup lineage is monophyletic and outside the ingroup;
  3. the outgroup root is compared with the midpoint and MAD roots (MAD: minimal ancestor deviation,
     with its ambiguity index), and any difference is explained. MAD assumes roughly clock-like
     evolution and tends to fall between paralog clades that evolve at different rates, so a
     disagreement is reported, not a reason to replace the outgroup root;
  4. the subfamily topology is compared with the same tips pruned from the whole-family trees, as a
     check, not a requirement.

## 4. Subfamily trees

For each subfamily, from full-length sequences:

1. **Sequences:** ingroup members and the orphan references that fall in the ingroup in the
   whole-family trees, plus the outgroup sample. Characterized references keep their tip IDs, so
   labels and iTOL files carry over.
2. **Alignment:** MAFFT L-INS-i up to about 2,500 sequences (all nine subfamilies); ClipKIT gappy 0.5
   as for the family, with the number of columns kept reported against the 455 of the whole family.
3. **Quick tree:** FastTree, for a first look at the ingroup and outgroup checks before the long run.
4. **Model and tree:** IQ-TREE 3 with ModelFinder restricted to the plausible matrices
   (`-mset Q.pfam,LG,WAG,JTT,Q.yeast`, as in the earlier aromatic trees, which saves most of the
   model-selection time), 1,000 UFBoot with NNI refinement and 1,000 SH-aLRT, `-safe`, fixed seed.
   Resources from the whole-family run: IQ-TREE chose 3–5 threads; 8 CPUs and 32 GB per job.
5. **Review:** the checks of §3, then iTOL files (taxonomy, anchors, functions) and a review note per
   subfamily, like `review_fasttree.md`.

Order: 1 (`hexose`, best anchored) first, as the test of the procedure; then 2–5; then 6–7; then 8. Tree 9 waits for the checks of §6.

## 5. What becomes generic

The aromatic families need the same steps (the class A flavin monooxygenases first), so they are
written for any family:

| Step | What | Kind |
|---|---|---|
| F10 | Choose subfamilies on the whole-family trees and record them in `families/<family>/subfamilies.tsv`: one row per set of tips, as an MRCA of named tips in a named tree, orthogroups, listed tips, or UniProt proteins, with the sample size for outgroup lineages | human decision, like `references.tsv` |
| F11 | `scripts/build_subfamilies.py`: each subfamily's tip list, outgroup sample and sequences, with a report (sizes, anchors, *D. squalens* proteins, orthogroups, tips in no subfamily) for review | script |
| F12 | Align and trim per subfamily | `align_family.sbatch` with a subfamily argument and L-INS-i |
| F13 | FastTree, then ModelFinder and IQ-TREE per subfamily | `family_tree.sbatch` with a subfamily argument and `-mset` |
| F14 | Ingroup and outgroup checks; outgroup, midpoint and MAD roots; iTOL files | script, then human review |

The whole-family tree is FastTree, and F7 draws it with its orthogroup, score and coverage rings for
choosing subfamilies. If a supported backbone is wanted, it comes from IQ-TREE on one
representative per MMseqs2 cluster (e.g. 50% identity, 1,000–2,000 sequences), as for OG0000018 in the
earlier aromatic analysis, rather than from all members.

## 6. Open points

- What to do with the tips in no subfamily (§2), and the scope of tree 9, which may widen to the whole
  lineage of 2,088 tips (0.97). OG0000544 (331 tips, 1.00), OG0000297 (405, 1.00), OG0000018 (450, 0.99),
  OG0000296 (310, 0.97) and OG0000015 (110, 0.79) meet the rule of §2 as FastTree clades, but none has a
  characterized member, and they may not be sugar porters at all. Two cheap checks on Roihu come first:
  the full run of F11, which shows whether Pho84, Git1 and Jen1 are family members, and the MFS_1 against
  PF00083 scores of §2, extended from the 335 proteins of tree 9 to all members.
