# Expectations for the sugar transporter tree

Written on 2026-10-23, before any tree of the family was built. They are checked by eye at the two
review points (FastTree tree, then the IQ-TREE tree) and are not tested by any script. A failed hard
check means the alignment, trimming or tree building should be examined before the tree is
interpreted; a failed soft check is a result to look into, not by itself a sign of an artifact.

Outgroups are listed in `outgroups.tsv` (14 sequences). Pairwise identities below are from global
alignments of the outgroup sequences (BLOSUM62), and only describe the sequences; the sources for
each expectation are given with it ([EK]: from expert knowledge, to be added to Zotero).

## Outgroup

**Hard checks**
1. The seven *A. thaliana* STPs form one clade. STPs are a plant-specific family [EK: Büttner 2010]
   and were the outgroup of Xu et al. [@xu2024Genomewidepredictiontranscriptome], where they fell
   outside the fungal sugar transporters. Identities: 47–65% among STPs, 27–37% to any other outgroup.
2. The four human class I GLUTs (GLUT1–4) form one clade [EK: class I of the SLC2 family].
   Identities: 53–66% among them, 27–35% to any other outgroup.
3. *E. coli* GalP and AraE are sister to each other. Identity 65%, against 29–37% to anything else.
4. The STPs and the class I GLUTs fall outside all cohort sugar porters, or at most next to a single
   fungal subfamily, as in Xu et al. for the STPs.

**Soft checks**
5. The three *E. coli* proteins form one clade. XylE is about 35% identical to every outgroup,
   bacterial or eukaryotic, so it may branch elsewhere.
6. The eukaryotic outgroups (STPs, GLUTs) join each other before either joins the bacteria.
7. Within class I, GLUT1 and GLUT3 are sisters [EK]. Identities alone do not support this (GLUT1 is
   65% identical to GLUT3 and 66% to GLUT4).

The tree is drawn rooted on the bacterial proteins (`root` = TRUE in `outgroups.tsv`); if check 5
fails, the root is placed on GalP + AraE.

## Ingroup

**Hard checks**: each characterized subfamily forms a clade, as its orthogroup suggests:
8. maltose and α-glucoside transporters (MAL, MalP, Mrt): OG0000174;
9. cellodextrin and lactose transporters (CDT, CltB, Lac12, DsST1–6): OG0000100;
10. Fsy fructose symporters (Fsy1, Frt1, DsST8–9): OG0000345;
11. myo-inositol transporters (ITR1, ITR2, Itr2, DsST7): OG0001120.

**Soft checks**
12. Within each subfamily, members follow the species tree at the level of classes, apart from
    lineage-specific expansions (for example the Hxt expansion of *Saccharomyces*).
13. The sensors (Snf3, Rgt2, RAG4, Hxs1, RCO-3, MoST1, MoHXT1) sit within the hexose transporters of their own
    lineages, not on long separate branches.
