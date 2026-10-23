# Review of the FastTree tree

Review 1 of `expectations.md`, on `09_families/sugar_transporters/tree/fasttree.nwk` (FastTree, LG +
Gamma; 8,681 tips: 8,640 members, 27 orphans, 14 outgroups), 2026-10-23. Support values are FastTree's
local supports. The tree was rooted on *E. coli* GalP + AraE, as `expectations.md` provides when the
three bacterial proteins do not form a clade.

## Expectations

| check | result |
|---|---|
| 1. *A. thaliana* STPs form a clade | met (support 1.0) |
| 2. Human class I GLUTs form a clade | met (1.0) |
| 3. GalP and AraE are sisters | met (1.0) |
| 4. STPs and GLUTs fall outside all cohort sugar porters, or next to one subfamily | not met: both fall within a 390-tip clade with the myo-inositol transporters (OG0001120) and OG0003609; STPs are sister to OG0003609 (0.55), GLUTs sit next to eight ITRs and *Geosiphon* GpMst1 (0.86) |
| 5. (soft) the three *E. coli* proteins form a clade | not met: XylE is sister to four proteins of anaerobic gut fungi (`Caecom1-444819`, `Anasp1-296144`, `Neosp1-707593`, `Pirfi3-586812`), support 1.0 |
| 6. (soft) STPs and GLUTs join each other before the bacteria | not met (see 4) |
| 7. (soft) GLUT1 and GLUT3 are sisters | met (0.79) |
| 8. MAL (OG0000174) forms a clade | nearly: 688 of 700 in one clade; 12 members from 11 genomes elsewhere |
| 9. CDT/lactose (OG0000100) forms a clade | met: all 1,066 in a 1,080-tip clade (0.89) |
| 10. Fsy (OG0000345) forms a clade | met: all 453 in a 462-tip clade (0.88), with the four *Saccharomyces* Fsy orphans |
| 11. ITR (OG0001120) forms a clade | nearly: 200 of 217 in one clade; the rest near OG0003609 |
| 12. (soft) subfamilies follow the species tree | not assessed |
| 13. (soft) sensors without long separate branches | met: terminal branches 0.10–0.22 (median of all tips 0.12); yeast Snf3, Rgt2, RAG4 and Hxs1 group together, as do RCO-3 and MoST1 |

## Decisions

- **Root.** The final tree is rooted on GalP + AraE only (`root` = TRUE in `outgroups.tsv`). XylE,
  the STPs and the GLUTs stay in the tree as landmarks: they nest within the fungal sugar porters, so
  they cannot root them. XylE's grouping with gut-fungal proteins may reflect horizontal gene
  transfer from bacteria [EK].
- **Alternative root, to be reported.** Midpoint rooting places the root between OG0000772 +
  OG0000102 (342 proteins) and all other members. These are the weakest PF00083 members (mean score
  127 and 38, against 264 for all members); OG0000102 may belong to another MFS family. Both are kept.
- **Next.** ModelFinder and IQ-TREE on the full alignment. Clade-level trees, with sister clades as
  local outgroups, are to be defined from the IQ-TREE tree.
