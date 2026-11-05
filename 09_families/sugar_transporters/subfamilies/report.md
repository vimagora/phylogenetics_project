# Subfamilies of sugar transporters

Written by `scripts/build_subfamilies.py` with `--lists-only`.

| Subfamily | Ingroup tips | Characterized | *D. squalens* | Outgroup tips | Main orthogroups |
|---|---|---|---|---|---|
| hexose | 690 | 77 | 4 | 16 | OG0000008 674, (none) 11, OG0002107 1 |
| mal | 903 | 12 | 0 | 16 | OG0000174 694, OG0000008 186, OG0016147 8 |
| cdt | 1080 | 20 | 7 | 16 | OG0000100 1066, OG0000102 3, OG0000395 2 |
| fsy | 462 | 14 | 2 | 16 | OG0000345 453, (none) 4, OG0091281 1 |
| itr | 390 | 7 | 2 | 10 | OG0001120 215, OG0003609 148, (none) 12 |
| galacturonate_quinate | 575 | 8 | 2 | 16 | OG0000008 572, (none) 1, OG0008922 1 |
| gut_fungi | 155 | 0 | 0 | 16 | OG0000008 151, OG0022899 4 |
| pentose | 2160 | 43 | 11 | 16 | OG0000008 2129, (none) 7, OG0092402 1 |
| weak | 335 | 1 | 3 | 20 | OG0000772 226, OG0000102 109 |

## hexose

| Set | Kind | Tree | Tips in the set | Support | Sample |
|---|---|---|---|---|---|
| ingroup | node 37 | fasttree | 690 | 1.000 | all |
| mal_group | node 28 | fasttree | 903 | 0.579 | 8 |
| galacturonate_quinate | node 44 | fasttree | 575 | 0.984 | 8 |

Orthogroups: OG0000008 674, (none) 11, OG0002107 1, OG0102364 1, OG0062669 1, OG0112303 1.

Characterized (72): AmMST1, BHT1, BHT2, BHT3, CgHxt1, CgHxt2, CgHxt4, CgHxt5, CgHxt6, Gal2, Ght1, Ght2, Ght3, Ght4, Ght5, Ght6, Glt1, Gxf1, Gxs1, HGT-1, HGT2, HXT1, HpGCR1, Hxs1, Hxt1, Hxt10, Hxt11, Hxt13, Hxt14, Hxt15, Hxt16, Hxt17, Hxt2, Hxt3, Hxt4, Hxt5, Hxt6, Hxt7, Hxt8, Hxt9, HxtB/MstC, HxtC, HxtD/MstA, HxtE, MST1.2, MST1.3, MoHXT1 (MoSnf3), MoHXT2, MoHXT3, MoST1 (MoRco3), MstA, MstC, MstE, MstG, MstH, Pshxt1, RAG1, RAG4, RCO-3, Rgt2, Snf3, Stp1, Str3, Sut1, Sut3, TmelHXT1, UmHxt1, XYT1, Xltr1p, YHT1, YHT2, YHT3.

Outgroup `mal_group`: Aspni_NRRL3_1-8653 (FsyX), Fusox3_1-2580761, Ustma-UmSrt1 (UmSrt1), Fusox3_1-2512555, Aspac1-53789, Kredeu1-490247, Fusox3_1-2030885, Sacce1-3598 (Mph3).

Outgroup `galacturonate_quinate`: Rhier1-738949, Rhimi59_2-1907711, Lepor2-622677, Trire_Chr-109101 (Gat1), Neucr2-10204 (qa-y), Arxad1-3619, Bcine1-11562 (GalA), Aspni_NRRL3_1-9845.

## mal

| Set | Kind | Tree | Tips in the set | Support | Sample |
|---|---|---|---|---|---|
| ingroup | node 28 | fasttree | 903 | 0.579 | all |
| hexose_group | node 37 | fasttree | 690 | 1.000 | 8 |
| galacturonate_quinate | node 44 | fasttree | 575 | 0.984 | 8 |

Orthogroups: OG0000174 694, OG0000008 186, OG0016147 8, OG0227384 1, OG0001222 1, OG0000527 1.

Characterized (13): Agt1/Mal11, FsyX, HpMAL2, Mal31, Mal61, MalP, Mbt1, Mlt1, Mph2, Mph3, Mrt, Tr69957, UmSrt1.

Outgroup `hexose_group`: Lipst1_1-27377, Fusox3_1-2474014, Trire_Chr-111552 (Stp1), Phaso-Pshxt1 (Pshxt1), Cryne_H99_1-3080, Plicr1-36565, Sacce1-1375 (Hxt7), Schpo1-4281 (Ght5).

Outgroup `galacturonate_quinate`: Rhier1-738949, Rhimi59_2-1907711, Lepor2-622677, Trire_Chr-109101 (Gat1), Neucr2-10204 (qa-y), Arxad1-3619, Bcine1-11562 (GalA), Aspni_NRRL3_1-9845.

## cdt

**Ingroup tips also in:** weak (3).

| Set | Kind | Tree | Tips in the set | Support | Sample |
|---|---|---|---|---|---|
| ingroup | node 22 | fasttree | 1080 | 0.888 | all |
| fsy_group | node 61 | fasttree | 462 | 0.875 | 8 |
| og0000544 | node 80 | fasttree | 331 | 1.000 | 8 |

Orthogroups: OG0000100 1066, OG0000102 3, OG0000395 2, (none) 2, OG0000079 1, OG0000527 1.

Characterized (21): CBT-1, CDT-1, CDT-2, CRT1 (T. reesei RUT-C30 strain), CdtC, CdtD, CdtG, CltA, CltB/LacpB, Crt1 (T. reesei QM6a), CtA, DhSyi1, LAC12, LacpA, Pc_ST, PiHXT5, ST16, ST3, ST6, ST7, StrC.

Outgroup `fsy_group`: Sacba_Sacpa-SbFsy1_SpaFsy1 (SbFsy1/SpaFsy1), Klula1-2859, Cenge3-614790, Cenge3-660017, PleosPC15_2-1093407, Debha1-1605 (DhSyl1), Capcor1-3620, Trire_Chr-115647 (Frt1).

Outgroup `og0000544`: Tulca1-163642, Stamag_NRRLY2024_1-2216, Rhimi59_2-1874264, Morel2-134289, Lepor2-231990, Dekbr2-8845, Calco1-476253, Aspor1-9098.

## fsy

| Set | Kind | Tree | Tips in the set | Support | Sample |
|---|---|---|---|---|---|
| ingroup | node 61 | fasttree | 462 | 0.875 | all |
| cdt_group | node 22 | fasttree | 1080 | 0.888 | 8 |
| og0000297 | node 64 | fasttree | 405 | 0.997 | 8 |

Orthogroups: OG0000345 453, (none) 4, OG0091281 1, OG0218944 1, OG0229254 1, OG0001033 1.

Characterized (13): CmFsy1, DhSgl1, DhSyl1, DhSyl2, FRT1, FST1, Frt1, Fsy1, Fsy1 (yeast-like), SbFsy1/SpaFsy1, SeFsy1, SuFsy1, ZrFsy1.

Outgroup `cdt_group`: Treme1-63529, Stano2-7930, Rhoto_IFO0880_4-9638, Trire_Chr-113691 (CRT1 (T. reesei RUT-C30 strain)/Crt1 (T. reesei QM6a)), Lipst1_1-120451, Fusox3_1-2583184, Debha1-5483 (DhSyi1), Aspnid-CltA (CltA).

Outgroup `og0000297`: Kredeu1-499307, Canpa1-3570, Morel2-132939, Lacind1-1943505, Fompi3-1120432, Amamu1-346138, Spipu1-8708, Pilcr1-69208.

## itr

| Set | Kind | Tree | Tips in the set | Support | Sample |
|---|---|---|---|---|---|
| ingroup | node 68 | fasttree | 390 | 0.840 | all |
| ecoli | tips |  | 2 |  | all |
| hexose_group | node 37 | fasttree | 690 | 1.000 | 8 |

Orthogroups: OG0001120 215, OG0003609 148, (none) 12, OG0002200 4, OG0003502 2, OG0156484 1.

Characterized (6): GpMst1, ITR1, ITR2, XP_012052055.1, XP_012053726.1, itr2.

Outgroup `ecoli`: Escco-GalP, Escco-AraE.

Outgroup `hexose_group`: Lipst1_1-27377, Fusox3_1-2474014, Trire_Chr-111552 (Stp1), Phaso-Pshxt1 (Pshxt1), Cryne_H99_1-3080, Plicr1-36565, Sacce1-1375 (Hxt7), Schpo1-4281 (Ght5).

## galacturonate_quinate

| Set | Kind | Tree | Tips in the set | Support | Sample |
|---|---|---|---|---|---|
| ingroup | node 44 | fasttree | 575 | 0.984 | all |
| hexose_group | node 37 | fasttree | 690 | 1.000 | 8 |
| pentose_group | node 10 | fasttree | 2160 | 0.920 | 8 |

Orthogroups: OG0000008 572, (none) 1, OG0008922 1, OG0064474 1.

Characterized (8): AHA12287.1, GalA, Gat-1, Gat1, Gat2, GatA, Qut-D, qa-y.

Outgroup `hexose_group`: Lipst1_1-27377, Fusox3_1-2474014, Trire_Chr-111552 (Stp1), Phaso-Pshxt1 (Pshxt1), Cryne_H99_1-3080, Plicr1-36565, Sacce1-1375 (Hxt7), Schpo1-4281 (Ght5).

Outgroup `pentose_group`: Pilcr1-94396, Dales1-253029, Pilcr1-78995, Aspni_NRRL3_1-4757 (LatA), Aspni_NRRL3_1-2351 (XltD), Triat2-265450, Trire_Chr-108806 (Str2), Trire_Chr-115600 (Trhxt1).

## gut_fungi

| Set | Kind | Tree | Tips in the set | Support | Sample |
|---|---|---|---|---|---|
| ingroup | node 150 | fasttree | 155 | 0.997 | all |
| galacturonate_quinate | node 44 | fasttree | 575 | 0.984 | 8 |
| pentose_group | node 10 | fasttree | 2160 | 0.920 | 8 |

Orthogroups: OG0000008 151, OG0022899 4.

Characterized (0): none.

Outgroup `galacturonate_quinate`: Rhier1-738949, Rhimi59_2-1907711, Lepor2-622677, Trire_Chr-109101 (Gat1), Neucr2-10204 (qa-y), Arxad1-3619, Bcine1-11562 (GalA), Aspni_NRRL3_1-9845.

Outgroup `pentose_group`: Pilcr1-94396, Dales1-253029, Pilcr1-78995, Aspni_NRRL3_1-4757 (LatA), Aspni_NRRL3_1-2351 (XltD), Triat2-265450, Trire_Chr-108806 (Str2), Trire_Chr-115600 (Trhxt1).

## pentose

| Set | Kind | Tree | Tips in the set | Support | Sample |
|---|---|---|---|---|---|
| ingroup | node 10 | fasttree | 2160 | 0.920 | all |
| hexose_group | node 37 | fasttree | 690 | 1.000 | 8 |
| galacturonate_quinate | node 44 | fasttree | 575 | 0.984 | 8 |

Orthogroups: OG0000008 2129, (none) 7, OG0092402 1, OG0156532 1, OG0123645 1, OG0000006 1.

Characterized (41): An25, AraT, CgMfs1, DhStl1, GT1, GlpA, GlpB, GlpC, Gtt1, HGT1, HXT1, HxtA, LAT-1, LAT1, LAT2, LatA, MST3.1, MtLAT-1, PcAraT, RTO4_9841, Stl1, Stl1 (Hgt10), Str1, Str2, Tmel2281, TrXlt1/Trire_104072, Trhxt1, WaStl1, XAT1, XltA, XltB, XltC, XltD, XtrD, Xut1, Xut3, Xylh, Xyp29, YHT4, YHT5, YHT6.

Outgroup `hexose_group`: Lipst1_1-27377, Fusox3_1-2474014, Trire_Chr-111552 (Stp1), Phaso-Pshxt1 (Pshxt1), Cryne_H99_1-3080, Plicr1-36565, Sacce1-1375 (Hxt7), Schpo1-4281 (Ght5).

Outgroup `galacturonate_quinate`: Rhier1-738949, Rhimi59_2-1907711, Lepor2-622677, Trire_Chr-109101 (Gat1), Neucr2-10204 (qa-y), Arxad1-3619, Bcine1-11562 (GalA), Aspni_NRRL3_1-9845.

## weak

**Ingroup tips also in:** cdt (3).

| Set | Kind | Tree | Tips in the set | Support | Sample |
|---|---|---|---|---|---|
| ingroup | orthogroups |  | 335 |  | all |
| sugar_porters | tips |  | 10 |  | all |

Orthogroups: OG0000772 226, OG0000102 109.

Characterized (1): Tr44175.

Outgroup `sugar_porters`: Sacce1-3057 (Hxt1), Aspni_NRRL3_1-3147 (MstA), Aspni_NRRL3_1-11715 (XltA), Picst3-33008 (Xut1), Sacce1-759 (Mph2), Aspor1-9642 (MalP), Neucr2-1907 (CDT-1), Klula1-1118 (LAC12), Sacce1-1538 (ITR1), Schpo1-573 (itr2).

- `other_pf00083`: Sacce-Pho84 (UniProt P25297): not checked (--lists-only)
- `other_pf00083`: Sacce-Git1 (UniProt P25346): not checked (--lists-only)
- `other_pf00083`: Sacce-Jen1 (UniProt P36035): not checked (--lists-only)
- `monocarboxylate`: Homsa-MCT1 (UniProt P53985): not checked (--lists-only)
- `monocarboxylate`: Sacce-Mch5 (UniProt Q08777): not checked (--lists-only)
- `distant_mfs`: Sacce-Flr1 (UniProt P38124): not checked (--lists-only)
- `distant_mfs`: Sacce-Tpo1 (UniProt Q07824): not checked (--lists-only)
- `distant_mfs`: Sacce-Atr1 (UniProt P13090): not checked (--lists-only)
- `distant_mfs`: Sacce-Dal5 (UniProt P15365): not checked (--lists-only)
- `distant_mfs`: Aspnid-CrnA (UniProt P22152): not checked (--lists-only)

## Tips in no subfamily

1934 of 8681 tips are in no ingroup, 12 of them *D. squalens*.

Orthogroups: OG0000018 476, OG0000297 464, OG0000544 327, OG0000296 310, OG0000015 115, OG0000079 45, OG0000118 37, OG0000008 19, OG0012023 15, OG0005150 13, OG0000065 11, OG0004810 9.

Characterized: none.
