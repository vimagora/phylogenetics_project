---
title: Cohort selection
date: 2026-10-19
tags:
  - cohort
  - lignocellulose
  - mycocosm
  - methods
bibliography: C:/Users/vimag/Documents/obsidian-vaults/Fungal_research/references.bib
csl: C:/Users/vimag/Documents/obsidian-vaults/Fungal_research/90 Templates/csl/apa.csl
---

# Cohort selection

## 1. Methods

### 1.1 Candidate genomes and sampling design

Genomes were drawn from the JGI MycoCosm catalogue of published, non-superseded fungal genomes (2,802 portals, retrieved 13 October 2026). The cohort serves the phylogenetic characterization of four families: sugar transporters, transcription factors, aromatic-catabolic enzymes and CAZymes. Function and substrate specificity are inferred from placement relative to characterized proteins.

A genome was considered only if it contains characterized protein (an anchor); it is a close relative needed to resolve orthology; it is situated around an anchor or a lineage-specific expansion; it represents one side of a lifestyle or decay-mode contrast; it is part of a lineage that would otherwise leave a long, unbroken branch. The coverage was checked across taxonomy, including outgroups for rooting; decay mode: white, brown and soft rot, and intermediate or uncertain decay; and lifestyle: wood, litter, soil and dung saprotrophs, ecto-, ericoid, orchid and arbuscular mycorrhizal fungi, endophytes, pathogens and obligate biotrophs, yeasts, and anaerobic gut fungi. Sister relationships and distances were taken from the kingdom-wide phylogeny of Li et al. [@li2021genomescalephylogenykingdom], using the published maximum-likelihood tree (1,672 taxa, 290 genes). In all, 186 candidate portals were downloaded and assessed. They comprise the cohort candidates. 

### 1.2 Data retrieval and assembly statistics

For each portal, the gene models (GFF3), the repeat-masked assembly, and the default protein and CDS sets were downloaded through the MycoCosm API. Comment lines at the start of the FASTA files were removed, and "/" in protein headers was replaced by an underscore, because BUSCO does not accept it.

Assembly statistics (sequences, total length, N50, GC content and ambiguous bases) were computed with SeqKit 2.13.0. The assemblies are soft-masked. Across the cohort and its shadows, ambiguous bases had a median of 0.24%. They exceeded 10% in nine assemblies, the highest being *Sphaerobolus stellatus* (30.7%) and *Conidiobolus coronatus* (20.5%). Cohort scaffold N50 had a median of 1.5 Mb (minimum 31 kb, *Caecomyces churrovis*).

### 1.3 Collapse of alternative and redundant gene models

MycoCosm GFF3 files store each alternative transcript, and some redundant models from different gene callers, as a separate gene. Each proteome was therefore reduced to one protein per locus. 
Models on the same strand were treated as one locus when their coding sequences shared at least 50% of the shorter model's CDS length, directly or through a chain of overlaps. The longest protein was kept, with ties broken by gene name and then by the lowest protein ID. Models without a protein (≤ 2% per genome) contained internal stop codons where checked. They were kept aside as candidate pseudogenes.

Collapsing removed 4,794 models from 69 of the 136 proteomes with recorded counts; counts for the remaining 50 are to be added. The largest reductions were in *Neurospora crassa* (1,055 models), *Puccinia triticina* (807), *P. graminis* (773), *Phlebia radiata* (549), and *Spizellomyces punctatus* (472).

### 1.4 Completeness

Completeness was assessed with BUSCO 6.1.0 in protein mode against fungi_odb12 (1,122 markers), on the original and the collapsed proteomes. Collapsing turned 583 duplicated markers into single-copy ones. The largest drops were in *Fusarium oxysporum* 4287 v2 (22.9% → 0.7%), *N. crassa* (8.7% → 0.1%) and *P. triticina* (8.6% → 3.1%). Three complete markers were lost, each because the longest model kept was a fragment or a probable read-through fusion. No genome's completeness changed by more than 0.1 percentage points.

### 1.5 Quality tiers

Since several lineages genuinely lack markers that are conserved in Dikarya, completeness was scored relative to each genome's lineage. In the lineage-relative completeness (LRC) score, a marker missing or fragmented in a genome counts against it only if more than a third of the genome's comparison group have it. That group is the genome's order, or its class or nearest relative when the order has fewer than three genomes. The cuts were set from the score distribution of the candidates: median 1,091, and a robust SD (1.5 × MAD) of 23.7.

| tier | rule | allowed use |
|---|---|---|
| **C**, count-grade | LRC ≥ 1,067.3 (median − 1 SD) | counts, expansions and losses |
| **C+**, count-grade with confirmation | 1,019.8 ≤ LRC < 1,067.3 | counts; every inferred absence is confirmed against the assembly |
| **P**, presence-only | LRC < 1,019.8 (median − 3 SD), or more than 5% fragmented markers | presences; absences coded as unknown |
| **L**, losses-only | more than 6% duplicated markers after step A | losses; no counts |

Genomes whose duplication reflects real genes are exempt from the L rule. The genomes with a duplication flag instead are *Saccharomyces cerevisiae* and *Rhizopus microsporus*, which have whole-genome duplications, *Linnemannia elongata*, where two independently annotated strains show the same duplication, the genus *Neocallimastix* where all three sequenced species show 36–41% duplicated markers, against 3.8–5.3% in other anaerobic gut fungi.

### 1.6 Shadow proteomes

Sixteen second strains, haplotypes, older annotations or congeners were kept as shadow proteomes. Each is searched with the orthogroup HMMs for presence only and never enters OrthoFinder, and species-level presence is coded as the union of the cohort genome and its shadow. They carry the protein IDs used in the literature, which are mapped onto the cohort genome by exact or reciprocal-best match: *F. oxysporum* v2, *Dichomitus squalens* LYAD-421, *Coprinopsis cinerea* Okayama 7 and *Cryptococcus deneoformans* JEC21. They separate alleles from paralogs in rust dikaryons. They confirm real duplications, in *Linnemannia* and *Neocallimastix*. They confirm losses, or give a second annotation of an anchor species.

## 2. Decisions by taxonomic group and lifestyle HERE

Higher ranks follow Li et al. [@li2021genomescalephylogenykingdom, journal p. 1653]. Chytridiomyceta holds Chytridiomycota, Monoblepharidomycota and Neocallimastigomycota, and Mucoromycota holds three subphyla, including Glomeromycotina. The tables keep MycoCosm's phylum label "Glomeromycota" for *Rhizophagus*. Saccharomycotina follows the class and order system of Groza et al. 2024 [EK] (not in vault; add to Zotero).
- Where a MycoCosm name is outdated, the current name is given in §3.

**Overall composition.**
- The cohort spans 9 phyla, 28 classes and 58 orders: 72 Basidiomycota, 67 Ascomycota and 11 early-diverging genomes.
- Fifteen genomes serve all four families as anchors or key contrasts: `Aspni_NRRL3_1`, `Aspnid1`, `Bcine1`, `Dicsqu464_2`, `Fusgr1`, `Lacbi2`, `Neucr2`, `Phchr2`, `Pirfi3`, `Rhimi59_2`, `Sacce1`, `Spipu1`, `Spoth2`, `Trire_Chr` and `Ustma2_2`.
- Abbreviations: ST, sugar transporter; TF, transcription factor.

### 2.1 Taxonomic groups

**Early-diverging lineages (11 genomes).**
- **Root.** The species tree is rooted on the zoosporic lineages.
  - *Spizellomyces punctatus* (`Spipu1`, Chytridiomycota) is the root for every family tree.
  - *Hyaloraphidium curvatum* (`Hyacur1`, Monoblepharidomycota, C+) adds a second zoosporic lineage. It is sister to the anaerobic gut fungi in concatenation (BS 85%) [@li2021genomescalephylogenykingdom, journal p. 1659].
- **Blastocladiomycota and Zoopagomycota** break the long branches to Dikarya.
  - Blastocladiomycota is represented by *Catenaria anguillulae* (`Catan2`, P), chosen over *Allomyces* (§2.3).
  - Zoopagomycota is represented by *Conidiobolus coronatus* (`Conco1`, P).
  - Both are presence-only and are each the only processed genome of their phylum.
  - The finding that Zoopagomycetes "only harbor STs from clade A and clade E" [@xu2025Expansionmaltosesucrose, p. 4] therefore cannot be retested with `Conco1`.
- **Mucoromycota.** *Rhizopus microsporus* (`Rhimi59_2`) and *Linnemannia elongata* (`Morel2`) polarize gains in Dikarya; both carry duplication flags. Mortierellomycotina is "sister to Mucoromycotina and Glomeromycotina (BS = 98%)" [@li2021genomescalephylogenykingdom, journal p. 1660].
  - *Rhizopus* sits at the origin of ST clades D and H.
  - With *Ustilago*, it is needed to root the fungal duplication into HmgA (clade I, tyrosine and phenylacetate catabolism) and MhdA (clade II, G-unit lignin monomers) [@li2025Fungalbioconversionligninderived, PDF p. 5; @lubbers2021Vanillicacidmethoxyhydroquinone, p. 7].
- **Arbuscular mycorrhiza.** *Rhizophagus irregularis* (`RhiirDAOM197198_2`) is the AM contrast. It has no GH6, GH7, PL1, PL3, PL4, PL9 or CBM1 [@tisserant2013Genomearbuscularmycorrhizal, PDF p. 4], and "no secreted invertase or sucrose transporter was identified" [@tisserant2013Genomearbuscularmycorrhizal, PDF p. 5].
- **Anaerobic gut fungi.** Four genera are sampled; see §2.2.

**Saccharomycotina and Taphrinomycotina (17 genomes).** Half of the characterized sugar transporters (87 of 171) are Saccharomycotina proteins, and ascomycete genomes average about 90 STs against 26 in basidiomycetes [@xu2025Expansionmaltosesucrose, p. 3]. The yeasts were therefore sampled as a backbone.
- ***Lipomyces starkeyi*** (`Lipst1_1`, C+) roots all yeast transporter, sensor and Mig1/Gal4 subfamilies, as the earliest Saccharomycotina lineage.
- **Dipodascales:**
  - *Yarrowia lipolytica* (`Yarli1`) carries 6 YHT anchors and the TCM1 cinnamate-hydroxylase-like anchor.
  - *Blastobotrys adeninivorans* (`Arxad1`) carries the Agdc1p gallate decarboxylase [@li2025Fungalbioconversionligninderived, pp. 12, 14].
- **Serinales** (5 genomes):
  - *Scheffersomyces stipitis* (`Picst3`), the xylose fermenter;
  - *Candida intermedia* (`Canint_NRRLY981_1`), the Gxf1/Gxs1 anchors;
  - *Debaryomyces hansenii* (`Debha1`), six polyol and inositol porters;
  - *C. parapsilosis* and *C. albicans*, both kept because they are not equivalent. *C. parapsilosis* (`Canpa1`) is the organism of the gentisate/hydroxyquinol cluster (Mnx1/2/3, Hdx1, Gdx1, Fph1) [@li2025Fungalbioconversionligninderived, pp. 5, 18]. *C. albicans* (`Canalb1`) lacks it, giving a sharp presence/absence contrast, and carries Hqd2 and ITR1.
- **Pichiales** (3 genomes):
  - *Ambrosiozyma monospora* (C+, the LAT1/LAT2 L-arabinose anchors) is the only count-grade genome.
  - *Ogataea polymorpha* (`Hanpo2`, P) carries HpGcr1 and HpMAL2.
  - *Brettanomyces bruxellensis* (`Dekbr2`, P) carries BHT1–3, DbPAD2 and vinylphenol reductase [@li2025Fungalbioconversionligninderived, PDF p. 7].
  - Pichiales transporter counts therefore rest on `Ambmon_NRRLY1484_1` alone.
- **Saccharomycetales:**
  - *S. cerevisiae* (`Sacce1`), the post-whole-genome-duplication anchor set (27 STs, Snf3/Rgt2 sensors, Fdc1/Pad1);
  - *K. lactis* (`Klula1`), chosen over *K. marxianus* as the pre-duplication relative;
  - *Zygosaccharomyces rouxii* (`Zygrou1`), the second pre-duplication genome and the ZrFsy1 anchor.
  
  *S. cerevisiae* "cannot metabolize pentose sugars and lacks the corresponding STs in its genome" [@xu2024Genomewidepredictiontranscriptome, p. 2].
- ***Starmerella magnoliae*** (Trigonopsidales, CmFsy1) is used as P (§1.5).
- **Taphrinomycotina:** *Schizosaccharomyces pombe* (`Schpo1`; 7 Ght/itr anchors) roots the Ascomycota, and *Saitoella complicata* (`Saico1`) shortens its branch.

**Pezizomycotina, early classes.**
- **Orbiliomycetes:** *Arthrobotrys oligospora* (`Artol1`) brackets the root of Pezizomycotina.
- **Pezizomycetes:**
  - *Tuber melanosporum* (`Tubme1v2`) is the ECM origin, with a restricted PCWDE set [@murat2018Pezizomycetesgenomesreveal, PDF p. 2, 4].
  - *Morchella importuna* (`Morco1`) is its saprotrophic comparator, and dates the origin of XlnR, CLR and GaaR.
  - *Ascobolus immersus* (`Ascim1`, C+) is the dung saprotroph and the exact strain of the CE15 AiGE1.

**Eurotiomycetes (12 genomes).**
- **Eurotiales** is sampled densely, because it holds most characterized regulators and enzymes. The anchor sets are:
  - *A. niger* NRRL3: 11 ST, 10 TF and more than 25 aromatic-enzyme anchors;
  - *A. nidulans*: 14 ST and 11 TF anchors;
  - *A. oryzae*: XlnR, AmyR, MalR, CreA and ManR;
  - *A. aculeatus*: ClbR and the only characterized AA16;
  - *A. fumigatus* (`Aspfu1`): the HppD–HmgA cluster and AfAA9/AA11;
  - *Penicillium oxalicum*: ClrB, CxrA and the CdtC/D/G transporters;
  - *P. chrysogenum* (`Pench1`): an anchor organism for all four families;
  - *Penicillium occitanis* (`Penoc1`): the closest genome to *Talaromyces cellulolyticus* (TacA, TctA);
  - *Thermoascus aurantiacus*: the thermophile.
- **Why the density.** Duplications of XlnR/AraR, AmyR/InuR/Ace3 and ClbR need it. XlnR orthologs are "present in almost all filamentous Ascomycota", but their target set "is highly species specific" [@benocci2017Regulatorsplantbiomass, p. 9].
  - *A. nidulans* cannot grow on cinnamic acid because its *padA* is truncated [@lubbers2025UpdatedPerspectiveAromatic, pp. 4–5, 7].
  - ManR "is no longer present in the latest version of the A. oryzae genome" [@benocci2017Regulatorsplantbiomass, p. 10], so absences in `Aspor1` (C+) are checked against the assembly.
- **Outside Eurotiales:**
  - *Coccidioides immitis* (`Cocim1`, Onygenales) is a control with a published CLR loss: *Coccidioides* "have lost their clr-1 homolog" [@coradetti2012Conservedessentialtranscription, p. 7400].
  - Two Chaetothyriales give the candidate VibMO1 donor clade: *Cladophialophora psammophila* and *Capronia coronata* (`Capcor1`, the CcFDC UbiD anchor [@li2025Fungalbioconversionligninderived, p. 9]).

**Dothideomycetes (8 genomes).**
- **Pathogens:**
  - the latent apoplastic pathogen *Zymoseptoria tritici* and its biotrophic sister *Cladosporium fulvum*;
  - the necrotroph *Cochliobolus sativus* (exact strain of BsGH7-3; proxy for the *Alternaria* CreA);
  - *Parastagonospora nodorum* (`Stano2`, P; SnGE1);
  - *Macrophomina phaseolina* (root-mycobiome isolate; Cel45A).
- **Others:**
  - the yeast-like saprotroph *Aureobasidium pullulans*;
  - the Gloniaceae pair, the ECM *Cenococcum geophilum* against the saprotroph *Glonium stellatum* (§2.2).

**Leotiomycetes (6 genomes).**
- **Necrotrophs:** *Botrytis cinerea* (`Bcine1`) is the C-tier Sclerotiniaceae genome. *Sclerotinia sclerotiorum* (`Sclsc1`, P) is kept for its CRE1 and SsShy1 anchors. Botrytis and Sclerotinia are "especially suited for pectin decomposition" (44 and 33 putative pectin-specific CAZymes) [@amselem2011GenomicAnalysisNecrotrophic, PDF p. 11–12].
- **Obligate biotroph:** *Blumeria hordei* (`Blugr2`).
- **Ericoid mycorrhizae:** *Oidiodendron maius* and *Rhizoscyphus ericae*.
- **Root endophyte:** *Leptodontidium* sp. (`Lepor2`), which replaces *Phialocephala*. It has about 164 STs, the largest repertoire known, and comes from the root-mycobiome study in which CAZyme repertoires were analysed by lifestyle [@mesny2021Geneticdeterminantsendophytism, PDF p. 5, 7].

**Sordariomycetes (20 genomes).**
- **Hypocreales:**
  - *Trichoderma reesei* (`Trire_Chr`): the core cellulase and Ace1/2/3–XYR1 anchors, and a reduced ST repertoire;
  - *T. atroviride*: the mycoparasite contrast;
  - *Fusarium graminearum* (`Fusgr1`);
  - *F. oxysporum* 4287 (`Fusox3_1`, with `Fusox2` as the anchor-ID shadow);
  - the insect pathogens *Metarhizium robertsii* (also a root endophyte; Mrt) and *Akanthomyces lecanii* (`Corco1`; proxy for CYP1 and DIO4).
  - The third *Trichoderma* genome was cut; the XYR1 anchor of *T. harzianum* is placed through `Trire_Chr`.
- **Sordariales:**
  - *Neurospora crassa* (13 ST and 9 TF anchors, seven AA9);
  - the thermophiles *Thermothelomyces thermophilus* (`Spoth2`, the largest AA9 set; MtLAT-1; MtCLR-2) and *Thermothielavioides terrestris*;
  - the dung saprotroph *Podospora anserina*;
  - *Chaetomium globosum* (`Chagl_1`, P), the classic soft-rot organism.
  - GH61 (AA9) is expanded in the Chaetomiaceae [@berka2011Comparativegenomicanalysis, PDF p. 2].
- **Order gaps filled:**
  - Glomerellales: *Colletotrichum graminicola* (hemibiotroph; 4 Hxt, Mbt1, four AA5_2 anchors) and *Verticillium dahliae* (C+; VdPelB);
  - Diaporthales: *Cytospora (Valsa) mali* (`Valma1`, VmHbh1);
  - Ophiostomatales: *Grosmannia clavigera* (`Grocl1`), a CLR control for which "a clr-2 homolog was not identified" [@coradetti2012Conservedessentialtranscription, p. 7400];
  - Magnaporthales: *Pyricularia oryzae* (`Magor1`; CreA, ARA1, XLR1).
- **Xylariales:** soft rot and endophytes, with four genomes (§2.2).

**Pucciniomycotina and Ustilaginomycotina.**
- **Rusts:**
  - *Melampsora larici-populina* (haplotype 0, with haplotype 1 as shadow);
  - *Puccinia graminis* (phased haplotype A, C+, with v2 as shadow).
  
  Rusts have a relatively small GH set, similar to that of *Laccaria* [@duplessis2011Obligatebiotrophyfeatures, PDF p. 5].
- **Yeasts:**
  - *Rhodotorula toruloides* (`Rhoto_IFO0880_4`) carries the RTO4_13270 galacturonate regulator. That Gal4-like protein is "not at all conserved beyond the Pucciniomycotina" [@protzko2019GenomewideEnzymaticAnalysis, p. 9]. The RTO4_9841 transporter "was also both transcriptionally induced and required for robust growth on D-galUA" [@protzko2019GenomewideEnzymaticAnalysis, PDF p. 5].
  - *Cystobasidium minutum* (`Rhomi1`) fills Cystobasidiomycetes and carries the P450rm (CYP53) anchor.
- **Ustilaginomycotina:** *Ustilago maydis* (`Ustma2_2`) is the biotroph with few secreted cell-wall-degrading enzymes [@kamper2006Insightsgenomebiotrophic, PDF p. 3]. It also carries the Shy1, Srt1 and Hxt1 anchors.
- **The base of Basidiomycota** cannot reject a hard polytomy [@li2021genomescalephylogenykingdom, journal p. 1660]. The relationships among the three subphyla are therefore treated as unresolved.

**Tremellomycetes and Dacrymycetes.**
- ***Cryptococcus*:** the cohort genome is *C. neoformans* H99 (`Cryne_H99_1`), not *C. deneoformans* JEC21 (a shadow), because the listed transporter accessions look like H99 RefSeq IDs [EK] and H99 is C.
  - It carries an ITR expansion of 11 genes, against 3 in its sister *Tremella mesenterica* (`Treme1`).
  - *Cutaneotrichosporon oleaginosum* (`Triol1`) is the nearest genome to *Trichosporon cutaneum* (PhhY, TcMLE, salicylate hydroxylase).
- **Dacrymycetes:** *Calocera cornea* holds the only large basidiomycete ST expansion known, 66 genes [@xu2025Expansionmaltosesucrose, p. 3–4]. It is bracketed by *Dacryopinax primogenitus*. Both are brown rots without PODs.

**Agaricomycetes (62 genomes, 16 orders).**
- **Early orders: the frame for the origin of class II peroxidases.** PODs are reconstructed in the ancestor of Auriculariales and the remaining Agaricomycetes, and the earliest sampled clades (Cantharellales, Sebacinales, Trechisporales, Dacrymycetes) lack them [@nagy2016ComparativeGenomicsEarlyDiverging, PDF p. 1]. The cohort samples:
  - **Cantharellales:** *Botryobasidium botryosum*, with 32 LPMOs and 3 CDH and no POD [@riley2014Extensivesamplingbasidiomycete, PDF p. 2]; the ECM *Hydnum rufescens*; the orchid mycorrhizal *Tulasnella calospora* (P).
  - **Sebacinales:** the endophyte *Serendipita indica* (`Pirin1`).
  - **Trechisporales:** *Sistotremastrum suecicum* (`Sissu1`).
  - **Phallomycetidae:** *Sphaerobolus stellatus* (`Sphst1`, P and L). Its oxidoreductase expansions track a random-gene control [@nagy2016ComparativeGenomicsEarlyDiverging, PDF p. 4], so its families are not counted.
  - **Auriculariales:** *Auricularia subglabra* and *Exidia glandulosa* (L).
  - **Corticiales:** *Punctularia strigosozonata* (C+; PsCE CE15 and PsAAOA).
  - **Atheliales:** *Piloderma olivaceum* (P and L; ECM).
  - **Amylocorticiales:** *Plicaturopsis crispa* (P and L).
- **Polyporales (18 genomes): the densest order.** It holds most white-rot and many aromatic anchors. The AA1, AA2 and AA14 paralog clades need dense sampling to resolve around their anchors [@hage2021Genefamilyexpansions, PDF p. 6–7]. The 14 white rots are:
  - phlebioid: *Phanerochaete chrysosporium*, *Phlebia radiata* and *Phlebiopsis gigantea*, plus *Phanerochaete carnosa* (P);
  - selective delignifiers: *Gelatoporia subvermispora* (P) and its C-tier sister *Obba rivulosa*;
  - *Dichomitus squalens* (the study's model organism) and *Trametes versicolor*;
  - *Trametes cinnabarina*, *T. coccinea*, *Fomes fomentarius* and *Ganoderma lucidum*, which date the FMO expansions;
  - *Panus rudis* (PanI, and the proxy for the *Lentinus similis* AA9);
  - *Bjerkandera adusta* (`Bjead1_1`; the only Polyporales VP anchors, and BaAAO).
  
  Four antrodia-clade brown rots complete the order (§2.2).
- **Hymenochaetales, Russulales, Gloeophyllales and Jaapiales.**
  - Hymenochaetales: two white rots, *Fomitiporia mediterranea* (C+) and *Xylodon paradoxus* (*Schizopora paradoxa* in MycoCosm).
  - Russulales: *Stereum hirsutum* (VibMO1 proxy; ShGE1) and *Heterobasidion irregulare*, both white rots; the ECM *Lactarius indigo* (L); and the saprotroph *Gloeopeniophorella convolvens*.
  - Gloeophyllales: two brown rots, *Gloeophyllum trabeum* and *Neolentinus lepideus* (NlePAD).
  - Jaapiales: *Jaapia argillacea*, with intermediate decay.
- **Boletales (7 genomes):**
  - the brown rots *Serpula lacrymans* S7.9, *Coniophora olivacea* and *Hydnomerulius pinastri*;
  - the ECM *Paxillus involutus* (C+, the strain of the Rineau 2012 study), *P. adelphus*, *Suillus luteus* and *Pisolithus tinctorius* (L).
- **Agaricales (16 genomes):**
  - white rots: *Pleurotus ostreatus* PC15 and *P. eryngii*, *Lentinula edodes*, *Armillaria mellea* and *Hypholoma sublateritium* (CE15);
  - brown rot: *Fistulina hepatica*;
  - intermediate decay: *Schizophyllum commune* and *Cylindrobasidium torrendii*;
  - litter and dung saprotrophs: *Agaricus bisporus* H97 (C+), *Gymnopus androsaceus* (C+), *Coprinopsis cinerea* AmutBmut (`Copcin2`, with `Copci1` as the CC1G_ shadow) and *Agrocybe pediades* (LiP-like POD anchors);
  - ECM: *Laccaria bicolor*, *Amanita muscaria* with its asymbiotic sister *A. thiersii*, and *Hebeloma cylindrosporum*.
  - *P. eryngii* is placed with the cohort's own species tree, because the Li et al. tip of that name sits in Russulales, probably because the genome behind it is misidentified.
- **Thelephorales:** the ECM *Thelephora terrestris*.

**Not sampled** (accepted gaps):
- Wallemiomycetes and Exobasidiomycetes;
- lichen-forming fungi;
- Microsporidia;
- Phallomycetidae beyond *Sphaerobolus* (Geastrales): Gomphales, Hysterangiales and Phallales. Atheliales rests on the presence-only *Piloderma*.

### 2.2 Lifestyles and decay modes

**White rot and class II peroxidases (PODs).**
- PODs are typed by catalytic residues, not by family membership.
- **Copy numbers.** White-rot species "have 5 to 26 copies (avg. 14) of genes encoding ligninolytic PODs, but all brown rot species lack these enzymes, as do the ectomycorrhizal (ECM) Laccaria bicolor, the soil saprotroph C. cinerea, and Schizophyllum commune" [@floudas2012Paleozoicoriginenzymatic, PDF p. 9].
- **Sampling.** The cohort samples PODs where they arose. The POD-free early lineages (`Dacsp1`, `Calco1`, `Botbo1`, `Pirin1`, `Sissu1`) are coded as ancestrally POD-free, not as losses. Above that node it samples these white rots:
  - Auriculariales;
  - two Hymenochaetales;
  - two Russulales;
  - Corticiales;
  - 14 Polyporales;
  - at least three Agaricales lineages.
- **Anchors:**
  - the *Pleurotus ostreatus* VP, genetically shown to matter for lignin degradation [@salame2014InactivationPleurotusostreatus, PDF p. 1–2];
  - *P. eryngii* VPL1-3;
  - the *Bjerkandera* VP and the *Agrocybe* LiP-like PODs, the main POD anchors outside Polyporales;
  - *P. chrysosporium* LiP/MnP [@martinez2004Genomesequencelignocellulose, PDF p. 3];
  - MnP6 of the selective delignifier *Gelatoporia* [@fernandez-fueyo2012ComparativegenomicsCeriporiopsis, PDF p. 1].
- **Counts.** The selective delignifier *Gelatoporia subvermispora* is P, so counts for that clade come from its C-tier sister *Obba rivulosa*.

**Brown rot.** Five origins are sampled, each with a sister of contrasting decay mode:

| origin | genomes |
|---|---|
| Dacrymycetes | `Calco1`, `Dacsp1` |
| Gloeophyllales | `Glotr1_1`, `Neole1` |
| Antrodia clade | `Fompi3`, `Pipbet1_1`, `PosplRSB12_1`, `Wolco1` |
| Boletales | `SerlaS7_9_2`, `Conol1`, `Hydpi2` |
| Agaricales | `Fishe1` |

- **What is lost, and what is kept.** Brown-rot Polyporales, Gloeophyllales and Dacrymycetales lack GH6 and GH7, whereas *Fistulina* and the saprotrophic Boletales keep complete or partial cellobiohydrolase sets [@floudas2015Evolutionnovelwood, PDF p. 11]. Crystalline-cellulose machinery is therefore coded as copy number, not presence.
- **Fistulina pseudogenes.** *Fistulina* shows ongoing pseudogenization [@floudas2015Evolutionnovelwood, PDF p. 1]. Its GFF models without a protein are kept as candidate pseudogenes and coded as fragments, never as absences.
- **Brown-rot model genomes:** *Postia placenta* (C+) lacks exo-cellobiohydrolases [@martinez2009Genometranscriptomesecretome, PDF p. 1, 5], and *Serpula* [@eastwood2011PlantCellWall, PDF p. 7–8].

**Soft rot and intermediate or uncertain decay.**
- **Soft rot in Ascomycota:**
  - *Kretzschmaria deusta* (`Kredeu1`, C) is the count-grade soft-rot genome;
  - *Daldinia eschscholtzii* and *Xylaria* cf. *heliscus* are Xylariales endophytes with soft-rot capacity;
  - *Chaetomium globosum* (P) gives presences only.
- **Intermediate decay in Basidiomycota:** *Schizophyllum commune*, *Cylindrobasidium torrendii* and *Jaapia argillacea*.
- **Uncertain decay:** *Botryobasidium*, *Plicaturopsis* and *Sistotremastrum*.
- ***Armillaria*** acquired at least 1,025 genes through 124 horizontal-transfer events, mainly from Ascomycota [@sahu2023Verticalhorizontalgene, PDF p. 2]. Domains of foreign origin are scored as originations.

**Litter and dung saprotrophs.**
- **Litter:**
  - *Agaricus bisporus*, the humic-litter decomposer with limited PODs [@morin2012Genomesequencebutton, PDF p. 1, 4];
  - *Gymnopus androsaceus*, whose litter-against-wood contrast within Omphalotaceae is made with *Lentinula*;
  - *Agrocybe pediades*, from grassland litter.
- **Dung:** *Coprinopsis cinerea*, *Podospora anserina* and *Ascobolus immersus*.
- **Soil:** *Morchella*, the Mortierellomycotina and the Mucorales give the soil saprotrophs.

**Ectomycorrhizal origins.**
- **The question.** ECM fungi have a reduced PCWDE set compared with their wood-decaying ancestors, but keep a unique array [@kohler2015Convergentlossesdecay, PDF p. 1]. Multiple ECM lineages lost most lignocellulose-degrading genes [@miyauchi2020Largescalegenomesequencing, PDF p. 2].
- **Pairing.** Each origin is paired with its nearest saprotroph, so that "loss in ECM" can be told apart from "absent in the clade":

| origin | ECM genome(s) | saprotrophic comparator |
|---|---|---|
| *Amanita* | `Amamu1` | `Amath1`, the asymbiotic sister [@hess2014TransposableElementDynamics, PDF p. 1–2] |
| *Laccaria* | `Lacbi2` | `Copcin2` (Psathyrellaceae) |
| *Hebeloma* | `Hebcy2` | Agaricales saprotrophs (`Agrped1`, `Hypsu1`) |
| Russulaceae | `Lacind1` (L) | `Glocon1`; ECM Russulaceae lost PCWDE genes against it [@looney2022Evolutionarytransitionectomycorrhizal, PDF p. 1] |
| Paxillineae | `Paxin1` (C+), `Paxru2` | `Hydpi2`; *Paxillus* keeps a trimmed brown-rot Fenton system [@rineau2012ectomycorrhizalfungusPaxillus, PDF p. 1–2] |
| Suillineae | `Suilu4` | Boletales brown rots |
| Sclerodermatineae | `Pisti2` (L) | Boletales brown rots |
| Thelephorales | `Theter1` | Polyporales, the nearest well-sampled saprotrophic order |
| Cantharellales | `Hydru2` (C+) | `Botbo1` |
| Atheliales | `Pilcr1` (P and L) | none; presences only |
| Pezizales | `Tubme1v2` | `Morco1` |
| Gloniaceae | `Cenge3` | `Glost2` |

- **Limits on the data:**
  - No counts come from `Lacind1`, `Pisti2` or `Pilcr1`, and `Pilcr1` is not a data point for ECM losses.
  - Absences in `Hydru2` and `Paxin1` need confirmation against the assembly.
- **Sugar transporters:** "The lack of invertase and sucrose transporter genes in most ectomycorrhizal fungi means that they are unable to use apoplastic sucrose released by their host plant" [@miyauchi2020Largescalegenomesequencing, PDF p. 12].

**Ericoid, orchid and arbuscular mycorrhizae, and endophytes.**
- **Ericoid mycorrhizae** (*Oidiodendron maius*, *Rhizoscyphus ericae*) keep a dual saprotrophic and biotrophic repertoire, making them a mycorrhizal contrast without reduction [@martino2018Comparativegenomicstranscriptomics, PDF p. 3].
- **Orchid mycorrhiza:** *Tulasnella calospora* (P) is the only one. *Tulasnella* keeps moderate to high numbers of wood-decay CAZymes [@nagy2016ComparativeGenomicsEarlyDiverging, PDF p. 5], so its presences are informative.
- **Arbuscular mycorrhiza:** *Rhizophagus* is the AM contrast (§2.1).
- **Endophytes:**
  - the members: `Pirin1`, `Lepor2`, `Dales1`, `Xylhel2`, `Macpha1` and `Metro1`;
  - *Serendipita indica* has an expanded set of CBM1-containing proteins [@zuccaro2011EndophyticLifeStrategies, PDF p. 5], and a cluster of three hexose transporters closest to those of *C. cinerea* [@zuccaro2011EndophyticLifeStrategies, PDF p. 6];
  - pectate lyases are enriched in detrimental root colonizers compared with beneficial ones [@mesny2021Geneticdeterminantsendophytism, PDF p. 7].

**Plant pathogens and obligate biotrophs.**
- **Necrotrophs and hemibiotrophs:**
  - `Bcine1`, `Sclsc1` (P), `Fusgr1` and `Fusox3_1`;
  - `Magor1`, `Colgr2` and `Verda1` (C+);
  - `Cocsa1`, `Stano2` (P), `Zymtr1`, `Clafu1`, `Macpha1`, `Valma1` and `Grocl1`.
- **Wood pathogens:** `Hetan2`, `Armmel1`, `Ganluc1`, `Fomme1` and `Kredeu1`.
- **Obligate biotrophs:**
  - *Blumeria*, in which "We found no canonical cellulose-, xylan-, or pectin-degrading enzymes" [@spanu2010GenomeExpansionGene, PDF p. 4];
  - the two rusts;
  - *Ustilago maydis*.
- **Losses** in `Blugr2` and absences in `Pgt_201_A1` (C+) need confirmation against the assembly.

**Animal and insect pathogens, and mycoparasites.** These were not a target lifestyle, but they come with anchors and help separate plant-biomass genes from general nutrition:
- animal pathogens: *Aspergillus fumigatus*, *Coccidioides immitis*, *Cryptococcus neoformans* and *Candida albicans*;
- insect pathogens: *Metarhizium robertsii*, *Akanthomyces lecanii* and *Conidiobolus coronatus*;
- mycoparasites: *Trichoderma atroviride* and *Tremella mesenterica*.

**Anaerobic gut fungi.**
- **Genomes.** Four genera are sampled:
  - *Piromyces finnis* (`Pirfi3`, C), the best-quality genome;
  - *Anaeromyces robustus* (C+);
  - *Caecomyces churrovis* (C+, non-rhizoidal);
  - *Neocallimastix californiae* (`Neosp1`), C+ with a duplication flag under the genus exemption (§1.5), with *N. cameroonii* and *N. lanati* as shadows.
- **Cellulosome biology.** Their dockerin and scaffoldin domains have no similarity to bacterial counterparts, although several catalytic domains came from gut bacteria by horizontal transfer [@haitjema2017partslistfungal, PDF p. 2]. Scaffoldins are conserved and unique to Neocallimastigomycota [@haitjema2017partslistfungal, PDF p. 2, 6].
- **Why *N. californiae* holds the slot.** It is the organism of the scaffoldin and dockerin anchors (MycoCosm IDs 673330 and 460692) [@haitjema2017partslistfungal, PDF p. 4, 6].
- **How they are searched and scored:**
  - Transporters are searched with the SWEET Pfam (PF03083) as well as the MFS sugar-porter family.
  - Domains of bacterial origin are pruned and scored as originations.
  - The anaerobes also serve as a negative control for oxygen-dependent ring cleavage.

**Yeasts.**
- **The transporter and sensor backbone** (§2.1). *Calocera* aside, the large Hxt/MAL radiations, the pentose-transporter gains and losses and the Snf3/Rgt2 sensors are all yeast or Pezizomycotina events.
- **Aromatic catabolism.** It is now in scope through *C. parapsilosis*, *Blastobotrys*, *Yarrowia*, *Brettanomyces* and *Cutaneotrichosporon*.
- **CAZymes.** The yeasts form the CAZyme-poor end of the gradient, and the Pucciniomycotina yeast *Rhodotorula* adds the pectin-sugar side.

### 2.3 Choices between alternatives

| question | choice | reason |
|---|---|---|
| Blastocladiomycota | *Catenaria anguillulae* over *Allomyces macrogynus* | one genome suffices for the phylum; *Allomyces* has 19,446 gene models against 12,804, consistent with its reported genome duplication [EK], which would weaken copy-number rooting |
| second genomes of *Kluyveromyces*, *Lactarius* and *Armillaria* | one genome each (`Klula1`, `Lacind1`, `Armmel1`) | the second genome had no role of its own |
| third or later genomes of *Penicillium*, *Aspergillus* and *Trichoderma* | not included | order-level gaps ranked first; the *T. harzianum* XYR1 anchor is placed through *T. reesei* |
| *Candida albicans* and *C. parapsilosis* | both | not equivalent: the gentisate/hydroxyquinol cluster is present in one and absent in the other |
| *Sclerotinia sclerotiorum* (P) beside *Botrytis cinerea* (C) | both | *Sclerotinia* carries its own CRE1 and SsShy1 anchors; counts come from *Botrytis* |
| *Paxillus* | both `Paxin1` (C+) and `Paxru2` (C) | *P. involutus* is the strain of the trimmed Fenton-system study; *P. adelphus* gives count-grade data for the ECM origin |
| *Neocallimastix* | *N. californiae* in the cohort; *N. cameroonii* and *N. lanati* as shadows | *N. californiae* is the organism of the scaffoldin and dockerin anchors |
| *Coprinopsis cinerea* | AmutBmut (C) over Okayama 7 (C+) | better tier; Okayama 7 kept as shadow for the CC1G_ IDs |
| *Cryptococcus* | *C. neoformans* H99 (C) over *C. deneoformans* JEC21 (C+) | the transporter accessions look like H99 IDs [EK]; better tier |
| *F. oxysporum* 4287 | v3 over v2 | same strain; v2 has about 10,000 more gene models and splits the XlnR anchor across several of them |
| *Blastobotrys adeninivorans* | reference annotation `Arxad1` over `Blaade_LS3_1` | same assembly, scores within two markers; `Arxad1` carries the Agdc1p annotation |
| *Dichomitus squalens* | CBS 464.89 over LYAD-421 | the study strain with the TF anchors; LYAD-421 kept as shadow for the CE15 anchors |
| *Trametes versicolor* | `Trave1` over ATCC 20869 | same tier; `Trave1` carries the anchor IDs |
| *Starmerella magnoliae* | used as P | its C+ score rests on a single comparator (§1.5) |
| presence- and losses-only genomes with unique roles (`Catan2`, `Conco1`, `Tulca1`, `Dekbr2`, `Hanpo2`, `Pilcr1`, `Sphst1`) | kept | no processed genome serves the same role; each is used only as its tier allows |

## 3. Summary by portal

- **Columns:**
  - C, D, F and M are BUSCO percentages (complete, duplicated, fragmented, missing) on the collapsed proteome.
  - LRC is the lineage-relative completeness score (§1.5).
- **Tier:** the use allowed in the study (C, C+, P, L, P+L). ᵈ marks a duplication flag: an exempt whole-genome or genus-wide duplication, or duplication confirmed in a second strain.
- **Comparator-sensitive:** the tier differs against a wider comparison group (§1.5).
- **Order:** cohort first, then shadows; within each, by phylum (early-diverging lineages first), class and order.
- **Role:** the genome's main role in the cohort, shortened.

| portal | species | phylum · class · order | lifestyle | tier | C | D | F | M | LRC | role / note |
|---|---|---|---|---|---|---|---|---|---|--- |
| `Spipu1` | *Spizellomyces punctatus DAOM BR117* | Chytridiomycota · Chytridiomycetes · Spizellomycetales | saprotroph (aquatic/soil chytrid) | **C** | 97.3 | 2.0 | 1.4 | 1.2 | 1102 | Chytrid root |
| `Hyacur1` | *Hyaloraphidium curvatum SAG235-1* | Monoblepharidomycota · Monoblepharidomycetes · Hyaloraphidiales | saprotroph (aquatic) [EK] | **C+** | 89.9 | 0.9 | 4.6 | 5.4 | 1044 | Monoblepharidomycetes |
| `Anasp1` | *Anaeromyces robustus* | Neocallimastigomycota · Neocallimastigomycetes · Neocallimastigales | anaerobic gut fungus | **C+** | 79.3 | 4.5 | 2.9 | 17.8 | 1054 | Anaerobic gut fungus genus 1 |
| `Caecom1` | *Caecomyces churrovis A* | Neocallimastigomycota · Neocallimastigomycetes · Neocallimastigales | anaerobic gut fungus | **C+** | 80.2 | 5.3 | 2.6 | 17.2 | 1062 | Fourth AGF genus (non-rhizoidal) |
| `Neosp1` | *Neocallimastix californiae G1* | Neocallimastigomycota · Neocallimastigomycetes · Neocallimastigales | anaerobic gut fungus | **C+**ᵈ | 80.1 | 36.3 | 3.1 | 16.8 | 1060 | Anaerobic gut fungus genus 2 |
| `Pirfi3` | *Piromyces finnis* | Neocallimastigomycota · Neocallimastigomycetes · Neocallimastigales | anaerobic gut fungus | **C** | 81.6 | 3.8 | 2.9 | 15.6 | 1077 | Anaerobic gut fungus genus 3 |
| `Catan2` | *Catenaria anguillulae PL171* | Blastocladiomycota · Blastocladiomycetes · Blastocladiales | zoosporic saprotroph/parasite | **P** | 72.9 | 4.4 | 9.2 | 17.9 | 881 | Blastocladiomycota outgroup (preferred over Allma1); comparator-sensitive |
| `Conco1` | *Conidiobolus coronatus NRRL28638* | Zoopagomycota · Entomophthoromycetes · Entomophthorales | insect pathogen; saprotroph | **P** | 80.5 | 1.3 | 5.6 | 13.9 | 926 | Zoopagomycota outgroup; comparator-sensitive |
| `Morel2` | *Linnemannia elongata AG-77* | Mucoromycota · Mortierellomycetes · Mortierellales | saprotroph (soil) | **C**ᵈ | 97.3 | 12.5 | 1.9 | 0.8 | 1097 | Mortierellomycotina soil saprotroph |
| `Rhimi59_2` | *Rhizopus microsporus ATCC11559* | Mucoromycota · Mucoromycetes · Mucorales | saprotroph; plant/animal pathogen | **C**ᵈ | 97.1 | 9.3 | 1.2 | 1.7 | 1096 | Mucorales |
| `RhiirDAOM197198_2` | *Rhizophagus irregularis DAOM197198* | Glomeromycota · Glomeromycetes · Glomerales | arbuscular mycorrhizal | **C** | 94.8 | 1.4 | 2.2 | 2.9 | 1085 | AM fungus |
| `Arxad1` | *Blastobotrys adeninivorans LS3* | Ascomycota · Dipodascomycetes · Dipodascales | yeast (soil) | **C** | 96.3 | 1.2 | 1.2 | 2.6 | 1112 | Agdc1p gallate decarboxylase and pyrogallol route |
| `Yarli1` | *Yarrowia lipolytica (strain CLIB122)* | Ascomycota · Dipodascomycetes · Dipodascales | yeast | **C+** | 92.0 | 1.2 | 2.2 | 5.8 | 1062 | 6 YHT anchors; comparator-sensitive |
| `Macpha1` | *Macrophomina phaseolina MPI-SDFR-AT-0080* | Ascomycota · Dothideomycetes · Botryosphaeriales | plant pathogen / root endophyte | **C** | 99.0 | 0.4 | 0.8 | 0.2 | 1111 | Botryosphaeriales |
| `AurpulNBB1` | *Aureobasidium pullulans NBB 7.2.1* | Ascomycota · Dothideomycetes · Dothideales | saprotroph / epiphyte | **C** | 99.7 | 0.4 | 0.1 | 0.2 | 1119 | Dothideales |
| `Clafu1` | *Fulvia fulva* | Ascomycota · Dothideomycetes · Mycosphaerellales | plant pathogen | **C** | 96.7 | 0.4 | 0.7 | 2.6 | 1085 | Apoplastic biotroph, sister to Zymoseptoria |
| `Zymtr1` | *Zymoseptoria tritici IP0323* | Ascomycota · Dothideomycetes · Mycosphaerellales | plant pathogen (latent phase) | **C** | 99.1 | 0.2 | 0.3 | 0.6 | 1112 | Latent apoplastic pathogen with reduced PCWDE repertoire |
| `Cenge3` | *Cenococcum geophilum 1.58* | Ascomycota · Dothideomycetes · Mytilinidiales | ECM | **C** | 98.0 | 1.3 | 1.1 | 0.9 | 1100 | Ascomycete ECM (Mytilinidiales) |
| `Glost2` | *Glonium stellatum CBS 207.34* | Ascomycota · Dothideomycetes · Mytilinidiales | saprotroph (wood/bark; decay ?) | **C** | 96.3 | 0.3 | 2.6 | 1.1 | 1081 | Saprotrophic sister of Cenococcum for ECM loss contrast |
| `Cocsa1` | *Bipolaris sorokiniana* | Ascomycota · Dothideomycetes · Pleosporales | plant pathogen | **C** | 98.7 | 0.3 | 0.8 | 0.5 | 1107 | Necrotroph |
| `Stano2` | *Parastagonospora nodorum* | Ascomycota · Dothideomycetes · Pleosporales | plant pathogen | **P** | 91.6 | 0.2 | 5.5 | 2.9 | 1028 | Necrotroph |
| `Capcor1` | *Capronia coronata CBS 617.96* | Ascomycota · Eurotiomycetes · Chaetothyriales | saprotroph (black yeast relative) | **C** | 98.6 | 0.0 | 0.9 | 0.5 | 1107 | CcFDC UbiD anchor |
| `Claps1` | *Cladophialophora psammophila CBS 110553* | Ascomycota · Eurotiomycetes · Chaetothyriales | saprotroph (black yeast) | **C** | 98.2 | 0.2 | 1.0 | 0.8 | 1103 | non-Eurotiales Eurotiomycete |
| `Aspac1` | *Aspergillus aculeatus ATCC16872* | Ascomycota · Eurotiomycetes · Eurotiales | saprotroph | **C** | 97.9 | 0.8 | 0.8 | 1.3 | 1101 | Only characterized AA16 (exact strain) |
| `Aspfu1` | *Aspergillus fumigatus Af293 from AspGD* | Ascomycota · Eurotiomycetes · Eurotiales | saprotroph; animal pathogen | **C** | 97.0 | 0.3 | 2.0 | 1.0 | 1091 | AfAA9, AfAA11A/B, GH7, GH28 anchors |
| `Aspni_NRRL3_1` | *Aspergillus niger NRRL3* | Ascomycota · Eurotiomycetes · Eurotiales | saprotroph | **C** | 99.6 | 0.5 | 0.2 | 0.3 | 1119 | Pectinase anchor set (Pga, Pel, Rhg), GH7, GOx |
| `Aspnid1` | *Aspergillus nidulans* | Ascomycota · Eurotiomycetes · Eurotiales | saprotroph | **C** | 98.0 | 0.2 | 1.1 | 1.0 | 1102 | AA13, AA9 C4, GH7, GH28, PL1 anchors |
| `Aspor1` | *Aspergillus oryzae RIB40* | Ascomycota · Eurotiomycetes · Eurotiales | saprotroph | **C+** | 92.0 | 1.2 | 2.2 | 5.8 | 1034 | MalP anchor |
| `Pench1` | *Penicillium chrysogenum* | Ascomycota · Eurotiomycetes · Eurotiales | saprotroph | **C** | 98.6 | 1.0 | 0.4 | 1.1 | 1109 | PrGE1 CE15 and GOx anchors |
| `Penoc1` | *Penicillium occitanis CL100* | Ascomycota · Eurotiomycetes · Eurotiales | saprotroph | **C** | 95.5 | 0.7 | 1.3 | 3.2 | 1073 | closest genome to T. cellulolyticus: TacA, TctA, XlnR |
| `Penox1` | *Penicillium oxalicum 114-2* | Ascomycota · Eurotiomycetes · Eurotiales | saprotroph | **C** | 99.4 | 0.5 | 0.3 | 0.4 | 1118 | GH7 Cbh1/Egl1 and GH45 anchors, exact strain 114-2 |
| `Theau2` | *Thermoascus aurantiacus ATCC26904* | Ascomycota · Eurotiomycetes · Eurotiales | saprotroph (thermophile) | **C** | 98.5 | 0.2 | 0.4 | 1.2 | 1108 | TaLPMO9A and TaCel7A anchors |
| `Cocim1` | *Coccidioides immitis RS* | Ascomycota · Eurotiomycetes · Onygenales | animal pathogen | **C** | 99.3 | 0.0 | 0.4 | 0.4 | 1116 | documented CLR-1 loss |
| `Bcine1` | *Botrytis cinerea B05.10* | Ascomycota · Leotiomycetes · Helotiales | necrotroph | **C** | 99.6 | 0.2 | 0.1 | 0.3 | 1118 | Necrotroph |
| `Blugr2` | *Blumeria hordei* | Ascomycota · Leotiomycetes · Helotiales | obligate biotroph | **C** | 98.4 | 0.8 | 1.0 | 0.6 | 1104 | Obligate biotroph without canonical PCWDEs |
| `Lepor2` | *Leptodontidium sp. MPI-SDFR-AT-0119* | Ascomycota · Leotiomycetes · Helotiales | root endophyte | **C** | 99.3 | 0.8 | 0.3 | 0.4 | 1114 | Helotiales root endophyte with CAZyme data (replaces Phisc1) |
| `Rhier1` | *Hyaloscypha hepaticicola* | Ascomycota · Leotiomycetes · Helotiales | ericoid mycorrhizal | **C** | 99.1 | 0.0 | 0.7 | 0.2 | 1112 | Second ERM fungus |
| `Sclsc1` | *Sclerotinia sclerotiorum* | Ascomycota · Leotiomycetes · Helotiales | necrotroph | **P** | 90.0 | 0.2 | 7.6 | 2.4 | 1010 | CRE1 anchor |
| `Oidma1` | *Oidiodendron maius Zn* | Ascomycota · Leotiomycetes · incertae sedis | ericoid mycorrhizal | **C** | 96.6 | 1.1 | 1.6 | 1.8 | 1084 | ERM fungus retaining saprotrophic CAZyme repertoire |
| `Lipst1_1` | *Lipomyces starkeyi NRRL Y-11557* | Ascomycota · Lipomycetes · Lipomycetales | yeast (soil, oleaginous) | **C+** | 90.3 | 0.6 | 5.0 | 4.7 | 1063 | Earliest Saccharomycotina lineage |
| `Artol1` | *Arthrobotrys oligospora ATCC 24927* | Ascomycota · Orbiliomycetes · Orbiliales | saprotroph; nematode-trapping | **C** | 97.6 | 0.1 | 1.6 | 0.8 | 1101 | Orbiliomycetes gap |
| `Ascim1` | *Ascobolus immersus RN42* | Ascomycota · Pezizomycetes · Pezizales | saprotroph (coprophilous) | **C+** | 94.5 | 1.6 | 2.2 | 3.3 | 1063 | Dung saprotroph |
| `Morco1` | *Morchella importuna CCBAS932* | Ascomycota · Pezizomycetes · Pezizales | saprotroph | **C** | 96.2 | 1.2 | 1.9 | 2.0 | 1082 | Pezizomycetes |
| `Tubme1v2` | *Tuber melanosporum Mel28* | Ascomycota · Pezizomycetes · Pezizales | ECM | **C** | 95.9 | 1.2 | 2.0 | 2.0 | 1079 | Truffle ECM with restricted PCWDE set |
| `Ambmon_NRRLY1484_1` | *Ambrosiozyma monospora NRRL Y-1484* | Ascomycota · Pichiomycetes · Pichiales | yeast | **C+** | 81.5 | 1.0 | 2.7 | 15.9 | 1061 | Pichiales |
| `Dekbr2` | *Brettanomyces bruxellensis CBS 2499 (MycoCosm: Dekkera)* | Ascomycota · Pichiomycetes · Pichiales | yeast | **P** | 73.6 | 2.2 | 7.1 | 19.3 | 985 | BHT1-3 anchors |
| `Hanpo2` | *Ogataea polymorpha NCYC 495 leu1.1* | Ascomycota · Pichiomycetes · Pichiales | yeast (methylotrophic) | **P** | 77.5 | 0.9 | 5.4 | 17.1 | 1015 | Pichiales |
| `Canalb1` | *Candida albicans SC5314* | Ascomycota · Pichiomycetes · Serinales | yeast; animal commensal/pathogen | **C** | 80.3 | 1.6 | 3.0 | 16.7 | 1068 | ITR1 anchor |
| `Canint_NRRLY981_1` | *Candida intermedia NRRL Y-981* | Ascomycota · Pichiomycetes · Serinales | yeast | **C** | 80.8 | 1.1 | 2.9 | 16.2 | 1076 | Gxf1/Gxs1 xylose transporter anchors |
| `Canpa1` | *Candida parapsilosis CDC317* | Ascomycota · Pichiomycetes · Serinales | commensal yeast | **C** | 80.6 | 1.8 | 2.9 | 16.5 | 1078 | Mnx1/2/3, Hdx1, Gdx1, Fph1 gentisate/hydroxyquinol cluster anchors |
| `Debha1` | *Debaryomyces hansenii CBS 767* | Ascomycota · Pichiomycetes · Serinales | yeast | **C** | 83.5 | 1.5 | 2.2 | 14.3 | 1104 | 6 polyol/inositol porter anchors |
| `Picst3` | *Scheffersomyces stipitis NRRL Y-11545* | Ascomycota · Pichiomycetes · Serinales | yeast (xylose-fermenting) | **C** | 80.7 | 1.6 | 3.7 | 15.5 | 1080 | 7 pentose/cellobiose ST anchors |
| `Klula1` | *Kluyveromyces lactis* | Ascomycota · Saccharomycetes · Saccharomycetales | yeast | **C** | 75.9 | 1.1 | 3.2 | 20.9 | 1090 | 5 ST anchors incl. LAC12 lactose and Rag4 sensor |
| `Sacce1` | *Saccharomyces cerevisiae S288C* | Ascomycota · Saccharomycetes · Saccharomycetales | yeast | **C**ᵈ | 75.7 | 6.2 | 2.8 | 21.6 | 1079 | Yeast PCWDE loss |
| `Zygrou1` | *Zygosaccharomyces rouxii NRRL Y-64007* | Ascomycota · Saccharomycetes · Saccharomycetales | yeast (osmotolerant) | **C** | 76.2 | 1.2 | 3.0 | 20.8 | 1092 | ZrFsy1 anchor |
| `Schpo1` | *Schizosaccharomyces pombe* | Ascomycota · Schizosaccharomycetes · Schizosaccharomycetales | yeast | **C** | 90.6 | 4.1 | 1.7 | 7.7 | 1088 | 7 Ght/itr anchors; comparator-sensitive |
| `Valma1` | *Valsa mali str. 03-8* | Ascomycota · Sordariomycetes · Diaporthales | plant pathogen (canker) | **C** | 97.4 | 0.3 | 0.6 | 2.0 | 1093 | VmHbh1 anchor |
| `Colgr2` | *Colletotrichum graminicola M1.001* | Ascomycota · Sordariomycetes · Glomerellales | plant pathogen (hemibiotroph) | **C** | 96.4 | 0.3 | 2.5 | 1.1 | 1082 | 4 Hxt plus Mbt1 maltose anchors |
| `Verda1` | *Verticillium dahliae VdLs.17* | Ascomycota · Sordariomycetes · Glomerellales | plant pathogen (vascular wilt) | **C+** | 93.4 | 0.3 | 3.5 | 3.1 | 1048 | Vascular wilt |
| `Corco1` | *Akanthomyces lecanii RCEF 1005 (MycoCosm: Lecanicillium)* | Ascomycota · Sordariomycetes · Hypocreales | insect pathogen | **C** | 98.8 | 0.4 | 0.4 | 0.9 | 1108 | Akanthomyces proxy (CYP1, DIO4) |
| `Fusgr1` | *Fusarium graminearum* | Ascomycota · Sordariomycetes · Hypocreales | plant pathogen | **C** | 97.1 | 0.4 | 2.0 | 0.9 | 1089 | Necrotroph |
| `Fusox3_1` | *Fusarium oxysporum f. sp. lycopersici 4287* | Ascomycota · Sordariomycetes · Hypocreales | plant pathogen (vascular wilt) | **C** | 98.8 | 0.6 | 0.4 | 0.9 | 1109 | F. oxysporum 4287: CRE1, XlnR anchors |
| `Metro1` | *Metarhizium robertsii ARSEF 23* | Ascomycota · Sordariomycetes · Hypocreales | insect pathogen; root endophyte | **C** | 98.2 | 0.9 | 0.7 | 1.1 | 1103 | Mrt sucrose/galactoside anchor |
| `Triat2` | *Trichoderma atroviride* | Ascomycota · Sordariomycetes · Hypocreales | saprotroph; mycoparasite | **C** | 98.5 | 0.4 | 1.1 | 0.4 | 1105 | XYR1 anchor |
| `Trire_Chr` | *Trichoderma reesei QM6a* | Ascomycota · Sordariomycetes · Hypocreales | saprotroph | **C** | 98.6 | 0.2 | 0.8 | 0.6 | 1106 | Core cellulase anchors GH7/GH6/GH45/AA9/CE15 |
| `Magor1` | *Pyricularia oryzae 70-15* | Ascomycota · Sordariomycetes · Magnaporthales | plant pathogen | **C** | 98.6 | 0.0 | 0.4 | 1.1 | 1106 | Hemibiotroph |
| `Grocl1` | *Grosmannia clavigera kw1407* | Ascomycota · Sordariomycetes · Ophiostomatales | plant pathogen (beetle-vectored) | **C** | 96.0 | 0.0 | 1.6 | 2.4 | 1077 | documented CLR-2 loss |
| `Chagl_1` | *Chaetomium globosum* | Ascomycota · Sordariomycetes · Sordariales | soft rot [EK] | **P** | 85.9 | 0.2 | 7.9 | 6.1 | 965 | Classic soft-rot organism (P tier) |
| `Neucr2` | *Neurospora crassa OR74A* | Ascomycota · Sordariomycetes · Sordariales | saprotroph | **C** | 99.7 | 0.1 | 0.2 | 0.1 | 1119 | Seven AA9, CDH-1/2, AA13, GH45 anchors |
| `Podan3` | *Podospora anserina S mat+* | Ascomycota · Sordariomycetes · Sordariales | saprotroph (coprophilous) | **C** | 98.0 | 0.1 | 1.3 | 0.7 | 1100 | Six PaLPMO9 and two CE15 anchors |
| `Spoth2` | *Thermothelomyces thermophilus* | Ascomycota · Sordariomycetes · Sordariales | saprotroph (thermophile) | **C** | 96.3 | 0.1 | 1.7 | 2.0 | 1082 | Largest AA9 anchor set |
| `Thite2` | *Thermothielavioides terrestris* | Ascomycota · Sordariomycetes · Sordariales | saprotroph (thermophile) | **C** | 97.1 | 0.3 | 1.5 | 1.4 | 1090 | Five TtLPMO9, Cel7A, CDH, GE anchors |
| `Dales1` | *Daldinia eschscholtzii FL0578* | Ascomycota · Sordariomycetes · Xylariales | endophyte (decay ?) | **C** | 99.3 | 0.0 | 0.4 | 0.3 | 1114 | Xylariales soft rot/endophyte |
| `HypEC38_3` | *Hypoxylon sp. EC38* | Ascomycota · Sordariomycetes · Xylariales | endophyte (decay ?) | **C** | 98.1 | 0.1 | 0.8 | 1.1 | 1101 | Xylariales |
| `Kredeu1` | *Kretzschmaria deusta IL1129* | Ascomycota · Sordariomycetes · Xylariales | soft rot; wood pathogen | **C** | 99.2 | 0.6 | 0.5 | 0.3 | 1113 | Xylariales soft-rot/white-rot-like wood pathogen |
| `Xylhel2` | *Xylaria cf. heliscus FL0509* | Ascomycota · Sordariomycetes · Xylariales | saprotroph | **C** | 99.1 | 0.3 | 0.4 | 0.4 | 1112 | Xylarialean endophyte |
| `Saico1` | *Saitoella complicata NRRL Y-17804* | Ascomycota · Taphrinomycetes · Taphrinales | yeast | **C** | 94.6 | 1.6 | 4.2 | 1.2 | 1068 | Second Taphrinomycotina; comparator-sensitive |
| `Stamag_NRRLY2024_1` | *Starmerella magnoliae NRRL Y-2024* | Ascomycota · Trigonopsidomycetes · Trigonopsidales | yeast (flower/insect) | **P** | 75.0 | 0.5 | 4.4 | 20.6 | 1021 | CmFsy1 anchor; used as P (LRC tier C+ rests on one comparator) |
| `Agabi_varbisH97_2` | *Agaricus bisporus var bisporus (H97)* | Basidiomycota · Agaricomycetes · Agaricales | litter decomposer (humicolous) | **C+** | 90.9 | 0.9 | 2.1 | 7.0 | 1020 | Humic litter decomposer with limited PODs |
| `Agrped1` | *Agrocybe pediades AH 40210* | Basidiomycota · Agaricomycetes · Agaricales | saprotroph (grassland litter) | **C** | 98.3 | 1.3 | 0.9 | 0.8 | 1103 | Agaricales LiP-like POD anchors, exact strain |
| `Amamu1` | *Amanita muscaria Koide* | Basidiomycota · Agaricomycetes · Agaricales | ECM | **C** | 96.9 | 1.4 | 1.9 | 1.2 | 1087 | ECM Amanita |
| `Amath1` | *Saproamanita thiersii* | Basidiomycota · Agaricomycetes · Agaricales | saprotroph (asymbiotic) | **C** | 97.7 | 0.6 | 1.3 | 1.0 | 1096 | Asymbiotic Amanita sister for ECM loss contrast |
| `Armmel1` | *Armillaria mellea ELDO17* | Basidiomycota · Agaricomycetes · Agaricales | white rot (soft-rot-like CAZymes); pathogen | **C** | 98.3 | 2.9 | 0.9 | 0.8 | 1103 | White-rot pathogen |
| `Copcin2` | *Coprinopsis cinerea AmutBmut* | Basidiomycota · Agaricomycetes · Agaricales | saprotroph (dung/litter) | **C** | 99.4 | 1.1 | 0.3 | 0.4 | 1115 | AA9 burst, AA12, generic peroxidase |
| `Cylto1` | *Cylindrobasidium torrendii FP15055* | Basidiomycota · Agaricomycetes · Agaricales | intermediate decayer | **C** | 98.9 | 1.8 | 0.6 | 0.4 | 1110 | Intermediate weak decay, soft-rot-like (Agaricales) |
| `Fishe1` | *Fistulina hepatica* | Basidiomycota · Agaricomycetes · Agaricales | brown rot | **C** | 97.5 | 2.5 | 1.2 | 1.2 | 1094 | Brown-rot origin in Agaricales |
| `Gyman1` | *Gymnopus androsaceus JB14* | Basidiomycota · Agaricomycetes · Agaricales | litter decomposer | **C+** | 94.7 | 4.6 | 2.8 | 2.6 | 1062 | Litter decomposer |
| `Hebcy2` | *Hebeloma cylindrosporum h7* | Basidiomycota · Agaricomycetes · Agaricales | ECM | **C** | 96.6 | 2.0 | 2.9 | 0.5 | 1084 | Third Agaricales ECM origin |
| `Hypsu1` | *Hypholoma sublateritium* | Basidiomycota · Agaricomycetes · Agaricales | white rot | **C** | 98.6 | 0.9 | 1.0 | 0.4 | 1106 | Strophariaceae white rot |
| `Lacbi2` | *Laccaria bicolor* | Basidiomycota · Agaricomycetes · Agaricales | ECM | **C** | 98.0 | 4.2 | 1.4 | 0.6 | 1099 | ECM model |
| `Lenedo1` | *Lentinula edodes Le(Bin) 0899 ss11* | Basidiomycota · Agaricomycetes · Agaricales | white rot | **C** | 98.8 | 2.1 | 0.7 | 0.4 | 1109 | Agaricales white rot (Omphalotaceae) |
| `Pleery1` | *Pleurotus eryngii ATCC 90797* | Basidiomycota · Agaricomycetes · Agaricales | white rot; weak root parasite [EK] | **C** | 98.8 | 1.4 | 0.7 | 0.4 | 1109 | VPL1-3, VPS1 and AAO anchors |
| `PleosPC15_2` | *Pleurotus ostreatus PC15* | Basidiomycota · Agaricomycetes · Agaricales | white rot | **C+** | 94.5 | 1.6 | 2.7 | 2.9 | 1060 | VP/MnP lineage |
| `Schco3` | *Schizophyllum commune H4-8* | Basidiomycota · Agaricomycetes · Agaricales | intermediate decayer | **C** | 97.2 | 1.2 | 1.8 | 1.0 | 1091 | Intermediate decay |
| `Plicr1` | *Plicaturopsis crispa* | Basidiomycota · Agaricomycetes · Amylocorticiales | wood decay (7 PODs) | **P+L** | 81.3 | 18.7 | 0.9 | 17.8 | 912 | Only Amylocorticiales |
| `Pilcr1` | *Piloderma olivaceum F 1598* | Basidiomycota · Agaricomycetes · Atheliales | ECM | **P+L** | 84.9 | 17.6 | 3.4 | 11.7 | 953 | Atheliales ECM |
| `Aurde3_1` | *Auricularia subglabra* | Basidiomycota · Agaricomycetes · Auriculariales | white rot | **C** | 97.6 | 3.1 | 1.8 | 0.6 | 1095 | Earliest class II peroxidase lineage |
| `Exigl1` | *Exidia glandulosa* | Basidiomycota · Agaricomycetes · Auriculariales | white rot | **L** | 96.5 | 16.7 | 2.3 | 1.2 | 1083 | Second Auriculariales |
| `Conol1` | *Coniophora olivacea MUCL 20566* | Basidiomycota · Agaricomycetes · Boletales | brown rot | **C** | 98.6 | 0.4 | 1.0 | 0.4 | 1107 | Boletales brown rot (better genome than Conpu1) |
| `Hydpi2` | *Hydnomerulius pinastri* | Basidiomycota · Agaricomycetes · Boletales | brown rot | **C** | 96.6 | 2.1 | 2.8 | 0.6 | 1085 | Brown-rot Boletales with full CBH set |
| `Paxin1` | *Paxillus involutus ATCC 200175* | Basidiomycota · Agaricomycetes · Boletales | ECM (Fenton chemistry) | **C+** | 93.0 | 2.9 | 4.9 | 2.0 | 1045 | ECM Paxillus |
| `Paxru2` | *Paxillus adelphus Ve08.2h10* | Basidiomycota · Agaricomycetes · Boletales | ECM | **C** | 98.3 | 0.7 | 1.4 | 0.3 | 1104 | ECM Paxillus (replaces Paxin1 on quality) |
| `Pisti2` | *Pisolithus tinctorius Marx 270* | Basidiomycota · Agaricomycetes · Boletales | ECM | **L** | 98.7 | 8.0 | 0.6 | 0.7 | 1107 | ECM Sclerodermatineae (L tier) |
| `SerlaS7_9_2` | *Serpula lacrymans S7.9* | Basidiomycota · Agaricomycetes · Boletales | brown rot | **C** | 97.6 | 0.5 | 2.1 | 0.3 | 1096 | Boletales brown rot (better strain than S7.3) |
| `Suilu4` | *Suillus luteus UH-Slu-Lm8-n1* | Basidiomycota · Agaricomycetes · Boletales | ECM | **C** | 98.3 | 1.2 | 0.8 | 0.9 | 1104 | ECM Suillineae |
| `Botbo1` | *Botryobasidium botryosum* | Basidiomycota · Agaricomycetes · Cantharellales | wood decay, no PODs | **C** | 98.6 | 0.6 | 0.9 | 0.5 | 1106 | 32 LPMOs and 3 CDH, no POD |
| `Hydru2` | *Hydnum rufescens s.l.* | Basidiomycota · Agaricomycetes · Cantharellales | ECM | **C+** | 93.0 | 1.2 | 4.1 | 2.9 | 1044 | ECM Cantharellales |
| `Tulca1` | *Tulasnella calospora AL13/4D* | Basidiomycota · Agaricomycetes · Cantharellales | orchid mycorrhizal | **P** | 84.3 | 3.9 | 10.6 | 5.1 | 946 | Only orchid-mycorrhizal genome |
| `Punst1` | *Punctularia strigosozonata* | Basidiomycota · Agaricomycetes · Corticiales | white rot | **C+** | 94.1 | 0.2 | 1.8 | 4.1 | 1056 | Corticiales white rot |
| `Sphst1` | *Sphaerobolus stellatus* | Basidiomycota · Agaricomycetes · Geastrales | saprotroph (wood/dung) | **P+L** | 89.8 | 6.9 | 4.8 | 5.3 | 1008 | Phallomycetidae saprotroph |
| `Glotr1_1` | *Gloeophyllum trabeum* | Basidiomycota · Agaricomycetes · Gloeophyllales | brown rot | **C** | 98.2 | 0.4 | 1.5 | 0.3 | 1102 | Gloeophyllales brown rot |
| `Neole1` | *Neolentinus lepideus* | Basidiomycota · Agaricomycetes · Gloeophyllales | brown rot | **C** | 98.7 | 1.1 | 0.9 | 0.4 | 1107 | NlePAD sinapate decarboxylase anchor |
| `Fomme1` | *Fomitiporia mediterranea* | Basidiomycota · Agaricomycetes · Hymenochaetales | white rot; plant pathogen | **C+** | 95.0 | 0.8 | 1.1 | 3.9 | 1066 | Hymenochaetales white rot |
| `Schpa1` | *Xylodon paradoxus KUC8140* | Basidiomycota · Agaricomycetes · Hymenochaetales | white rot | **C** | 98.0 | 1.5 | 1.0 | 1.0 | 1100 | Hymenochaetales white rot |
| `Jaaar1` | *Jaapia argillacea* | Basidiomycota · Agaricomycetes · Jaapiales | wood decay, no or 1 POD | **C** | 98.8 | 1.3 | 0.6 | 0.6 | 1108 | Intermediate decay |
| `Bjead1_1` | *Bjerkandera adusta* | Basidiomycota · Agaricomycetes · Polyporales | white rot | **C** | 98.1 | 0.5 | 1.0 | 0.9 | 1101 | Only Polyporales VP anchors |
| `Cersu1` | *Gelatoporia subvermispora* | Basidiomycota · Agaricomycetes · Polyporales | white rot (selective) | **P** | 89.9 | 1.2 | 1.1 | 9.0 | 1009 | Selective delignifier |
| `Dicsqu464_2` | *Dichomitus squalens CBS 464.89* | Basidiomycota · Agaricomycetes · Polyporales | white rot | **C** | 99.4 | 1.6 | 0.6 | 0.0 | 1115 | User model |
| `Fomfom1` | *Fomes fomentarius CIRM-BRFM 1821* | Basidiomycota · Agaricomycetes · Polyporales | white rot | **C** | 98.9 | 1.3 | 0.7 | 0.4 | 1110 | Largest Polyporaceae FMO expansion (x11) |
| `Fompi3` | *Fomitopsis schrenkii* | Basidiomycota · Agaricomycetes · Polyporales | brown rot | **C** | 98.2 | 0.6 | 1.3 | 0.4 | 1102 | Antrodia-clade brown rot |
| `Ganluc1` | *Ganoderma lucidum G.260125-1* | Basidiomycota · Agaricomycetes · Polyporales | white rot; plant pathogen | **C** | 95.7 | 1.5 | 0.5 | 3.7 | 1074 | sister of Dichomitus |
| `Obbri1` | *Obba rivulosa 3A-2* | Basidiomycota · Agaricomycetes · Polyporales | white rot (selective) | **C** | 99.4 | 0.4 | 0.4 | 0.2 | 1115 | MnpA/B anchors |
| `Panru1` | *Panus rudis PR-1116 ss-1* | Basidiomycota · Agaricomycetes · Polyporales | white rot | **C** | 98.1 | 1.3 | 1.4 | 0.4 | 1101 | Panus/Lentinus clade |
| `Phaca1` | *Phanerochaete carnosa HHB-10118-Sp* | Basidiomycota · Agaricomycetes · Polyporales | white rot | **P** | 85.8 | 0.8 | 4.2 | 10.0 | 963 | PcGE1 CE15 anchor (P tier) |
| `Phchr2` | *Phanerochaete chrysosporium RP-78* | Basidiomycota · Agaricomycetes · Polyporales | white rot | **C** | 97.2 | 0.4 | 2.0 | 0.8 | 1091 | LiP/MnP, CDH, GLX, P2O, GH7, GH45 anchors |
| `Phlgi1` | *Phlebiopsis gigantea* | Basidiomycota · Agaricomycetes · Polyporales | white rot | **C** | 96.3 | 0.3 | 2.0 | 1.8 | 1080 | Phlebioid conifer pioneer |
| `Phlrad1` | *Phlebia radiata Fr. (isolate 79, FBCC0043)* | Basidiomycota · Agaricomycetes · Polyporales | white rot | **C** | 98.6 | 0.9 | 0.8 | 0.6 | 1106 | Lip1 anchor |
| `Pipbet1_1` | *Fomitopsis betulina CIRM-BRFM 1772* | Basidiomycota · Agaricomycetes · Polyporales | brown rot | **C** | 97.4 | 2.3 | 0.6 | 2.0 | 1093 | Second Antrodia-clade brown rot |
| `PosplRSB12_1` | *Rhodonia placenta* | Basidiomycota · Agaricomycetes · Polyporales | brown rot | **C+** | 94.1 | 0.3 | 3.5 | 2.4 | 1056 | Brown-rot model lacking exo-cellobiohydrolases |
| `Pyccin1` | *Trametes cinnabarina* | Basidiomycota · Agaricomycetes · Polyporales | white rot | **C** | 98.4 | 0.7 | 0.7 | 0.9 | 1104 | Vanillin/ferulate routes |
| `Pycco1` | *Trametes coccinea* | Basidiomycota · Agaricomycetes · Polyporales | white rot | **C** | 99.3 | 0.7 | 0.4 | 0.3 | 1114 | AA14 PcAA14A/B anchors, exact strain |
| `Trave1` | *Trametes versicolor* | Basidiomycota · Agaricomycetes · Polyporales | white rot | **C** | 95.5 | 0.9 | 0.5 | 3.9 | 1072 | Lpg1 LiP anchor |
| `Wolco1` | *Wolfiporia cocos MD-104 SS10* | Basidiomycota · Agaricomycetes · Polyporales | brown rot | **C** | 98.0 | 1.1 | 1.3 | 0.7 | 1099 | Third Antrodia-clade brown rot |
| `Glocon1` | *Gloeopeniophorella convolvens OM19405* | Basidiomycota · Agaricomycetes · Russulales | saprotroph (intermediate / ?) | **C** | 96.9 | 1.2 | 1.8 | 1.3 | 1089 | Saprotrophic counterpart of ECM Russulaceae |
| `Hetan2` | *Heterobasidion irregulare* | Basidiomycota · Agaricomycetes · Russulales | white rot; conifer pathogen | **C** | 96.4 | 0.8 | 2.2 | 1.3 | 1084 | Conifer white-rot pathogen |
| `Lacind1` | *Lactarius indigo 2018DUKE089* | Basidiomycota · Agaricomycetes · Russulales | ECM | **L** | 96.4 | 8.1 | 2.0 | 1.6 | 1084 | ECM Russulaceae |
| `Stehi1` | *Stereum hirsutum FP-91666 SS1* | Basidiomycota · Agaricomycetes · Russulales | white rot | **C** | 95.7 | 0.9 | 1.0 | 3.3 | 1076 | Proxy for Boreostereum VibMO1 |
| `Pirin1` | *Serendipita indica* | Basidiomycota · Agaricomycetes · Sebacinales | root endophyte | **C+** | 93.5 | 0.9 | 4.7 | 1.8 | 1049 | Root endophyte with CBM1 expansion |
| `Theter1` | *Thelephora terrestris UH-Tt-Lm1* | Basidiomycota · Agaricomycetes · Thelephorales | ECM | **C** | 99.2 | 1.6 | 0.4 | 0.4 | 1113 | ECM Thelephorales |
| `Sissu1` | *Sistotremastrum suecicum* | Basidiomycota · Agaricomycetes · Trechisporales | saprotroph (wood; decay uncertain) | **C** | 98.1 | 0.6 | 1.1 | 0.8 | 1101 | Trechisporales, POD-less early Agaricomycete |
| `Rhomi1` | *Cystobasidium minutum MCA 4210* | Basidiomycota · Cystobasidiomycetes · Cystobasidiales | yeast | **C** | 95.9 | 0.1 | 1.9 | 2.2 | 1087 | R. minuta P450rm (CYP53) |
| `Calco1` | *Calocera cornea* | Basidiomycota · Dacrymycetes · Dacrymycetales | brown rot | **C** | 95.3 | 0.4 | 2.0 | 2.8 | 1084 | 66-ST expansion, the largest known in Basidiomycota |
| `Dacsp1` | *Dacryopinax primogenitus DJM 731 SSP1* | Basidiomycota · Dacrymycetes · Dacrymycetales | brown rot (no PODs) | **C+** | 93.1 | 0.4 | 2.5 | 4.4 | 1060 | Dacrymycetes brown rot |
| `Rhoto_IFO0880_4` | *Rhodotorula toruloides* | Basidiomycota · Microbotryomycetes · Sporidiobolales | yeast | **C** | 96.8 | 0.5 | 1.7 | 1.5 | 1094 | RTO4_9841 galacturonate transporter (genetic evidence) |
| `Mellp3_H0` | *Melampsora larici-populina 98AG31 haplotype 0* | Basidiomycota · Pucciniomycetes · Pucciniales | obligate biotroph | **C** | 94.1 | 1.5 | 3.5 | 2.4 | 1078 | Obligate-biotroph rust |
| `Pgt_201_A1` | *Puccinia graminis f. sp. tritici 21-0 haplotype A* | Basidiomycota · Pucciniomycetes · Pucciniales | obligate biotroph | **C+** | 90.7 | 2.1 | 3.9 | 5.3 | 1039 | Proxy for rust HXT1 (Uromyces, P. striiformis) |
| `Cryne_H99_1` | *Cryptococcus neoformans var. grubii H99* | Basidiomycota · Tremellomycetes · Tremellales | yeast; animal pathogen | **C** | 93.9 | 0.4 | 2.8 | 3.3 | 1070 | ITR expansion (11 genes) and 2 anchors; comparator-sensitive |
| `Treme1` | *Tremella mesenterica Fries* | Basidiomycota · Tremellomycetes · Tremellales | mycoparasite | **C+** | 93.5 | 0.5 | 3.0 | 3.5 | 1065 | Sister for the Cryptococcus ITR expansion (3 vs 11 genes); comparator-sensitive |
| `Triol1` | *Cutaneotrichosporon oleaginosum IBC0246 (MycoCosm: Trichosporon)* | Basidiomycota · Tremellomycetes · Trichosporonales | yeast (oleaginous) | **C** | 96.4 | 0.5 | 1.7 | 1.9 | 1098 | Nearest portal to T. cutaneum (PhhY, TcMLE, salicylate hydroxylase) |
| `Ustma2_2` | *Ustilago maydis 521* | Basidiomycota · Ustilaginomycetes · Ustilaginales | plant pathogen (biotroph) | **C** | 98.3 | 0.1 | 0.7 | 1.0 | 1109 | Biotroph with few PCWDEs |
| `Neocon1` | *Neocallimastix cameroonii var. constans G3* | Neocallimastigomycota · Neocallimastigomycetes · Neocallimastigales | anaerobic gut fungus | shadow | 84.4 | 39.6 | 1.5 | 14.1 | 1106 | second Neocallimastix species |
| `Neolan1` | *Neocallimastix lanati* | Neocallimastigomycota · Neocallimastigomycetes · Neocallimastigales | anaerobic gut fungus | shadow | 84.8 | 41.3 | 0.9 | 14.3 | 1111 | third Neocallimastix species |
| `Morel_U14_1` | *Linnemannia elongata* | Mucoromycota · Mortierellomycetes · Mortierellales | root-associated | shadow | 96.2 | 13.8 | 1.8 | 2.0 | 1090 | Independent annotation supporting the real-duplication exemption |
| `Blaade_LS3_1` | *Blastobotrys adeninivorans LS3* | Ascomycota · Dipodascomycetes · Dipodascales | yeast (soil) | shadow | 96.6 | 1.3 | 1.2 | 2.2 | 1114 | Second annotation of the Arxad1 assembly |
| `Blugra3_16` | *Blumeria graminis f. sp. tritici 96224* | Ascomycota · Leotiomycetes · Helotiales | obligate biotroph | shadow | 86.0 | 0.1 | 5.3 | 8.6 | 965 | Confirms ST losses in Blumeria |
| `Fusox2` | *Fusarium oxysporum f. sp. lycopersici 4287* | Ascomycota · Sordariomycetes · Hypocreales | plant pathogen (vascular wilt) | shadow | 96.5 | 0.7 | 2.1 | 1.3 | 1084 | Older annotation carrying anchor IDs |
| `Agabi_varbur_1` | *Agaricus bisporus var. burnettii JB137-S8* | Basidiomycota · Agaricomycetes · Agaricales | litter decomposer (humicolous) | shadow | 91.2 | 0.9 | 2.0 | 6.9 | 1023 | Second Agaricus variety |
| `Copci1` | *Coprinopsis cinerea* | Basidiomycota · Agaricomycetes · Agaricales | saprotroph (dung/litter) | shadow | 95.0 | 0.8 | 2.8 | 2.2 | 1066 | Okayama 7 annotation carrying the CC1G_ IDs of the literature [EK] |
| `PleosPC9_1` | *Pleurotus ostreatus PC9* | Basidiomycota · Agaricomycetes · Agaricales | white rot | shadow | 93.7 | 1.5 | 4.2 | 2.1 | 1051 | Strain of CRE1, VP and MnP anchor IDs |
| `SerlaS7_3_2` | *Serpula lacrymans S7.3* | Basidiomycota · Agaricomycetes · Boletales | brown rot | shadow | 91.8 | 0.7 | 2.4 | 5.8 | 1031 | Confirms brown-rot ST absences |
| `Dicsq1` | *Dichomitus squalens LYAD-421 SS1* | Basidiomycota · Agaricomycetes · Polyporales | white rot | shadow | 95.5 | 0.7 | 0.8 | 3.7 | 1071 | Strain of DsGE1/2 and Floudas 2012 annotation |
| `Pycci1` | *Trametes cinnabarina* | Basidiomycota · Agaricomycetes · Polyporales | white rot | shadow | 86.5 | 0.3 | 6.3 | 7.2 | 970 | P tier |
| `Traver1` | *Trametes versicolor ATCC 20869* | Basidiomycota · Agaricomycetes · Polyporales | white rot | shadow | 98.2 | 0.7 | 1.0 | 0.8 | 1102 | Second genome of Trave1 |
| `Mellp3_H1` | *Melampsora larici-populina 98AG31 haplotype 1* | Basidiomycota · Pucciniomycetes · Pucciniales | obligate biotroph | shadow | 94.2 | 1.3 | 3.4 | 2.4 | 1079 | Second haplotype |
| `Pucgr2` | *Puccinia graminis f. sp. tritici* | Basidiomycota · Pucciniomycetes · Pucciniales | obligate biotroph | shadow | 89.1 | 8.0 | 3.4 | 7.5 | 1028 | Presence checks in a dikaryon |
| `Cryne_JEC21_1` | *Cryptococcus deneoformans* | Basidiomycota · Tremellomycetes · Tremellales | yeast; animal pathogen | shadow | 91.7 | 0.5 | 3.5 | 4.8 | 1044 | JEC21 reference annotation; second species of the ITR comparison |
