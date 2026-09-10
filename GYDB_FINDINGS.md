# GyDB Gag PrLD analysis — findings and interpretations

**Date:** 2026-08-09
**Status:** exploratory. None of this is in `claude/DECISIONS.md`, so under this
project's own convention nothing here is citable as established. Every number is
reproducible from the scripts listed at the end.

This is a side analysis on GyDB, not the REXdb pipeline. It exists because GyDB
publishes whole-element protein cores across many superfamilies, which lets the
same question be asked outside Viridiplantae Ty1/copia.

---

## 1. Scope and dataset

The GyDB *cores* database holds **2,636 protein cores in 74 classes**, covering
**778 distinct element names** — 440 classified elements, 97 superlineage
(`CopiaSL_`/`GypsySL_`) entries, and 242 accession-named sequences that are not
classified elements.

**333** elements carry both a GAG and an RT core. **272** of those have no `X` in
either, and are the analysis set throughout.

| superfamily | elements |
| --- | --- |
| Ty3/Gypsy | 107 |
| Ty1/Copia | 93 |
| Retroviridae | 49 |
| Bel/Pao | 23 |

Superfamily assignments come from the GyDB elements sheet
(`GyDB Elements - Superfamily & Host Database`). 231 of 272 matched by name; the
remaining 41 are `CopiaSL_`/`GypsySL_` entries that declare their own
superfamily, plus `GALV`/`GaLV` (case) and `Calypso5-1`/`Calypso` (version
suffix). `FFV` is filed as `Retrovirus` in the sheet and was merged into
`Retroviridae`. Caulimoviridae never enter the set — pararetroviruses carry
`GAGCOAT`, not `GAG`.

**Design principle, unchanged from the pipeline:** RT is aligned to build the
phylogeny, Gag is scored to make the trait, and the two are joined by element
name. The tree is never built from Gag.

## 2. Methods

| stage | tool and settings |
| --- | --- |
| align | MAFFT v7.525 L-INS-i (`--localpair --maxiterate 1000`) |
| trim | trimAl 1.5 `-automated1` — 678 → **131 columns**, 129 parsimony-informative |
| model | IQ-TREE ModelFinder over LG, WAG, JTT, VT, rtREV, HIVb, HIVw, Q.pfam, Dayhoff, Blosum62 → **LG+R7** by BIC |
| tree | IQ-TREE 3.1.3 maximum likelihood; logL −43,971.62 |
| support | 1000 UFBoot + 1000 SH-aLRT |
| prion | PLAAC, α = 0.5 and α = 1.0, core length 60 (L-4, L-5) |
| disorder | metapredict 3.0.2, its own IDR decomposition (L-6) |
| zinc knuckle | regex `C-X2-C-X4-H-X4-C`, the pattern used in F-2 |

ModelFinder was restricted to a curated candidate set because the full 1,232-model
sweep made no progress in 12 minutes. rtREV was included as the domain-appropriate
candidate (estimated from retroviral RT); LG still won.

**Tree quality.** Median SH-aLRT 86.1, median UFBoot 95; **43%** of 269 internal
nodes clear both SH-aLRT ≥ 80 and UFBoot ≥ 95. Terminal clades are well resolved;
the deep backbone is not, at 131 columns. Read clade membership, not deep
branching order.

**Two structural results from the tree itself.** All four superfamilies come out
as contiguous blocks without being told they exist, which is independent evidence
the thin alignment is recovering real signal. And **Ty3/Gypsy is paraphyletic**
under every rooting tested — Retroviridae and Bel/Pao nest inside it — the
conventional view of retroviral origins, recovered here independently.

---

## 3. Findings

### 3.1 Prion-like character is confined to Gag

All **26** called PrLDs in the entire GyDB core proteome fall in `GAG`. Every
other protein class is empty, including 2,172 non-Gag X-free proteins.

Length-stratified, so the confound is controlled:

| class | 0–200 aa | 200–400 | 400–600 | 600+ |
| --- | --- | --- | --- | --- |
| GAG | 0/32 | 3/188 | **16/89** | **7/27** |
| INT | 0/23 | 0/334 | 0/75 | 0/5 |
| RT | 0/24 | 0/334 | 0/7 | — |
| AP | 0/423 | 0/45 | 0/29 | 0/3 |
| ENV | 0/2 | 0/32 | 0/20 | 0/16 |

At 400–600 aa, Gag is 18% and the 159 non-Gag proteins in the same band are 0%.

Gag's *median* LLR is not elevated — it sits mid-pack, below ENV. What separates
Gag is the right tail.

### 3.2 Domain position splits by superfamily

Of the 25 domains inside the 272-element set:

| superfamily | tips | PrLDs | N-terminal | middle | C-terminal | rate |
| --- | --- | --- | --- | --- | --- | --- |
| Ty3/Gypsy | 107 | 17 | **0** | 8 | 9 | 15.9% |
| Ty1/Copia | 93 | 3 | **2** | 0 | 1 | 3.2% |
| Retroviridae | 49 | 5 | **0** | 0 | 5 | 10.2% |
| Bel/Pao | 23 | 0 | 0 | 0 | 0 | 0% |

**Both N-terminal domains in the database are Ty1/copia** — `Ty1B` (midpoint 25%)
and `Tkm1` (32%), adjacent tips in the yeast Pseudovirus clade. All five
retroviral domains are C-terminal at 81–95%. The third copia hit, `CoDi5.6`, is
C-terminal at 68.8%, 1.8 points past the boundary and effectively a borderline
call.

Position is the midpoint of the called domain as a fraction of the Gag core:
N-terminal < 33%, middle 33–67%, C-terminal > 67%.

### 3.3 PLAAC detects Q/N *in the absence of charge*

Composition of the 26 called domains against the other 624,793 residues in the
database:

| | PrLD | rest | diff |
| --- | --- | --- | --- |
| Q | 16.9% | 4.5% | +12.4 |
| N | 14.4% | 4.6% | +9.7 |
| P | 14.6% | 4.9% | +9.6 |
| G | 8.9% | 5.7% | +3.2 |
| Y | 6.4% | 3.3% | +3.1 |
| **Q+N** | **31.3%** | **9.1%** | **+22.2** |
| D+E | 1.7% | 11.2% | −9.4 |
| K+R | 9.1% | 12.9% | −3.9 |

Charge, not Q/N abundance, is what separates a call from a non-call. Decomposing
each Pseudovirus element's best window:

| | window | gain from Q+N | cost from K/R/D/E | total LLR |
| --- | --- | --- | --- | --- |
| Ty1B | 80–139 | +17.5 | −0.7 | **+16.9** |
| Ty2 | 73–132 | +17.5 | −2.7 | +13.0 |
| Tkm1 | 110–169 | +14.5 | −3.9 | +13.0 |
| Tse1 | 78–137 | +22.2 | −7.0 | +5.2 |
| Ty4 | 350–409 | +17.8 | **−26.4** | **−18.3** |
| Tdh2 | 270–329 | **+28.6** | **−25.5** | −6.2 |
| pCal | 180–239 | +12.0 | −13.0 | −20.5 |

`Tdh2` has the highest Q+N of all seven and still fails. PLAAC's per-residue
log-odds make this unavoidable: N +1.57, Q +1.28, but E −2.09, K −1.49, D −1.20.
One glutamate cancels 1.6 asparagines.

The N-terminal 120 aa of `Ty4` and `Tdh2` are 28–30% ILVFM with only 4–5% P+G —
folded, charged, helical-looking, not low-complexity. The four positive elements
are 15–20% P+G with 3–7% K+R.

**A second, independent failure mode:** a domain is called only when the HMM's
Viterbi parse yields a contiguous prion-state run of at least the 60 aa core
length. `Ty2` scores +13.0, identical to `Tkm1`, but its longest run is 54 — so
no call. This makes the presence/absence split a threshold effect, not a clean
binary in the underlying data.

### 3.4 Positional conservation without sequence conservation

The four positive Pseudovirus elements all place their LLR maximum in the same
relative window (24–33% of the protein, residues 73–169), and the three negatives
all peak at 65–92%. Best N-half vs best C-half window: `Ty1B` +16.9 vs −4.3,
`Ty2` +13.0 vs −5.8, `Tkm1` +13.0 vs −8.1, `Tse1` +5.2 vs −15.3.

But pairwise identity over a MAFFT alignment is only 46% for `Ty1B`/`Ty2` and
**19–28%** for everything else — barely above the twilight zone. Same location,
different sequence. `Ty1B` and `Ty2` open `MESQQLS/MESQQLH…`; `Tkm1` opens
`MASNDIIST…` with no shared N-terminal motif.

### 3.5 Zinc knuckles and PrLDs are anti-correlated — REVISED TWICE (2026-08-13)

> **The finding stands; the numbers below do not.** The knuckle calls in this
> subsection came from a regex, which over-calls. A first revision rescored them
> with Pfam gathering thresholds and concluded the effect was unsupported — that
> revision was itself wrong, because Pfam GA misses 39% of real knuckles and the
> attenuation dragged the odds ratio toward 1. The current statement is the
> **second revision**, on InterPro calls, at the end of this section. Both earlier
> versions are kept so the sequence is auditable. **Cite only the second
> revision.**

CCHC (`C-X2-C-X4-H-X4-C`) is present in 149 of 272 Gag cores. Zero C2H2 zinc
fingers anywhere.

| superfamily | n | CCHC | PrLD | both |
| --- | --- | --- | --- | --- |
| Ty3/Gypsy | 107 | 44 (41%) | 17 | 3 |
| Ty1/Copia | 93 | 46 (49%) | 3 | 0 |
| Retroviridae | 49 | 40 (82%) | 5 | 2 |
| Bel/Pao | 23 | 19 (83%) | 0 | 0 |
| **total** | **272** | **149 (55%)** | **25** | **5** |

CCHC is in 58% of Gags without a PrLD and 20% of those with one.

**Adjusting for length strengthens the association.** Length is strongly related
to PrLD (Wilcoxon p = 2×10⁻⁷; median core 470 aa vs 334) but barely to CCHC
(p = 0.052; 359 vs 347 aa), so it was masking the effect rather than creating it.

| approach | effect | p |
| --- | --- | --- |
| raw 2×2 | — | 0.0003 (Fisher) |
| Mantel–Haenszel across length quartiles | common OR **0.144** (95% CI 0.039–0.430) | 0.0001 |
| logistic, PrLD ~ CCHC | OR 0.179 (0.065–0.492) | 0.0009 |
| logistic, + log(length) | OR **0.100** (0.032–0.312) | 0.0001 |
| logistic, + log(length) + superfamily | OR 0.100 (0.029–0.345) | 0.0003 |

Within-stratum, at 354–461 aa: 1/38 of CCHC-bearing Gags have a PrLD against 9/30
without. Length itself is a strong positive predictor (OR 70.8 per e-fold), which
is the F-18 architecture effect reappearing.

Only 1 of the 5 co-occurring cases (`ASSBSV`) has its knuckle inside the called
domain; in the other four they are spatially separate.

#### Revision 1 (2026-08-13): rescored with Pfam models — ITSELF SUPERSEDED

> Kept for the record. The conclusion drawn here — that the effect does not
> survive — was an artefact of the Pfam gathering threshold missing 39% of
> knuckles. See revision 2.

**What was wrong.** `C-X2-C-X4-H-X4-C` is essentially the classic PROSITE-style
consensus — fourteen residues, four of them fixed. In a 300–500 aa protein it
matches by chance often enough to swamp the real signal. Rescored against the
six zinc-knuckle families of Pfam clan CL0511 (`hmmsearch --cut_ga`):

| superfamily | n | regex | Pfam GA |
| --- | --- | --- | --- |
| Ty1/Copia | 93 | 46 (49%) | **19 (20%)** |
| Ty3/Gypsy | 107 | 44 (41%) | **20 (19%)** |
| Retroviridae | 49 | 40 (82%) | **38 (78%)** |
| Bel/Pao | 23 | 19 (83%) | **1 (4%)** |
| **total** | **272** | **149 (55%)** | **78 (29%)** |

Retroviridae barely moves — retroviral NC knuckles are canonical, so the profile
finds them. Bel/Pao falls from 83% to 4%. That contrast is the clearest evidence
that the regex was matching noise rather than divergent knuckles.

Sensitivity was checked rather than assumed: the whole clan was used, not
PF00098 alone, and only PF00098 (plus 7 hits from PF13696) fires anywhere.
Relaxing the threshold to E < 1e-3 changes no conclusion below.

**The association after rescoring.**

| test | regex (as published above) | Pfam GA |
| --- | --- | --- |
| all 272, raw 2×2 | OR 0.179, p = 0.0003 | OR 0.445, **p = 0.17 (n.s.)** |
| all 272, + log(length) | OR 0.100 (0.032–0.312), p = 7e-05 | OR 0.161 (0.046–0.559), p = 0.004 |
| Ty3/Gypsy, + log(length) | OR 0.142 (0.031–0.655), p = 0.012 | OR 0.210 (0.032–1.373), **p = 0.10 (n.s.)** |
| Retroviridae, + log(length) | OR 0.032, p = 0.031 | OR 0.066 (0.004–0.978), p = 0.048 — 5 events |
| Ty1/Copia | 3 events, nothing | 3 events, nothing |

The 20%-vs-58% contrast that headlined the finding becomes 16% vs 30%. Only the
pooled length-adjusted model survives, and it pools four superfamilies whose
knuckle rates run from 4% to 78% — which is a between-group comparison, not a
within-protein trade-off. The gypsy result that was actually carrying the section
drops to p = 0.10. The one clade large enough to test internally (n = 32, 5
domains) gives p = 0.29.

**Independent check on REXdb.** The same rescoring was run on 13,572 plant Gags
(`scripts/nc_scan.py`, `scripts/nc_prld.py`). Knuckle calls fell the same way —
copia 70% by regex against 21% by profile — and the pattern is identical: a
strong pooled effect (copia OR 0.079, p = 4e-13) that vanishes within a lineage
(Ale OR 0.339, p = 0.14; Ivana OR 0.860, p = 0.78), and that reverses sign
between superfamilies (gypsy OR 3.04). Knuckle presence there is close to a
lineage property: 0% in Ogre, Retand and Reina, 56% in Tekay, 37% in Tork.

**What stands.** Nothing at the level the section originally claimed. There is a
pooled association in both databases; in neither does it survive controlling for
clade, and the two disagree on direction in gypsy. It is not refuted either —
within Ale only ~31 elements carry a profile-grade knuckle, so a twofold effect
would be undetectable. The status is *unsupported by any analysis that controls
for relatedness*, and the phylogenetic test is the only thing that can settle it.

#### Revision 2 (2026-08-13): all 336 cores through InterProScan — CURRENT

All 336 Gag cores were submitted to InterProScan in four batches
(`data/gag_plaac/interpro/part{1..4}.tsv`). Knuckle calls are counted only from
domain-level signatures of IPR001878 / IPR036875 — **PF00098, SM00343, PS50158,
SSF57756** — with PANTHER family names and PF19317 (the p24 capsid C-terminal
domain, whose InterPro title contains the word *nucleocapsid*) deliberately
excluded. 327 of 336 cores returned hits; only 3 carry any Pol signature.

**Why revision 1 was wrong.** Against InterPro, each cheaper definition scores:

| definition | calls | sensitivity | specificity |
| --- | --- | --- | --- |
| regex `C-X2-C-X4-H-X4-C` | 57% | 99% | 72% |
| Pfam clan, gathering threshold | 25% | **61%** | 100% |
| Pfam clan, E < 1e-3 | 36% | 87% | 99% |
| InterPro, domain-level | **40%** | — | — |

Pfam GA misses four knuckles in ten. Missing a feature roughly at random does not
push an odds ratio in a random direction — it drags it toward 1. Revision 1 read
that attenuation as absence of effect. 75 of the 135 InterPro calls are supported
by all four methods and 99 by at least three, so the calls it adds are not noise.

**The association on the best available calls** (250 elements of the 272 set, 25
domains):

| definition | knuckle in PrLD+ | in PrLD− | OR | p |
| --- | --- | --- | --- | --- |
| regex | 5/25 (20%) | 140/225 (62%) | 0.152 | 0.0001 |
| Pfam GA | 4/25 (16%) | 74/225 (33%) | 0.389 | 0.11 (n.s.) |
| Pfam E < 1e-3 | 4/25 (16%) | 103/225 (46%) | 0.226 | 0.0049 |
| **InterPro** | **4/25 (16%)** | **113/225 (50%)** | **0.189** | **0.0012** |
| InterPro + log(length) | | | **0.081** (0.023–0.287) | **0.0001** |

Length remains the dominant predictor (OR ~102 per e-fold), and adjusting for it
strengthens the knuckle effect rather than removing it — the original claim in
this section, now on calls that can be defended.

**It holds within superfamilies**, which is what revision 1 could not find:

| superfamily | knuckle in PrLD+ | in PrLD− | OR | p |
| --- | --- | --- | --- | --- |
| Ty3/Gypsy | 2/17 (12%) | 36/85 (42%) | 0.181 | 0.026 |
| Retroviridae | 2/5 (40%) | 38/44 (86%) | 0.105 | 0.037 |
| Ty1/Copia | 0/3 | 34/73 (47%) | — | 0.25 (3 events) |

**Disagreement with REXdb, unresolved.** REXdb plant copia agrees in direction
(pooled OR 0.079–0.080), but REXdb gypsy runs the other way (OR 2.5–3.0), and
that reversal was traced to lineage composition — knuckle presence there is close
to a lineage property (0% in Ogre, Retand and Reina; 56% in Tekay), so the
per-element odds ratio moves with which lineages are in the set. Whether GyDB's
gypsy result would survive the same lineage stratification is not known; its 17
domains sit in a handful of clades.

**What stands.** Gags carrying a prion-like domain are about a third as likely to
carry a zinc knuckle, at matched length, in a comparison that holds within two
superfamilies independently. The remaining exposure is phylogenetic
non-independence, untreated here as everywhere in this document, and it is not a
small caveat: 17 gypsy domains in a few clades could be a handful of events. The
finding should be stated as *supported but not yet controlled for relatedness*.

### 3.6 Disorder is widespread and graded

metapredict, fraction of the Gag core inside a predicted IDR:

| superfamily | n | median frac IDR | median N-term 40 | median longest IDR |
| --- | --- | --- | --- | --- |
| Ty3/Gypsy | 107 | 0.37 | 0.35 | 96 aa |
| Retroviridae | 49 | 0.31 | 0.25 | 104 aa |
| Ty1/Copia | 93 | 0.29 | 0.28 | 72 aa |
| **Bel/Pao** | 23 | **0.10** | 0.17 | 30 aa |

Overall median 0.30, IQR 0.17–0.42; 28 elements have no IDR at all.

Disorder and prion-likeness travel together, as they must if both tools are
working: elements with a called domain sit at 0.50 of the core disordered against
0.29 for the rest, and 0.48 vs 0.28 on the N-terminal 40 residues.

`Ty1B` and `Tkm1` are the extremes — N-terminal-40 disorder **0.92–0.93**, 0.55–0.59
of the core in an IDR, single IDRs of 174 and 189 aa.

### 3.7 Bel/Pao is a clean negative on both measures

23 elements, zero PrLDs, and three to four times less disordered than any other
superfamily. Two independent measures agreeing on the same clade.

### 3.8 Within Ty1/copia

Mapped onto the published 69-tip GyDB Ty1/copia RT tree (Branch 1 / Branch 2):
3 PrLDs, **all in Branch 1** — the fungal and diatom clade. Branch 2, which
contains every plant element (*Oryza*, *Vitis*, *Arabidopsis*, *Zea*,
*Nicotiana*, *Solanum*, *Glycine*, *Populus*), scores zero, with LLRs between −10
and −30. This reproduces decision F-10 in an independent database and on an
independent tree.

### 3.9 The apparent conservation of the PrLD is composition, not ancestry

**Added 2026-08-10.** §3.4 showed low pairwise identity among the four Pseudovirus
elements. This tests the same claim across all 24 α = 1.0 hits, and adds the null
control that §3.4 lacked.

All-vs-all over the 24 PrLD-bearing Gags, scored three ways — whole protein, PrLD
alone, and the Gag with the PrLD excised ("flanks"). **176 of 276 pairs cannot be
aligned over even 50% of the shorter protein** and are dropped; a global aligner
with free terminal gaps otherwise reports things like `ASSBSV`/`GypsySL_monotypic`
at "100% identity" over 0.3% coverage.

Of the 100 pairs that survive, the median whole-Gag identity is **20.0%** — the
twilight zone. Only four exceed 40% at near-full coverage, and all four are
relationships predictable without any analysis: `297`/`Tom` 58.2, `BFV`/`EFV` 45.5,
`HMS-Beagle`/`Yoyo` 42.0, `412`/`Mdg1` 40.6.

The PrLD looks *more* conserved than its own flanks — 27.8% against 19.8%, with 84
of 100 pairs favouring the PrLD. **That is entirely an artefact of composition.**
Against a shuffle that destroys order but preserves exact amino acid content:

| | real | shuffled | excess |
| --- | --- | --- | --- |
| PrLD pairs | 27.8 | 27.1 | **+0.7** |
| Flank pairs | 19.8 | 20.8 | −1.0 |

Two sequences that are both ~34% Q+N (§3.3) align well on that alone. Median excess
over the null is zero, so **there is no general sequence conservation of the domain
across the set.**

Six pairs clear 2 SD past their own shuffle distribution at ≥50% PrLD coverage:

| pair | Gag | flank | PrLD | PrLD cov | excess | z |
| --- | --- | --- | --- | --- | --- | --- |
| `297` / `Tom` | 58.2 | 64.6 | 41.6 | 95.3 | +13.7 | 5.5 |
| `EFV` / `FFV` | 36.7 | 37.9 | 36.4 | 78.8 | +6.8 | 3.2 |
| `BFV` / `FFV` | 36.1 | 36.7 | 37.6 | 91.7 | +7.2 | 2.2 |
| `17.6` / `Woot` | 19.3 | 18.6 | 31.8 | 89.1 | +6.2 | 2.6 |
| `Cyclops-2` / `Tom` | 18.8 | 28.6 | 36.6 | 88.7 | +7.2 | 2.2 |
| `297` / `Athila4-1` | 18.2 | 16.3 | 33.8 | 77.5 | +7.0 | 2.2 |

Only `297`/`Tom` is convincing. The other five sit at z = 2.2–3.2 on 20 shuffle
draws, which is marginal, and four more pairs were rejected for thin PrLD coverage
despite high z — `412`/`Burdock` scores +33.3 excess over 14.9% coverage, which is
noise. The coverage trap applies to the domain alignment, not only the protein one.

**Where real Gag homology exists, the PrLD is the least conserved part of it:**
`297`/`Tom` 41.6 against 64.6 in the flanks, `BFV`/`EFV` 33.1 against 52.1,
`HMS-Beagle`/`Yoyo` 25.0 against 43.1. The cleanest single case is
`HMS-Beagle`/`Yoyo` — 42% identical across the whole Gag, but their PrLDs score
**−2.5 against shuffle**, no detectable shared ancestry at all. The domain has
turned over completely while the protein around it stayed recognisably homologous.

This extends §3.4's "same location, different sequence" from Pseudovirus copia to
gypsy and Retroviridae, and strengthens it: not merely low identity, but low
identity relative to the same protein's own flanks, against a null.

**Method caveat.** Pairwise BLOSUM62 global alignment (gap −11/−1, free terminal
gaps), not the MAFFT L-INS-i used elsewhere in this document. It was run with
Biopython under the Windows `prioneers` env; the `copia-prld` env of §7 is intact
and holds MAFFT 7.525, so this is reproducible with L-INS-i and should be redone
that way if any of §3.9 is promoted to a decision. Adequate for a twilight-zone
yes/no; the identity values should not be quoted to a decimal. n = 24 and
phylogenetically non-independent, as everywhere else here.

---

## 4. Interpretations

Ordered by how much weight the evidence supports.

**Gag is special, and that supports the project's premise.** 26/26 in one protein
class with length controlled is a clean result. Prion-like character is not a
background rate across retroelement proteins.

**The Layer A / Layer B split is holding up.** Disorder is common, graded and
present in all four superfamilies; Q/N-richness is rare and phylogenetically
clumped. That is the shape of "conserved disordered substrate with compositional
character turning over on top of it". Consistent with H1 — not yet distinguished
from a few independent gains.

**A structural trade-off is a plausible mechanism — REINSTATED (2026-08-13).**
Gag appears to solve nucleic-acid binding either with a structured zinc knuckle
or with a prion-like low-complexity region, and less often with both. Withdrawn
earlier the same day on the Pfam rescoring; that withdrawal was wrong, because
the threshold used had 61% sensitivity and the resulting attenuation was mistaken
for absence. On InterPro calls (§3.5 revision 2) the effect is present at OR
0.081 adjusted for length, and holds within Ty3/Gypsy and Retroviridae
separately. This gives F-2's compensation idea a direction and a testable
prediction. What is still owed is the phylogenetic test — the domains cluster in
few clades, so the effective number of independent events may be small. REXdb
gypsy disagrees in sign for reasons traced to lineage composition.

**The copia N-terminal PrLD may be a yeast phenomenon rather than a copia one.**
Both N-terminal domains are in *Saccharomyces* and *Kluyveromyces*; the entire
plant radiation scores zero. But "Pseudovirus clade" and "fungal host" are the
same variable in this dataset, and PLAAC is a yeast-trained HMM scored against a
*S. cerevisiae* background. **These data cannot separate a biological signal from
a detector calibrated on the organisms it was built from.** This is the most
consequential open question here, because it decides whether the project's
central trait exists in plants at all.

**Gag PrLDs are probably not ancestral to LTR elements.** Bel/Pao has neither
disorder nor prion-likeness, which argues against a single inherited innovation.

---

## 5. Limitations

**Phylogenetic non-independence, untreated everywhere.** Every p-value above
treats 272 related elements as independent draws. The 17 gypsy PrLDs sit in a
handful of clades; a single ancestral gain would be counted 17 times. The CCHC
result, the Pseudovirus elevation and the per-superfamily rates are all exposed.
The fix is a phylogenetic logistic regression or a correlated-evolution test on
this tree. Per L-8 the threshold should be fixed **before** it runs, and that is
now compromised for anyone who has read this document.

**GyDB Gag cores are cut per element with no common rule.** Core length, CCHC
presence and disorder fraction all depend on where GyDB chose to start and stop.
This is exactly the failure mode that destroyed F-2's per-lineage CCHC figures
when they were rescanned on full-length Gag. **No absolute rate in this document
should be quoted** until it is checked against full-length sequence.

**α = 0.5 scores are batch-dependent.** They come from scoring all 2,636 cores
at once; scoring a subset alone changes calls (4 of 366 Gag entries flipped in a
direct test). At α = 1.0 the same elements are called bar two, so no conclusion
here rests on the choice — but the numbers are conditional on the input set.

**The alignment is thin.** 131 columns for 272 tips spanning four superfamilies.
Clade membership is reliable; deep branching order is not. Ten tips have RT cores
under 150 aa and are weakly placed.

**Small n where it matters most.** The Pseudovirus elevation rests on 7 elements
in one clade; the copia N-terminal result rests on 2.

---

## 6. Suggested next steps

1. **Re-derive full-length Gag** for these 272 elements and re-run PLAAC, the
   CCHC scan and metapredict on it. This retires the largest methodological
   caveat and re-tests findings 3.3, 3.5 and 3.6 simultaneously.
2. **Phylogenetic correlated-evolution test** for CCHC vs PrLD, and a
   phylogenetic signal test for PrLD presence, with thresholds fixed in advance
   by someone who has not seen §3.5.
3. **Resolve the detector/organism confound** — score a set of non-fungal,
   non-plant Gags of matched length, or calibrate PLAAC against a
   plant-background control, before treating the plant null as biological.
4. **Log the settled parts** as F-19 onward in `claude/DECISIONS.md`. Nothing
   here is citable under this project's own rules until that happens.

---

## 7. Files

All under `data/gag_plaac/`, which is gitignored.

| file | contents |
| --- | --- |
| `prep_272.py` | builds the set → `rt_272.faa`, `gag_272.faa`, `traits_272.tsv` |
| `disorder_272.py` | metapredict → `disorder_272.tsv` |
| `cchc_model.R` | length-stratified and logistic CCHC analysis |
| `figures_272.R` | daylight tree, LLR-by-protein-class violin |
| `figures_circular.R` | circular phylograms and cladograms, PrLD |
| `figures_disorder.R` | circular cladograms, disorder |
| `tree_272/` | `rt.aln`, `rt.trim`, `rt.treefile`, `rt.iqtree`, `rt.log` |
| `traits_272.tsv` | one row per element: superfamily, host, clade, LLR, PrLD call and position, CCHC |
| `all_cores_a05.tsv` | PLAAC over all 2,636 cores at α = 0.5 |
| `class_llr.tsv` | per-protein-class LLR table |
| `copia_tree_tips.txt` | the 69 Ty1/copia tips, transcribed by hand from the GyDB tree figure |
| `figures/` | all PNG and PDF output |

Added 2026-08-10, under `data/gydb/`, also gitignored:

| file | contents |
| --- | --- |
| `scripts/extract_prld_fasta.py` | cuts the 24 α = 1.0 hits out of the Gag cores |
| `tables/gydb_gag_prld.faa` | the called domains, one record per hit |
| `tables/gydb_gag_prld_parents.faa` | the parent Gag cores, same IDs and order |
| `tables/gydb_gag_prld_coords.tsv` | all four PLAAC spans per hit (PRD, core, LLR window, max window), plus midpoint and zone |
| `scripts/conservation.py` | §3.9 — pairwise identity with coverage gate and shuffle null |
| `tables/gydb_gag_prld_conservation.tsv` | one row per coverage-passing pair |

PLAAC reports four different spans per hit and they do not always agree: `Retrosat-2`
has its PRD at 285–357 but its best LLR window at 44–103. §3.2's positions use PRD.

`r-ggrepel` and `r-phangorn` were added to the `copia-prld` conda environment for
label repulsion and midpoint rooting.
