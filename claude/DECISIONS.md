# Decisions log

Every scientific parameter, with the reasoning behind it and its status.

Nothing is entered here until Cameron has signed off on it explicitly. A prior
version of this file recorded decisions from a build that was discarded; it was
cleared on 2026-08-04 so that every entry below is one he made himself.

**Status conventions**

- `CONFIRMED` — signed off explicitly, with the date. Results derived from it
  can be reported.
- `PROVISIONAL` — a value is in use so a stage can run, but no result derived
  from it may be reported until it is confirmed.

Anything not listed below is undecided. An undecided parameter is never given a
silent default: the stage that needs it refuses to run and names the missing
value.

---

## Confirmed decisions

### F-1 — Release scope: REXdb Viridiplantae v4.0 only — CONFIRMED (2026-08-04)

Metazoa v3.1 is not ingested. It would have added 177 Ty1/copia elements
carrying both an RT and a Gag slice (245 copia rows total; 197 with RT, 216 with
Gag) — about 3.3% of the combined copia set.

Three properties of the v3.1 release rule it out:

1. **No lineage assignment.** Column 5 of `Metazoa_v3.1.classification` is empty
   for all 245 copia rows. Viridiplantae v4.0 assigns 20 lineages (Ale 1,790;
   Ivana 852; SIRE 738; Tork 604; …). Lineage is the sampling frame for L-2 and
   the grouping variable for the comparative analysis, so these tips would be
   unassignable.
2. **No element DNA locally**, so metazoan Gag slices could not be extended to
   the ORF start if F-4 goes the DNA route. Their trait would be scored over the
   capsid core while plant tips were scored over core-plus-N-terminus — a
   methodological difference aligned exactly with the deepest split in the tree,
   which would manufacture a plant/animal trait difference by construction.
3. **Separate accession namespace** (NCBI `NC_######.#` vs `REXdb_ID####`).

Not established: whether a Metazoa species table and DNA file exist upstream and
were merely not downloaded. It does not matter for this decision — the empty
lineage column is intrinsic to the classification file and blocks inclusion on
its own. Revisit only if metazoan tips become scientifically necessary, and then
as a deliberate reopening.

Retained: 18,126 elements in the Viridiplantae v4.0 release, of which 5,606 are
Ty1/copia and 5,164 of those carry both an RT and a Gag slice. These counts are
before any ambiguity or length filtering (F-6, F-8).

### F-2 — Ty3/gypsy retained, written to a separate file — CONFIRMED (2026-08-04)

Both superfamilies are fetched. Copia (5,606 elements) is the study set; gypsy
(8,749) is written to its own file and **never merged** into it.

**Justification is functional contrast, not rooting.** The original rationale
was that Ty1/copia Gag largely lacks canonical CCHC zinc-knuckle nucleocapsid
motifs while Ty3/gypsy carries them, so gypsy would serve as a
structured-NC control.

**That claim was wrong and is retracted (2026-08-04).** Measured on the
extracted Gag sequences with the pattern `C-X2-C-X4-H-X4-C`: **copia 51.1%**
(2,464/4,821) carry at least one, **gypsy 23.9%** (1,909/7,985). Copia carries
them *more* often, not less. The decision to retain gypsy stands, but this
particular argument no longer supports it — see the within-copia contrast below,
which is a better test anyway.

Every hit in both superfamilies falls in the **downstream extension**, none in
the published core or the upstream extension. REXdb's Gag slice therefore omits
the nucleocapsid region entirely — independent validation that the +120 aa
offset is recovering real functional Gag sequence.

**The better contrast is within copia**, where CCHC presence varies sharply by
lineage: Angela 83.2%, TAR 81.8%, Tork 74.5%, Ikeros 70.0%, Ale 50.3%, Ivana
40.0%, SIRE 11.2%, Bianca 1.8%. Testing whether lineages lacking a structured NC
compensate with more disordered or prion-like N-terminal character is a
same-ruler, same-profile comparison with the phylogeny available to control for
relatedness — free of every cross-superfamily artifact that troubles the
copia/gypsy comparison.

**Tested, and the lineage differences were an artifact — retracted.** Rescanning
the full Gag→PROT region instead of a fixed +120 aa window: Bianca 3.3% → 84.6%,
SIRE 25.0% → 62.5%, Ale 44.8% → 74.7%, while Angela (88.2%) and Tork (78.4%)
were unchanged. The apparent CCHC deficit in the long-Gag lineages was entirely
a consequence of the extraction window, not biology. **Do not cite the
per-lineage CCHC figures.** This is what prompted dropping the downstream
extension in F-4.

The copia 51.1% vs gypsy 23.9% comparison was measured under the same window
limitation and is therefore also **provisional** — it would need re-running on
full-length Gag before it could be trusted. The current extraction is
upstream-only, so it cannot support any CCHC claim at all; a separate
full-length extraction would be required.

Secondary: gypsy anchors the work to the domesticated-Gag condensate literature
(Arc / dArc, PEG10, PNMA), which is a gypsy-lineage literature, not a copia one.

**Rooting is explicitly NOT the justification, and is deferred to L-3.** REXdb
cuts the two superfamilies with different HMMs — every copia RT slice comes from
a `Ty1-RT` profile (median 256 aa, IQR 256–258), every gypsy RT slice from a
`Ty3-RT` profile (median 173 aa, IQR 173–174). The ~83 aa difference is a
profile-boundary artifact, not biology; RT is a conserved fold and gypsy RT is
not a third shorter than copia RT. Aligning them jointly would put an ~83-column
gap block in every gypsy tip, so either trimAl strips those columns and costs a
third of the copia RT signal, or they are retained and all gypsy tips share an
identical missing-data pattern that drives root placement.

**The same ruler problem applies to Gag** — copia Gag from `Ty1-GAG` (median
91 aa), gypsy from `Ty3-GAG` (median 125 aa) — so a naive cross-superfamily
trait comparison would partly measure profile boundaries.

Two things keep this decision safe:

1. The trait contrast does **not** require a joint alignment. Gag can be scored
   per superfamily and the distributions compared without ever aligning gypsy
   to copia.
2. If F-4 re-derives slices from element DNA and extends Gag to the ORF start,
   both superfamilies become defined by the same rule rather than by two
   profiles, which dissolves the Gag ruler problem. F-2 is therefore provisional
   in its *usefulness* even though the fetch decision itself is settled.

Separate files are the mechanism that stops the merged-set artifact from
happening by accident.

### F-4 / F-5 — Gag is extracted from element DNA by fixed in-frame offsets — CONFIRMED (2026-08-04)

REXdb's published Gag slice is the capsid core only (copia median 91 aa) and
does not contain the N-terminal region the project is about. Gag is therefore
re-derived from `Viridiplantae_v4.0_ALL_DNA.fasta`.

**The rule.** Locate the published GAG slice in a six-frame translation of the
element's own DNA. That fixes the reading frame. Then extend **in that frame**:

    75 aa upstream   (clipped if an in-frame stop arrives first)
    [ GAG slice = CA core, ~91 aa ]

Nothing is taken downstream. Nominal product ~166 aa. F-5 is subsumed: the
boundary *is* the offset rule.

**Upstream only.** The hypothesis concerns the disordered N-terminal region
upstream of the structured capsid core, so the C-terminal half of Gag is not
what is being measured. A 120 aa downstream offset was implemented first and
then removed, because it introduced a **lineage-correlated coverage bias**:
lineages with a longer Gag had C-terminal features fall outside the window.
Measured on the CCHC nucleocapsid motif, which sits near Gag's C-terminus —
detection within +120 aa vs within the full Gag→PROT region was Angela
88.2%/88.2%, Tork 78.4%/78.4%, Ale 44.8%/74.7%, SIRE 25.0%/62.5%, **Bianca
3.3%/84.6%**. Median motif position tracks the Gag→PROT gap (Tork 79/135,
Angela 87/137, Ale 98/171, SIRE 113/193, Bianca 122/195), so a fixed downstream
cut necessarily truncates long-Gag lineages preferentially.

Upstream has no equivalent problem because it is bounded by a **real biological
boundary — the in-frame stop** — rather than by an arbitrary offset landing
mid-domain.

Dropping downstream improved every uniformity measure. Copia length CV 0.110
(was 0.139 with 75/120; true-ORF extraction was 0.28); gypsy CV 0.080 (was
0.125). Within ±10% of the median: copia 85.2%, gypsy 88.8% (both 75.9% before).
Elements receiving the full offset: copia 76.7%, gypsy 83.0% (was 62.9%/64.0%
for both offsets). Per-lineage medians span 4 aa — Ale 167, Ivana 168, SIRE 168,
Tork 165, Angela 165, Ikeros 165, Bianca 164, TAR 165 — so the coverage bias is
gone, not merely reduced.

It is also simpler: no PROT lookup, no cross-frame logic, half the locate calls.

**Structural justification.** Retroviral Gag is MA–CA–NC. The CCHC nucleocapsid
motif was found exclusively *downstream* of REXdb's published slice — never
inside it, never upstream — which identifies the published slice as
approximately the **CA core**, and everything upstream of it as the
MA-equivalent N-terminal region. That is precisely the target region, confirmed
empirically rather than assumed.

**Known limitation, carried forward.** At ~75 aa the recovered N-terminal region
is barely larger than PLAAC's 60 aa default core window, so a window-based prion
detector is being asked to characterise ~1.25 windows of sequence. This makes
the PLAAC core-length parameter (L-5) unusually consequential and strengthens
the case for the Layer A / Layer B split, where disorder and charge do work
PLAAC cannot. Independent of this decision — it holds for any extraction rule.

`up_got` doubles as the index where the CA core begins in the emitted sequence,
so the N-terminal region and the core can be scored separately.

**Result:** copia 4,821 of 5,164 located (93.4%), median length 165 aa, IQR
164–168. Gypsy 7,985 of 8,582 (93.0%), median 198 aa, IQR 195–201. No extracted
sequence contains a stop codon. 675 copia and 1,993 gypsy carry an ambiguous
residue — flagged, not filtered (F-6).

**Why fixed offsets rather than true ORF boundaries.** Cameron's reasoning, and
it is correct: the reading frame comes from the GAG slice itself, so no start
codon has to be inferred at all. Getting the true start exactly right does not
matter for a compositional trait provided the same rule is applied to every
element — the error becomes a systematic offset that largely cancels in a
comparative analysis rather than random per-element noise.

It also removes a confound found while testing an earlier proposal: true-ORF
extraction gives Gag lengths varying 1.8× between lineages (SIRE median 564 aa
vs Tork 308; pooled CV 0.28, only 71.8% within ±20% of the median). PLAAC's LLR
scales with the length scanned, so variable-length Gag would hand longer-Gag
lineages more opportunity to score, and that would map onto the tree as false
phylogenetic signal. Uniform extraction defuses this at the source.

**Why the offsets are 75 and 120 rather than larger.** Measured on 2,270 copia
elements. Clean upstream reading before an in-frame stop: median 88 aa, IQR
74–101. Fraction of elements reaching a given upstream offset — 60 aa: 83.6%;
**75 aa: 74.4%**; 87 aa: 55.9%; 120 aa: 9.3%; 170 aa: **3.2%**.

An upstream offset of 170 was proposed and rejected on this evidence. Past the
in-frame stop the extension is translating 5′UTR/LTR as though it were protein —
a different sequence class, landing in exactly the N-terminal window the PrLD
analysis examines, and the most likely source of a spurious Q/N signal. The
product would also contain literal `*` characters that PLAAC cannot score.

Downstream, from the GAG slice end to the PROT slice start: median 157 aa, IQR
137–199. Fraction staying inside Gag — 85 aa: 99.9%; **120 aa: 98.5%**;
157 aa: 50.2%. 120 captures ~76% of the median Gag C-terminal region while
almost never reaching protease.

**Clipping is recorded, not silent.** Every element carries the requested and
achieved offsets and what stopped the extension, so the uniformity assumption is
auditable rather than assumed.

**Supporting structure, verified in the data.** Domain order is diagnostic and
consistent — copia `GAG-PROT-INT-RT-RH` (163/201 sampled), gypsy
`GAG-PROT-RT-RH-INT` (164/206), the textbook INT-position distinction. Domain
lengths are tight (copia PROT 69–72, RT 256–258, INT 195–200). Inter-domain
spacing is not. PROT is present on all 5,606 copia elements, so the downstream
anchor exists for every element that has a Gag.

**Rejected alternatives.** A fixed 450 aa window from a 5′ anchor: 85.7% of
copia Gags are shorter than 450, so it would read into Pol, and there is no
usable 5′ anchor — LTR terminal repeats were detectable in only 3.8% of records,
and the GAG slice sits 383–6,559 nt into the record (median 901).

**Probe findings that stand behind this.** Gag slices are locatable in their own
element DNA for 93.4% of copia and 93.0% of gypsy — not the 99.3% claimed by a
300-element sample in the discarded build. All hits are on forward frames, so
REXdb DNA records are already oriented. 7.6% of located copia elements have no
ATG anywhere in the open upstream region, which is why an ATG-based rule was
avoided independently of the arguments above.

**Deferred:** whether extracted regions containing ambiguous residues are
filtered (F-6), how RT-only elements are handled (F-7), length filtering (F-8),
and the output schema (F-9). Extraction records these as flags and filters
nothing.

### Framing note — the hypothesis splits into two layers (2026-08-04)

Arising from the F-2 discussion, recorded because it changes what the detectors
are for. Not itself a fetch decision.

If the Gag N-terminus is conserved because it does **RNA chaperone work**, what
is conserved is *disorder plus net positive charge* (K/R-rich), not Q/N
richness. PLAAC is a yeast prion HMM trained on Q/N-rich domains, so a
conserved, functionally essential, basic disordered region can score as **not
prion-like**. A low PLAAC LLR across the superfamily would then not be evidence
against the hypothesis — it would be the detector being off-target.

So the hypothesis is better stated in two layers:

- **Layer A — functional, expected conserved:** disorder + basicity. Measured
  with IUPred3/metapredict and charge composition, *not* PLAAC.
- **Layer B — compositional, expected to turn over:** Q/N character
  specifically. Measured with PLAAC.

H1 becomes: *Layer A is conserved; Layer B turns over on top of it.* This makes
the confound check (PrP, Luminidependens, FUS-family) load-bearing rather than
optional, since it calibrates how much of Layer B a yeast-trained predictor
misses outside fungi.

Precedent: yeast prion domains, Sup35 especially, are documented to evolve
rapidly in primary sequence while retaining compositional bias — compositional
turnover under conserved character, in another system.

**Not yet verified against the literature** (flagged, load-bearing): the exact
internal start position of the Ty1 p22/p18 restriction factor and which residues
it lacks; whether a direct N-terminal deletion series shows loss of
transposition; and the specific P-body dependencies of Ty1 retrosome formation.

### F-3 — No separate GyDB fetcher; REXdb is the only source — CONFIRMED (2026-08-04)

GyDB is already integrated into REXdb. 98 elements in the Viridiplantae v4.0
release carry `Data.source == "http://gydb.org/index.php/Main_Page"` (48 copia,
50 gypsy); the bulk of the release is phytozome (14,483) and genbank (2,636).

Decisive beyond the overlap: a GyDB-only sequence would come without element
DNA, so the F-4 rule could not recover its Gag N-terminus. Its Gag would be
capsid-core-only while every other tip carried the N-terminal region — the same
non-comparability that ruled out Metazoa in F-1.

Sources are therefore REXdb Viridiplantae v4.0, and nothing else.

### F-10 — Reference sequences extracted at full length — CONFIRMED (2026-08-04)

The five *S. cerevisiae* Ty1/copia elements with both RT and Gag — `GEN_Ty1B`
(5405), `GEN_Ty2` (5406), `YARCTy1-1` (5412), `YBL100W-B` (5413), `YLRWTy1-3`
(5415) — are written **twice**: to `gag_reference.faa` at their full open
upstream length, and to the study set under the uniform 75 aa rule so Ty1
remains a tree tip. (`GEN_Ty4` 5407 and `YHLWTy4-1` 5414 have RT but no Gag
slice.)

**Why full length is needed for these.** Ty1 is the canonical positive control,
and the 75 aa rule truncates its prion-like domain. Its open upstream region is
**223 aa**, beginning `MESQQLS…` — the canonical Ty1 Gag start — with the
Q/N/S/T/P-rich stretch (`QQLSQ…SQQPTTPPSSA…PQQ…QQAN…QSQFTQYPQY`) sitting
*distally*, exactly where the 75 aa window cuts. Q+N is 16.1% across the full
region but only 13.3% in the 75 aa retained, and the retained part is the bland
linker adjacent to the core. All five agree: 215–227 aa open upstream, Q+N
13.5–16.3% full versus 10.7–14.7% windowed.

**The 75 aa rule was NOT changed for the study set, and this is the important
finding.** Plant copia is genuinely short and genuinely not Q/N-rich — median
open upstream 89 aa, only 5.6% reach 150 aa and 1.6% reach Ty1's 223. Mean Q+N
is flat across every length bucket (≤100 aa: 8.8%; 101–150: 8.9%; 151–200: 8.0%;
>200: 8.5%), and among the 155 elements with ≥150 aa upstream the distal 75 aa
we discard is 8.9% Q+N against 7.9% proximal, higher in only 38% of cases. So no
prion-like domain is being truncated in plants; there is none there to truncate.
~8.8% is close to background protein composition.

**Consequence to carry forward:** the Q/N-rich N-terminal domain is real in yeast
Ty1 and absent in plant copia at the compositional level. A near-null PLAAC
result on the study set should be expected, and would reflect a detector/organism
mismatch rather than absence of an N-terminal trait. This is direct support for
the Layer A / Layer B split. Cameron is running PLAAC, charge, disorder and LLPS
himself; this decision covers only how the regions are pulled.

### F-11 — Exploratory scripts deleted — CONFIRMED (2026-08-04)

`scripts/probe_f4_gag_nterm.py` and `scripts/extract_gag.py` are deleted, along
with their outputs (`probe_f4_gag_nterm.tsv`, `gag_coords.tsv`, `gag_copia.faa`,
`gag_gypsy.faa`).

Both were scaffolding for reaching F-4 and F-10, not pipeline stages. Every
number they produced is recorded in this file, so the scripts carry no
information that would be lost — and keeping a superseded extractor beside the
real one invites the wrong file being run. `extract_gag.py` is absorbed into
`scripts/fetch.py`.

The repo keeps one script per stage. Analyses that only exist to settle a
decision are deleted once the decision is logged.

### F-12 — `fetch` stays standard-library only; no Biopython — CONFIRMED (2026-08-04)

Biopython would replace the codon table, `translate`, reverse-complement and the
FASTA reader — about 20 lines of 170. Rejected because:

- the script is otherwise dependency-free, which matters for reproducibility;
- `SeqIO.parse` builds a `SeqRecord` per entry, slower over the 134 MB DNA file
  than a plain line reader;
- `Seq.translate` warns on partial codons, so the sequence needs trimming to a
  multiple of three anyway.

The argument *for* it was that a hand-rolled genetic code can hide a silent bug.
That was tested, not assumed: the table has 64 codons, stops are exactly
`TAA/TAG/TGA`, and it matches 23 independently known codon assignments with no
mismatches.

No library covers the part that is actually complex — locating a published
protein slice with `X` wildcards across six reading frames.

Revisit if a later stage needs Biopython for its own reasons; there is no point
carrying the dependency for `fetch` alone.

**Also fixed in the same pass:** the Gag output was sorted with
`order.index(...)` inside the sort key, a linear scan per item (~115 M
operations for gypsy). Now sorted on the numeric element ID. Output is
byte-identical apart from `gag_reference.faa`, which is now in numeric rather
than file order.

### F-13 — Reverse reading frames removed — CONFIRMED (2026-08-05)

`fetch` searches only the three forward frames. All 12,806 located elements
matched forward (+1: 4,362, +2: 4,376, +3: 4,068); zero matched on a reverse
frame. REXdb DNA records are already oriented.

Provably output-neutral: the 940 elements that fail to locate were already being
searched against all six frames and failed, so no reverse-findable element
exists in this release. Output is byte-identical after the change.

`translate` was renamed `translate_codons` and `reading_frames` to
`forward_frames` — the old name collided with Python's `str.translate`, which
the reverse-complement line used to call three lines away.

**Assumption this introduces:** a future release with differently oriented
records would silently locate fewer elements. The `gag_located` count in
`elements.tsv` is the check.

### F-14 — 75 aa kept; the ORF start is not recorded — CONFIRMED (2026-08-05)

The extraction rule is unchanged. A `gag_orf_offset` column was added and then
**removed**: it altered no sequence, and re-running `fetch` regenerates it in
~25 s, so carrying it was not worth a column. To restore it, record
`len(open_region) - open_region.find("M")` where `open_region` is
`extend_upstream(protein, index, None)`.

The measurements below are the reason the rule was left alone, and are worth
keeping even though the column is gone.

**Why the question was reopened.** The in-frame stop is not the 5′UTR boundary.
It is the first *random* stop somewhere inside the UTR — it marks where the ORF
cannot extend past, not where it begins. The start codon is the real boundary.

**What the 75 aa rule actually does**, measured on 4,821 located copia elements:
2,540 (52.7%) have their ORF start further back than what was kept, so the
window is entirely coding and real Gag is truncated; 1,914 (39.7%) have it
inside, so some untranslated sequence is included (median 11 aa, IQR 4–25, p90
56, max 74); 367 (7.6%) have no ATG at all.

**Why cutting at the ATG was rejected.** It doubles length variation — CV 0.216
against 0.110, and elements within ±10% of the median fall from 85.2% to 55.1%.
Worse, the variation is lineage-correlated: Bianca 149, Tork 154, Ivana 166,
Angela 169, SIRE 180, Ale 180, TAR 182 — a 22% spread, against 2% for the fixed
window (164–168 across every lineage). That is the same length confound F-4 was
built to remove, reintroduced.

So both rules carry a lineage-correlated error. Recording the offset keeps the
sequence files uniform for scoring while leaving the ORF-start version one
subtraction away, with no need to re-run the DNA pass. Which is right is then an
empirical question — score both and see whether it moves anything — rather than
a decision made on principle now.

Reference sequences are unaffected: they are already uncapped (F-10). `GEN_Ty1B`
has its ORF start 219 aa before the core against the 75 kept, so 144 aa of
coding is truncated in its study-set representation — which is why the
full-length reference file exists.

### F-6 — No ambiguity filter; counts are recorded instead — CONFIRMED (2026-08-05)

Sequences are never excluded for unresolved residues. `rt_ambiguous` and
`gag_ambiguous` carry the counts so scoring can exclude or weight.

**Reached in three passes, each time because the ceiling was measured rather
than assumed.** The first decision was RT ≤ 2 / Gag = 0, argued from the
principle that a single `X` silently shifts a composition fraction. True, but
the magnitude was never checked. It is trivial:

- Of 675 copia Gag sequences carrying any `X`, **573 (85%) have exactly one**;
  median density 0.6% of the sequence.
- **The worst sequence in the entire release** has 13 `X` in 256 aa of copia RT
  (5.1%) and 6 in 164 aa of copia Gag (3.7%).

For scale, IQ-TREE routinely handles alignments with 20–30% gaps, and 6 unknown
residues out of 164 caps a Q+N bias at 3.7 pp — about 1.3 SD, and only for a
single element. Filtering at ≤ 2 excluded 247 RT and 210 Gag sequences to avoid
that.

Since ambiguity plausibly tracks element decay, any threshold removes old
elements preferentially — the lineage-biased loss this project keeps trying to
avoid — for a bounded and small gain. **That correlation was never measured**; it
is the reason the filter was dropped, not an established result.

Incidentally, the F-7 linkage removes the extreme RT cases anyway: the 13-`X`
element has a failing Gag, so the highest ambiguity reaching the tree is 7.

**Where the ambiguity comes from.** Almost entirely REXdb, not our translation.
The published protein FASTA carries `x` (18,769) and `X` (6,302) marking codons
REXdb could not resolve; `fetch` uppercases on read, so the lowercase majority
is caught. Our own `translate_codons` adds an `X` for any codon containing a
non-ACGT base, and the DNA release does hold IUPAC codes (Y 3,430, R 3,171,
M 1,400, K 851, W 700, S 242, N 29) — but these rarely land in the 75 aa window.

Splitting copia Gag by which part the `X` sits in: 663 in the published core
(REXdb), 8 in the upstream extension (ours), 4 in both, 4,146 clean. RT is used
as published, so all RT ambiguity is REXdb's.

**Known gap, accepted.** The filter counts only `X`, but the protein FASTA also
uses `B` (Asx) and `J` (Xle). In the copia/gypsy slices actually used that is
2 residues, both in RT, none in Gag — inside RT's ≤2 allowance already. Would
need fixing if either ever appeared in Gag, where tolerance is zero.

### F-7 — Missing traits are labelled, not deleted — CONFIRMED (2026-08-05)

**Every element with an RT enters the tree** — 5,606 copia, 8,749 gypsy. Where
the Gag is missing, `gag_status` records why:

| status | copia | gypsy | meaning |
| --- | --- | --- | --- |
| `ok` | 5,042 | 8,486 | extracted and written |
| `no_slice` | 442 | 167 | REXdb published no Gag — genuine absence |
| `unplaceable` | 89 | 85 | slice exists, not findable in its own DNA |
| `no_nterm` | 33 | 11 | placed, nothing upstream of the core (F-8) |

**Revised from an earlier form that deleted instead.** The original concern was
Cameron's and it was right: a tip with no trait must not silently mix genuine
absence with our own quality control. The first fix was to drop such elements
from the tree entirely, which cost 375 copia and 608 gypsy tips whose RT was
perfectly good — and since the failures track element decay, that loss was
age- and lineage-biased in exactly the way F-6 was relaxed to avoid.

Labelling solves the conflation properly. Deleting removes the ambiguity by
removing the data; a status column removes it by recording the distinction, and
downstream can then decide whether a pseudogenized Gag belongs in its analysis.

Gag still requires a usable RT — an element absent from the tree has nowhere to
carry a trait.

### F-15 — Gag located by anchoring on the core's first 30 aa — CONFIRMED (2026-08-05)

`find_slice` tries the whole published slice, then falls back to its first
30/25/20 aa. `gag_core_exact` records which: 4,789 copia matched in full, 253
were head-anchored (7,974 and 512 for gypsy).

**Why the full slice often fails — this was diagnosed, not assumed.** Of the 343
copia that previously failed: 328 (96%) had *some* 20-aa chunks present in the
DNA but not all, and mapping chunks to frames showed the signature directly —
e.g. `REXdb_ID8` reads `+3 +3 +3 none none +2 +2 +2`. The Gag starts in one
frame and continues in another. These are **frameshifted, pseudogenized Gag
ORFs**, with REXdb's annotation stitching a protein across the break. 146 showed
an outright frame change, 187 a mismatch without one, 10 had nothing matching.
Ambiguous bases were not the cause — 106 failures had them, but wildcard
matching rescued zero.

**Why the head is enough.** The core is only an anchor for walking upstream, and
the break is almost always downstream of the anchor. Head-anchoring recovered
245 of 343 copia with a median 75 aa upstream — the full window, identical in
quality to any other element.

**Specificity was checked before adopting it.** A 30 aa head anchor placed
uniquely on all 4,580 control elements and all 201 rescued ones, with **zero**
multiple matches anywhere.

Result: `unplaceable` fell from 343 to 89 (copia) and 597 to 85 (gypsy).

**The caveat `gag_core_exact` exists for.** In a head-anchored element the core
diverges from its DNA after ~30 aa, so that half of the emitted sequence is
REXdb's reconstruction of a broken region rather than a translation of this
element's DNA. The upstream half is cleanly derived either way. Scoring can use
the flag to exclude the core portion where it matters.

### F-8 — Drop only elements with a 0 aa N-terminal region — CONFIRMED (2026-08-05)

22 copia and 7 gypsy, after the F-6 filters. These have the capsid core and
nothing upstream, so scoring them against elements with 75 aa of N-terminus
measures window length rather than biology.

No larger threshold. ≥25 aa would drop 304 copia and ≥50 aa would drop 571, but
any such cutoff is invented, and clipping tracks decay — so the loss would be
lineage-biased in exactly the way F-6 avoids.

**Where filtering happens.** In `fetch`, not a separate stage. The FASTA files
carry only passing sequences; `elements.tsv` retains all 14,355 rows with
`in_tree` and `gag_status`, so every drop and its reason stay on record.
Reference sequences (F-10) bypass the filters — they are controls.

---

### F-16 — Elements whose two release files disagree are dropped — CONFIRMED (2026-08-05)

`gag_status == "record_mismatch"`: no 20 aa stretch of the published Gag slice
appears anywhere in that element's DNA. These are removed from the RT set as
well — the only elements dropped from the tree.

**They are all *Malus domestica*** — 10 copia (8 Ale, 2 Bianca) and 1 gypsy, all
phytozome-sourced, RT lengths normal at 252–262 aa. A systematic protein/DNA
record mismatch for one species' batch in REXdb, not a property of the elements.

**Distinct from the frameshift cases (F-15).** A frameshifted Gag still has most
of itself in the DNA, just not contiguously — 79 copia and 84 gypsy remain
`unplaceable` on that basis and **keep their place in the tree**, since their RT
is sound and their Gag absence is biology. Nothing matching at all means the two
files describe different sequences, which makes the record itself untrustworthy.

**"No trait" was not the criterion.** 564 copia tips have no trait; only these
10 leave. The tree is built from every available RT and pruned to trait-complete
tips for the comparative analysis, which estimates topology and branch lengths
from more taxa and lets the pruned tree inherit that. `gag_status` makes the
pruning a one-line filter.

Tree: copia 5,606 → 5,596, gypsy 8,749 → 8,748.

### F-17 — `search_frames` fast path removed; `rt_pass` renamed `in_tree` — CONFIRMED (2026-08-05)

`search_frames` had a separate exact-match branch for queries without an `X`.
It was redundant: a query with no `X` splits to itself, so the anchor *is* the
whole query and the wildcard comparison reduces to an exact one. Verified on all
13,746 Gag slices — identical frame and index in every case. Removing it drops
7 lines and one code path.

`rt_pass` became `in_tree`. It no longer means "passed a quality filter" — since
F-6 dropped the ambiguity filter, the only exclusion is `record_mismatch`, so
the column simply records whether the element's RT was written to `rt_*.faa`.

The other six columns that duplicate facts recoverable from the FASTA files
(`rt_len`, `rt_ambiguous`, `gag_len`, `gag_ambiguous`, `gag_core_len`, and
`in_tree` itself) were **kept deliberately**: `elements.tsv` is the audit record
for this stage and should be filterable without loading sequence files.

### F-9 — `REXdb_ID` is the join key — CONFIRMED (2026-08-05)

Settled by implementation rather than discussion. Every FASTA header is the bare
`REXdb_ID`, and `elements.tsv` has one row per element keyed on the same, so the
tree and the trait join on it directly. All 18,126 classification rows have a
matching info row and protein entry, so the key is complete with no fallback
needed. `elements.tsv` carries 20 columns and every copia/gypsy element,
including those contributing no sequence to either FASTA.

---

## `fetch` is complete

F-1 to F-15 are all decided, implemented and reflected in `scripts/fetch.py`.
Nothing in the fetch stage is open.

Two things recorded but never established, both flagged where they matter:

- **Ambiguity is assumed to track element decay.** It is the reason F-6 dropped
  its filter and F-7 stopped deleting, but the correlation was never measured.
- **10 copia elements have a Gag slice with no part present in their own DNA.**
  Distinct from the frameshift cases diagnosed in F-15, and not investigated.
  They fall under `unplaceable`.

## Decision queue — later stages

Not reached yet, recorded so they are not forgotten.

| # | decision | stage |
| --- | --- | --- |
| L-1 | MMseqs2 identity and coverage thresholds | cluster |
| L-2 | Sampling rule — elements per lineage per host taxon | sampling |
| L-3 | Whether the RT tree is actually rooted on gypsy | tree |
| L-4 | PLAAC α (0 / 0.5 / 1.0), and whether all three run as a sensitivity analysis | score |
| L-5 | PLAAC core length | score |
| L-6 | Disorder tool: IUPred3 or metapredict | score |
| L-7 | catGRANULE version: original or 2.0 / ROBOT | score |
| L-8 | Phylogenetic signal threshold — must be fixed before the tree is examined | signal |
| L-9 | Transition-count threshold — same pre-registration logic | ancestral |
