# Build Brief — Ty1/copia Gag PrLD Phylogenetics Pipeline

You are helping me scaffold a **new, clean, reproducible** GitHub repository from scratch. Read this entire brief before writing any code. **Do not start generating files until you have asked me the clarifying questions in Section 9 and I have answered them.** This is a scientific pipeline for a PhD project; correctness and reproducibility matter more than speed.

---

## 1. What this repository does (the science, in one paragraph)

This pipeline tests whether a disordered N-terminal region in the Gag protein of Ty1/copia LTR retrotransposons is **conserved across the superfamily**, while its **prion-like composition varies independently** (compositional turnover) — as opposed to the prion-like domain (PrLD) being a single inherited innovation (acquisition by descent). It does this by: (1) pulling Ty1/copia element sequences from curated databases, (2) building a phylogenetic tree from the **reverse transcriptase (RT)** domain, which is alignable superfamily-wide, and (3) scoring the **Gag N-terminal region** for prion-like / disordered / condensation-prone character, then (4) mapping that Gag trait onto the RT tree as tip data for comparative analysis.

**The single most important design principle:** the tree and the trait come from **two different regions of the same element**. RT is *aligned* to build the tree. Gag is *scored* to make the trait. They are joined by element accession. Never build the tree from Gag (it is too divergent to align across the superfamily), and never treat a per-sequence score as phylogenetic input.

---

## 2. Hard requirements for the repository

- **Reproducible by a stranger.** Someone who clones the repo and follows the README must be able to reproduce the analysis end-to-end with pinned software versions and no manual fiddling beyond documented config.
- **Clean and concise.** No dead code, no exploratory notebooks committed, no hardcoded absolute paths, no machine-specific paths (I have Windows-style `Rscript` paths in my *old* pipeline — do NOT carry those over; detect or configure instead).
- **One command per stage, plus one command to run everything.** A user should be able to run the whole pipeline with a single entry point, or run any stage independently.
- **Deterministic where possible.** Fixed random seeds; pinned tool versions; config-driven parameters (no magic numbers buried in functions).
- **Config-first.** All scientific parameters (see Section 6) live in a single human-readable config file, not scattered through the code. Changing a parameter must never require editing source.
- **Auditable.** Every stage logs what it did, writes intermediate outputs to disk, and can be re-run from cached intermediates without re-hitting the network.

---

## 3. Languages and tooling

- **Pipeline orchestration language:** Python (I am an intermediate Python user — favour clear, well-commented, standard-library-leaning code over clever abstractions).
- **Downstream statistical / phylogenetic comparative analysis:** R (I am more proficient in R; the phylogenetic signal and ancestral-reconstruction steps should be R scripts called from the pipeline).
- **External bioinformatics tools** the pipeline will orchestrate (all must be pinned in the environment file, and the README must explain how to install them):
  - **HMMER** (`hmmsearch`) — domain-gating fallback only: locate RT and Gag domains within each protein by scanning against **two targeted retroelement-specific profiles** (GyDB/REXdb-derived), never against the full Pfam database. Used only where the database does not already provide domain coordinates.
  - **MAFFT** — multiple sequence alignment of the RT region.
  - **trimAl** — alignment trimming.
  - **IQ-TREE** — maximum-likelihood phylogenetic inference from the RT alignment.
  - **PLAAC** — prion-like amino acid composition scoring of the Gag region (Java jar; called via subprocess).
  - **MMseqs2** — sequence-identity clustering to remove copy-number pseudoreplication.
  - Trait detectors beyond PLAAC (see Section 6 — some are web tools or separate installs; the pipeline should be structured so these plug in cleanly even if I run one or two of them manually at first).
- **No UniProt/NCBI client and no genome-mining tools** (LTR_retriever/LTRharvest) — sequences come from GyDB/REXdb only (Section 4), so do not add those dependencies.
- **Environment management:** provide a pinned `environment.yml` (conda) as the primary path, since most of these tools install cleanly from bioconda. Also provide a `requirements.txt` for the pure-Python deps. Pin versions.

---

## 4. Data sources (where sequences come from)

**Sequences are pulled from GyDB and REXdb ONLY.** Do not add any other sequence source — no UniProt-by-Pfam pull, no NCBI/Entrez, no genome mining (LTR_retriever/LTRharvest). If you think another source is needed to make the analysis viable, stop and ask me rather than adding it. The fetching layer must be modular per-source and independently cacheable to disk (cache first, network only on cache miss — same pattern throughout).

- **GyDB (Gypsy Database)** — curated Ty1/copia RT and Gag, already extracted and lineage-classified. Source for the RT backbone and broad copia sampling.
- **REXdb** — curated retrotransposon protein-domain reference, lineage-classified. Second source for the same, for breadth and cross-checking.

Build the fetcher so these two are the only configured sources, but the source interface is clean enough that a third could be added later *if I explicitly decide to* — do not build speculative hooks for sources I have not asked for.

- **Special tips that need provenance care** (do NOT pull these blindly — flag them for me to supply exact identifiers, and only if the element is actually present in GyDB/REXdb):
  - **EVD / EVADÉ** (*A. thaliana*, ATCOPIA93 family). If included, it must be the specific evdGAG sequence characterized by Polkhovskiy et al. 2025, matched via the TAIR locus — NOT any ATCOPIA93 copy (ATTRAPE is a near-identical non-functional sibling and must not be confused for it). If GyDB/REXdb does not carry the correct copy, flag this to me rather than substituting a different one.
  - **Ty1** (*S. cerevisiae*, Gag annotated as **TYA**). Flag for me to confirm the exact GyDB/REXdb entry rather than assuming.

---

## 5. Pipeline stages (the actual flow)

Structure the repo so each stage is a discrete, independently runnable module with clear inputs/outputs on disk. Proposed stages — **but if you see a cleaner decomposition, propose it and ask before committing to it**:

1. **`fetch`** — pull element sequences from GyDB and REXdb (the only two sources); cache raw results to disk. Output: per-element records (accession, full protein sequence, source, lineage/taxonomy).
2. **`domain_gate`** — locate the RT region and the Gag region *within* each already-fetched protein, then slice them out. **Method, in strict preference order:**
   1. **Database-native domain labels/coordinates first.** If GyDB/REXdb already provides the element pre-sliced into domains (REXdb is organized by domain and lineage; GyDB ships per-domain sequences), use those labels directly — do not re-derive what the curated source already gives you.
   2. **Targeted profile search as the fallback**, for anything not pre-sliced. Scan against a **small, fixed set of two profiles only** — one RT profile and one Gag profile — NOT the full Pfam database. Prefer **retroelement-specific HMMs** (GyDB's copia/gypsy domain HMMs, or REXdb-derived profiles via the DANTE approach) over generic pan-life Pfam models, because they are trained on LTR-element alignments and carry far fewer spurious hits. Use HMMER (`hmmsearch`, one profile against all sequences) with **curated/strict cutoffs** (Pfam `--cut_ga` gathering thresholds if a Pfam-derived model is ever used, or a strict E-value otherwise).
   - **Do NOT scan against all of Pfam.** A full ~20,000-model `hmmscan` is exactly the noisy, false-positive-prone approach this design avoids (it is what produced spurious Gag hits on unrelated host genes in my earlier exploratory work). Only the two target profiles are in scope.
   - **Pfam IDs are not forbidden as objects** — a profile is a profile regardless of catalogue number — but they must be used as *targeted domain-locating rulers*, never as a *pull/search filter*, and never as part of an all-of-Pfam sweep. Retroelement-specific profiles are strongly preferred over Pfam ones.
   - **Gag N-terminal rule (critical):** the PrLD sits in the *disordered N-terminal region upstream of the structured capsid (CA) core*, and that N-terminus is exactly what domain models bound poorly. So do NOT clip the Gag slice at the annotated CA/Gag domain boundary. Extract the Gag slice **from the start of the Gag ORF through the annotated core**, so the disordered N-terminus is retained. Record both the raw domain coordinates and the extended slice boundaries. If the ORF start is ambiguous, flag it to me rather than guessing.
   - Filter out elements lacking an intact, translatable Gag ORF (drop pseudogenized/frameshifted copies), and log every drop.
   - Output: two FASTA sets (RT, Gag) keyed on the same element accessions, plus a coordinates table.
3. **`cluster`** — MMseqs2 identity clustering (to collapse near-identical copies / pseudoreplication). Retain both raw and clustered views non-destructively (tag every record with cluster ID + size; never silently delete). Clustering decisions (identity/coverage thresholds) come from config.
4. **`tree`** — MAFFT-align the RT slices → trimAl → IQ-TREE. Output: the RT phylogeny + support values.
5. **`score`** — score the Gag slices for the trait, tiered (see Section 6). Output: one trait table, one row per element, joined to accessions.
6. **`signal`** (R) — phylogenetic signal on the trait over the RT tree: Pagel's λ (continuous trait, e.g. `phytools::phylosig`) and Fritz & Purvis D (binary trait, e.g. `caper::phylo.d`). **This is a gate:** low signal means ancestral reconstruction is uninterpretable and the pipeline should say so rather than proceed blindly.
7. **`ancestral`** (R, conditional) — stochastic character mapping (e.g. `phytools::make.simmap`) under multiple models (ER/SYM/ARD) with model comparison and node uncertainty. Only meaningful if Stage 6 shows signal.
8. **`report`** — assemble the tree figure (e.g. `ggtree`) with the trait mapped on, plus a summary of signal/reconstruction results.

Also (can be a standalone utility, not necessarily inline): a **confound check** that scores a curated set of experimentally confirmed non-Q/N prions (e.g. mammalian PrP, plant Luminidependens, FUS-family) with PLAAC to quantify the false-negative rate of a yeast-trained predictor outside fungi. Structure the repo so this can be added; ask me whether to include it in v1.

---

## 6. Scientific parameters that MUST be config-driven — and DECISIONS I have to make

Put all of these in the config file. **Several are unresolved and I must decide them before the relevant stage runs. Where a value is marked `[DECISION NEEDED]`, do not guess a default — stop and ask me (see Section 9).**

- **PLAAC alpha (α):** `[DECISION NEEDED]`. Options are 0, 0.5, 1.0. For a comparative design a fixed reference (α = 1.0, fixed yeast background) is batch-independent and likely correct, but I need to confirm. The pipeline must also support running all three (0/0.5/1) as a sensitivity analysis. Do not hardcode 0.5.
- **PLAAC core length:** default 60 (standard); expose in config.
- **Protein length filter:** min/max aa bounds for pulls; expose in config (my old pipeline used 50–3000).
- **"Low phylogenetic signal" threshold:** `[DECISION NEEDED]`. The λ and/or D value below/above which I will declare signal absent. **This must be fixed BEFORE the tree is examined (pre-registration)** — the pipeline should read it from config and report the outcome against it, never let me pick it after seeing results.
- **"Many transitions" threshold** for ancestral reconstruction: `[DECISION NEEDED]`, same pre-registration logic.
- **MMseqs2 clustering thresholds:** identity and coverage (my old work used 95% identity / 80% coverage per domain) — expose in config, confirm with me.
- **Trait tiers** (which detectors are primary/secondary/tertiary):
  - Primary (continuous): PLAAC LLR / NLLR at the fixed α.
  - Secondary (compositional): Q/N fraction; LCD-Composer.
  - Tertiary (orthogonal): disorder (IUPred3 or metapredict — `[DECISION NEEDED]` which one), and LLPS propensity via **catGRANULE** (note: there are two versions, original and catGRANULE 2.0 / ROBOT — `[DECISION NEEDED]` which one; this affects reproducibility and overlap with the disorder layer).
  - Binary PrLD call reported as a *derived* summary of the continuous score, with α-sensitivity.
- **Random seed:** fixed, in config.
- **Domain-gating profiles:** the source of the RT and Gag HMM profiles (GyDB copia/gypsy HMMs vs REXdb/DANTE-derived) should be config-driven and documented. Retroelement-specific profiles preferred; full-Pfam scanning is forbidden (see Stage 2). If the correct profile source is unclear for a given lineage, `[DECISION NEEDED]` — ask me.
- **Taxonomic sampling rule:** `[DECISION NEEDED]`. How many elements per copia lineage per host taxon, and the sampling frame (GyDB/REXdb lineage classification). This is a scientific sampling decision, not a default you should invent.
- **Scope / outgroup:** whether to include Ty3/gypsy as an outgroup comparison — `[DECISION NEEDED]`.
- **Unit of analysis for any Saccharomyces sub-analysis:** strains of one species vs. species across a clade — `[DECISION NEEDED]` if we include that arm.

---

## 7. Repository structure and hygiene

- Sensible layout: e.g. `src/` (or a package dir) for Python modules, `R/` for the R comparative scripts, `config/` for the config file(s), `data/` (gitignored, with a documented sub-structure for raw/cache/interim), `results/` (gitignored outputs), `envs/` for environment files, `docs/` if needed. Propose the exact layout and ask if unsure.
- **`.gitignore`** that excludes data, caches, results, large binaries, and the PLAAC jar — with a `data/README.md` explaining how to obtain/regenerate what's gitignored.
- **`README.md`** as the front door: the research question, the two hypotheses (H1 turnover vs H0 acquisition), a diagram or ordered list of the stages, install instructions (conda env + external tools), a quickstart (how to run one stage and how to run everything), and a clear statement of which config decisions the user must make before running.
- **`environment.yml`** (pinned) + **`requirements.txt`** (pinned).
- A **`LICENSE`** (ask me which; likely MIT) and a **`CITATION.cff`** stub.
- **No secrets, no absolute paths, no committed data.** Any contact email used in an HTTP User-Agent header for GyDB/REXdb requests should come from config, not be hardcoded.
- Logging to stdout with levels; a `--verbose` flag; a `--dry-run` mode that prints planned actions (esp. external tool commands) without executing.
- Each stage idempotent and re-runnable from cache with a `--skip-fetch`-style option where network is involved.

### Git / commit rules

- **Do NOT add yourself as an author, co-author, or committer on any commit.** All commits are authored solely by me. Do not add `Co-Authored-By: Claude` (or any Claude/Anthropic identity) trailer, do not add "Generated with Claude Code" or any similar line to commit messages, and do not set author/committer metadata to anything other than my existing git config. Commit messages should describe the change only, with no attribution to an AI tool.
- If you are unsure whether an action would add such attribution, don't do it — ask me.

---

## 8. What NOT to do

- Do **not** build the tree from Gag. RT builds the tree; Gag is scored. If you find yourself aligning Gag across the whole superfamily for the backbone, stop — that's the error this whole design exists to avoid.
- Do **not** scan against the full Pfam database in `domain_gate`. Use only the two targeted RT/Gag profiles (retroelement-specific preferred), and prefer database-native domain coordinates over any HMM search. Do **not** use Pfam (or any profile) as a *pull/search filter* — profiles are domain-locating rulers only.
- Do **not** clip the Gag slice at the CA domain boundary — extract from the ORF start so the disordered N-terminal PrLD region is retained.
- Do **not** hardcode α = 0.5, or any `[DECISION NEEDED]` value. Ask.
- Do **not** carry over Windows-specific or machine-specific paths.
- Do **not** silently drop data during clustering/dedup — tag and retain, keep it auditable.
- Do **not** let ancestral reconstruction run as if valid when phylogenetic signal is absent — the signal stage gates it.
- Do **not** commit large data, caches, or the PLAAC jar.
- Do **not** invent scientific defaults for sampling, thresholds, or detector choices — those are my calls.
- Do **not** pull sequences from anywhere other than GyDB and REXdb. No UniProt, NCBI, or genome mining.
- Do **not** add any AI/Claude attribution to git commits — no co-author trailer, no "generated with" line, no modified author metadata. Commits are mine alone.

---

## 9. BEFORE YOU BUILD: ask me these (and anything else unclear)

Do not scaffold the repo until we've resolved the questions below. Ask them a few at a time, wait for my answers, and if my answers surface new ambiguities or tradeoffs, raise those too. I would rather answer questions now than refactor later. Treat any `[DECISION NEEDED]` in Section 6 as a required question. At minimum, confirm with me:

1. **Repository name and license** — I'm deciding between a self-documenting name (e.g. `copia-gag-prld-phylo`) and a project-scoped name (e.g. `prioneers-phylo`). Ask which, and which license.
2. **Scope of v1** — do we build all 8 stages now, or a runnable core (fetch → domain_gate → tree → score) first, with signal/ancestral/report and the confound check as a documented second phase? I lean toward a solid core first — confirm.
3. **The `[DECISION NEEDED]` scientific parameters** in Section 6, especially: PLAAC α; the pre-registered signal and transition thresholds; MMseqs2 thresholds; disorder tool (IUPred3 vs metapredict); catGRANULE version; taxonomic sampling rule; Ty3/gypsy outgroup yes/no. If I say "I haven't decided yet" for any of these, DO NOT default them — instead, scaffold the config with a clearly-marked `TODO: DECISION REQUIRED` placeholder that makes the pipeline refuse to run that stage until I fill it in.
4. **Sources** — confirmed: GyDB and REXdb ONLY. Do not propose or add any other source. If during the build you hit something that seems to *require* another source (e.g. a lineage or a specific tip not present in either DB), stop and ask me rather than adding one.
5. **Special tips** — do you want me to hold slots for EVD and Ty1 and prompt you for their exact GyDB/REXdb identifiers (and flag if the correct copy isn't present) rather than auto-pulling or substituting them? I think yes — confirm.
6. **Environment strategy** — conda (`environment.yml`) as primary acceptable? Any tools I should expect to install/run manually at first (e.g. catGRANULE as a web tool) that the pipeline should accommodate as an external input rather than orchestrate directly?
7. **Domain-gating profiles** — confirm the RT and Gag profile source (GyDB copia/gypsy HMMs vs REXdb/DANTE-derived). If GyDB/REXdb already provides pre-sliced domains for most elements, confirm you should use those directly and treat `hmmsearch` as a fallback only. Do not scan full Pfam under any circumstances.
8. **Anything in this brief that is scientifically or technically ambiguous, internally inconsistent, or that you think is the wrong call** — tell me before building, not after. Push back if a design choice here looks weak.

Once I've answered these, propose the final repo structure and the config schema for my sign-off, and only then start writing files.
