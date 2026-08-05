# Ty1/copia Gag PrLD phylogenetics

Does the disordered N-terminal region of Ty1/copia Gag stay **conserved across
the superfamily** while its **prion-like composition turns over independently**
(H1), rather than the PrLD being a single inherited innovation (H0)?

**The design principle:** the tree and the trait come from two different regions
of the same element. **RT is aligned** to build the phylogeny — it is the only
part alignable superfamily-wide. **Gag is scored** to make the trait. They are
joined by element accession. The tree is never built from Gag.

## Status

The pipeline is being rebuilt one decision at a time, with each scientific
choice reasoned from the data and signed off before code is written against it.
See `claude/DECISIONS.md` for what has been settled and what is still open.

`fetch` is implemented. Nothing downstream is.

## Layout

```
data/REXdb/    the release files (not committed — see below)
data/          generated sequence sets
scripts/       the pipeline
claude/        AI instruction docs and the decisions log — deletable after dev
```

## Data

Download from <https://github.com/repeatexplorer/rexdb> into `data/REXdb/`:

- `Viridiplantae_v4.0.fasta` — protein domain slices
- `Viridiplantae_v4.0.classification` — superfamily and lineage
- `Viridiplantae_v4.0_ALL_info-species-source_taxid` — species and NCBI taxid
- `Viridiplantae_v4.0_ALL_DNA.fasta` — full element nucleotide sequences

Citation: Neumann et al., *Mobile DNA* 2019, doi:10.1186/s13100-018-0144-1.

Note that REXdb's Gag slices are the capsid core, not the whole ORF. The
disordered N-terminus is **not** in the protein FASTA — it is recovered from the
element DNA by `fetch`.

## Install

The bioinformatics tools have no Windows builds, so on Windows this runs under
WSL. PLAAC is built separately — see `plaac/README.md`.

```bash
conda env create -f environment.yml
conda activate copia-prld
```

## Run

```bash
python  scripts/fetch.py     # sequence sets from REXdb        (stdlib only)
python  scripts/score.py     # trait table: composition, PLAAC, disorder
python  scripts/tree.py      # sample -> MAFFT -> trimAl -> IQ-TREE
Rscript scripts/plots.R      # figures
```

`fetch` reads `data/REXdb/` and writes to `data/`:

| file | contents |
| --- | --- |
| `elements.tsv` | one row per copia/gypsy element, metadata and flags |
| `rt_copia.faa`, `rt_gypsy.faa` | RT slices as published — these build the tree |
| `gag_copia.faa`, `gag_gypsy.faa` | Gag N-terminal region + capsid core |

RT is used as REXdb publishes it. Gag is re-derived from element DNA, because
the published Gag slice is the capsid core only and does not contain the
N-terminal region this project is about: the slice is located in a six-frame
translation to fix the reading frame, then extended upstream to the in-frame
stop. There is no length cap: capping at 75 aa removed exactly the long
N-terminus elements that carry the signal (F-18).

The later stages add `traits.tsv` (one row per Gag sequence: composition,
charge, PLAAC at three alpha values, metapredict disorder), `tree/rt.treefile`
and `figures/`.

Every element with an RT enters the tree, so the tree is larger than the trait
set — prune it on `gag_status` for the comparative analysis. Nothing is filtered
on unresolved residues; the counts go in `rt_ambiguous` and `gag_ambiguous` so
scoring can decide. Copia and gypsy are never merged.

Current result, from 14,355 elements:

| | elements | RT | Gag | `no_slice` | `unplaceable` | `no_nterm` | `record_mismatch` |
| --- | --- | --- | --- | --- | --- | --- | --- |
| copia | 5,606 | 5,596 | 5,042 | 442 | 79 | 33 | 10 |
| gypsy | 8,749 | 8,748 | 8,486 | 167 | 84 | 11 | 1 |

`gag_status` says why a trait is missing. `no_slice` — REXdb published no Gag.
`unplaceable` — it did, but the slice isn't contiguous in the element's own DNA,
usually a frameshifted pseudogene. `record_mismatch` — no part of the slice is
in the DNA at all, so the two release files disagree; these are the only
elements also dropped from the tree. `gag_core_exact` flags cores placed by
their first 30 aa rather than in full.

## Later stages

Not built: cluster (MMseqs2), sampling, tree (MAFFT → trimAl → IQ-TREE), score
(PLAAC, charge, disorder, LLPS), phylogenetic signal, ancestral reconstruction,
report.
