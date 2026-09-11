# Ty1/copia Gag PrLD phylogenetics

Data is GyDB. The REXdb arm of the project is parked; its raw files are still
here so it can be picked up again.

## Install

Windows. R comes from the environment too, so `ggtree` and the fonts are pinned
rather than borrowed from a system install.

```bash
conda env create -f environment.yml
conda activate copia-prld
pip install --no-deps metapredict==2.65.1
```

PLAAC is a Java jar built from source: see `plaac/README.md`.

## Pipeline

```
cores-database ──┐
                 ├─> prep_272.py ──> rt_272.faa ──> [MAFFT → trimAl → IQ-TREE] ──> rt.treefile
superfamily.csv ─┘                   gag_272.faa                                        │
                                     traits_272.tsv <───────────────────────────────────┘

cores-database ──> plaac_all_cores.py ──> all_cores_a05.tsv
               └─> disorder_*.py       ──> disorder_272.tsv, disorder_all.tsv
```

`prep_272.py` takes every element carrying a clean GAG **and** RT core (no `X`) —
272 of them — and writes the tree input, the Gag set, and the tip table. The six
figure scripts read `traits_272.tsv` and `rt.treefile`.

The tree build itself is the one step Windows can't do: `mafft`, `trimal` and
`iqtree` have no win-64 builds. Versions and commands are at the top of
`environment.yml`, to run under WSL. `rt.treefile` is committed, so this is only
for rebuilding the phylogeny from scratch.

```bash
python  scripts/plaac_all_cores.py        # PLAAC a=0.5 over all 2,636 cores
python  scripts/prep_272.py               # the 272 set
python  scripts/disorder_272.py           # metapredict, the 272 Gag cores
python  scripts/disorder_all_cores.py     # metapredict, every core
Rscript scripts/fig_gydb_rt_layers.R      # then the five other fig_*/figures_* scripts
```

## Files

**`data/raw/gydb/`** — the source. `cores-database` is GyDB's 2,636 protein
cores, named `DOMAIN_element` (`GAG_17.6`, `RT_17.6`), covering 779 elements.
`gydb_superfamily_host.csv` is a hand-curated sheet mapping element to
superfamily, host and clade — not a GyDB download, and the only copy anywhere.

**`data/raw/rexdb/`** — REXdb v4.0 release files. Parked, unread by any script.

**`data/processed/gydb/`**

| file | what it is |
| --- | --- |
| `rt_272.faa` | RT cores — the tree is built from these |
| `gag_272.faa` | the matching Gag cores, scored but never aligned |
| `traits_272.tsv` | the tip table: superfamily, host, clade, length, LLR, PrLD coords, knuckle |
| `copia_tree_tips.txt` | which elements also sit on the REXdb copia tree |
| `tree_272/` | `rt.treefile` plus the alignment, trimmed alignment and IQ-TREE log behind it |
| `plaac_data/all_cores_a05.tsv` | PLAAC at α = 0.5 over all cores |
| `disorder_272.tsv`, `disorder_all.tsv` | metapredict, per Gag core and per class |
| `InterProScan/part1–4.tsv` | full domain calls for 327 elements; the zinc-knuckle source |

**`scripts/`** — four Python stages that build tables, six R scripts that draw
figures. **`figures/`** — the six PNGs those produce. **`plaac/`** — the jar and
the source to rebuild it.

**Superfamily is only partly curated.** The sheet covers 234 of the 336
Gag-bearing elements. The rest are `CopiaSL*` / `GypsySL*` entries typed from
their own filename by `classify()` in `prep_272.py` — sound, but not curation,
and those tips carry no host or clade. Of the 272 tips, 231 are curated matches
and 39 are name-assigned; the 39 are exactly the ones with a blank `host`.

Of the 333 elements with both cores, the 61 that don't reach the tree were all
dropped for one reason: an `X` in the sequence.

Citations: GyDB, Llorens et al. *NAR* 2011. PLAAC, Lancaster et al.
*Bioinformatics* 2014. metapredict, Emenecker et al. *Biophys J* 2021.
