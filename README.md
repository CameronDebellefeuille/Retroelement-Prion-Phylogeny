# LTR Retrotransposon Gag PrLD phylogenetics

Data is GyDB. The REXdb analysis is in progress.

## Install

```bash
conda env create -f environment.yml
conda activate copia-prld
pip install --no-deps metapredict==2.65.1
```

PLAAC is vendored unmodified from <https://github.com/whitehead/plaac> at tag
`v1.1.0b7`; `plaac/README.md` is upstream's own. Rebuild with
`cd plaac && PLAAC_VERSION=v1.1.0b7 ./build_plaac.sh`, then move
`target/plaac.jar` to `plaac/plaac.jar`, where the scripts look for it.
Three example proteomes (27 MB) are not kept; `plaac/CITATION.cff` is from
upstream master, as it postdates the tag.

## Pipeline

```
cores-database ──> plaac_all_cores.py ──> all_cores_a05.tsv ──┐
                                                              │
cores-database ──┐                                            │
superfamily.csv ─┴─────────> prep_272.py <────────────────────┤
                                  │                           │
                                  ├─> rt_272.faa ──> MAFFT → trimAl → IQ-TREE ──> rt.treefile
                                  ├─> gag_272.faa ──> disorder_272.py ──> disorder_272.tsv
                                  └─> traits_272.tsv          │
                                                              │
cores-database ───────────────────────────────────────────────┴─> disorder_all_cores.py ──> disorder_all.tsv
```
## Scripts
```bash
python  scripts/plaac_all_cores.py        # PLAAC a=0.5 over all 2,636 GyDB protein cores
python  scripts/prep_272.py               # creates the faa for the 272 gag-rt pairs for building the tree and PrLD mapping
python  scripts/disorder_272.py           # metapredict, the 272 Gag cores
python  scripts/disorder_all_cores.py     # metapredict, every core
Rscript scripts/fig_gydb_rt_layers.R      # then the five other fig_*/figures_* scripts
```

## Files

- **`data/raw/gydb/`**: (the source).
- **`cores-database`**: is GyDB's 2,636 protein cores
- `gydb_superfamily_host.csv`: is a hand-curated sheet mapping element to
superfamily extracted from the GyDB website
- **`data/raw/rexdb/`**: REXdb v4.0 release files (not used)
- **`data/processed/gydb/`**:

| file | what it is |
| --- | --- |
| `rt_272.faa` | RT cores: the tree is built from these |
| `gag_272.faa` | the matching Gag cores, scored but never aligned |
| `traits_272.tsv` | the tip table: superfamily, host, clade, length, LLR, PrLD coords, knuckle |
| `copia_tree_tips.txt` | which elements also sit on the REXdb copia tree |
| `tree_272/` | `rt.treefile` plus the alignment, trimmed alignment and IQ-TREE log behind it |
| `plaac_data/all_cores_a05.tsv` | PLAAC at α = 0.5 over all cores |
| `disorder_272.tsv`, `disorder_all.tsv` | metapredict, per Gag core and per class |
| `InterProScan/part1–4.tsv` | full domain calls for 327 elements; the zinc-knuckle source |

## Figures

```
rt.treefile + traits_272.tsv                     ──> fig_gydb_rt_layers.R
rt.treefile + traits_272.tsv + disorder_272.tsv  ──> figures_combined.R
traits_272.tsv + InterProScan/                   ──> fig_knuckle_superfamily_dark.R
traits_272.tsv                                   ──> fig_prld_position_dark.R
all_cores_a05.tsv                                ──> fig_gydb_classes_dark.R
disorder_all.tsv                                 ──> figures_disorder_classes.R
```
## Notes & Sources 

**Superfamily is only partly curated.** The sheet covers 234 of the 336
Gag-bearing elements. The rest are assigned `CopiaSL*` / `GypsySL*` in the GyDB.
In my data I've labeled them by `classify()` in `prep_272.py`, and those tips carry no information for host or clade.
Of the 272 tips, 231 are curated matches and 39 are name-assigned; the 39 are the ones with a blank `host`.

Of the 333 elements with both cores, the 61 that don't reach the tree were all
dropped for one reason: an `X` in the sequence.


Citations: GyDB, Llorens et al. *NAR* 2011. PLAAC, Lancaster et al.
*Bioinformatics* 2014. metapredict, Emenecker et al. *Biophys J* 2021.

## AI Disclaimer

Claude Code was used to assist with writing and debugging analysis scripts.
I designed the analyses, verified all outputs, and am responsible for the interpretations presented here.
