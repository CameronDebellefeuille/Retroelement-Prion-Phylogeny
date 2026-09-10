# Who supplies the trade-off cell? The same bands, broken down by lineage.
#
#   fig_knuckle_composition.{png,pdf}
#
# Same panels and the same length windows as fig_knuckle_matched.R, but the
# question is different. There the bars asked "what share of this group carries a
# knuckle". Here each bar is the share of all Gags in that length band falling
# into one of the two prion-like-domain cells, and the stack shows which lineages
# supply it:
#
#   left bar   domain, no knuckle   -- the cell the trade-off hypothesis is about
#   right bar  domain and knuckle   -- the cell it says should be rare
#
# The point of the stacking is the objection that the effect is one clade counted
# many times. If the left bar were solid Ale, the anti-correlation would be a fact
# about Ale rather than about Gag. Reading the composition is the answer to that,
# and it is a partial answer only: lineages are not independent of each other
# either, and nothing here corrects for the phylogeny.
#
# Bars do not sum to 100%: Gags with no domain at all are the remainder and are
# not drawn, because the domain-free majority would flatten everything else.
#
# Run from the repo root, in the conda env:
#   Rscript scripts/fig_knuckle_composition.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork)
})

OUT <- "data/gag_plaac/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"; GUIDE <- "#e6e4df"
OTHER <- "#a9a8a3"
# four hues validated all-pairs for protan/deutan separation, plus grey for the tail
HUES <- c("#D2551B", "#0E8A5F", "#7B2D9E", "#2D7DD2")

MIN_CELL <- 20
BREAKS <- c(0, 300, 400, 550, 750, Inf)
BAND_LABELS <- c("<300", "300-400", "400-550", "550-750", ">750")
TOP_N <- 4

cls <- strsplit(readLines("data/REXdb/Viridiplantae_v4.0.classification"), "\t")
cls <- cls[vapply(cls, function(x) length(x) >= 5 && x[3] == "LTR", logical(1))]
lineage <- tibble(id = vapply(cls, `[`, "", 1),
                  lin = vapply(cls, function(x) x[length(x)], ""))

rex <- read_tsv("data/gag_full/gag_full_plaac.tsv", show_col_types = FALSE) %>%
  select(id = rexdb_id, superfamily, gag_len, has_prd) %>%
  inner_join(read_tsv("data/gag_full/nc_hits.tsv", show_col_types = FALSE) %>%
               select(id = rexdb_id, nc_combined), by = "id") %>%
  left_join(lineage, by = "id") %>%
  mutate(len = as.numeric(gag_len), prd = has_prd == 1, knuckle = nc_combined == 1,
         band = cut(len, BREAKS, labels = BAND_LABELS, right = TRUE),
         set = ifelse(superfamily == "Ty1/copia", "REXdb  Ty1/copia",
                      "REXdb  Ty3/gypsy")) %>%
  filter(!is.na(band), !is.na(lin), lin != "")

# the bands that fig_knuckle_matched.R was able to compare, so the two figures
# show the same windows and nothing appears here that was untestable there
usable <- rex %>%
  group_by(set, band, prd) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(set, band) %>% filter(n() == 2, min(n) >= MIN_CELL) %>%
  ungroup() %>% distinct(set, band)

dat <- rex %>% semi_join(usable, by = c("set", "band"))

# the lineages that actually supply domain-bearing Gags, per superfamily
top <- dat %>% filter(prd) %>% count(set, lin, sort = TRUE) %>%
  group_by(set) %>% slice_head(n = TOP_N) %>% ungroup()

dat <- dat %>%
  mutate(lineage_shown = ifelse(paste(set, lin) %in% paste(top$set, top$lin),
                                lin, "other lineages"))

totals <- dat %>% count(set, band, name = "band_n")

comp <- dat %>%
  filter(prd) %>%
  mutate(cell = ifelse(knuckle, "domain and knuckle", "domain, no knuckle")) %>%
  count(set, band, cell, lineage_shown, name = "k") %>%
  left_join(totals, by = c("set", "band")) %>%
  mutate(share = k / band_n,
         cell = factor(cell, levels = c("domain, no knuckle", "domain and knuckle")))

# One panel per superfamily, each with its own palette and legend. A shared
# legend would reuse the same four hues for eight different lineages -- eight
# hues cannot be told apart reliably, so the split is the honest way to keep the
# colours meaningful.
panel_for <- function(which_set) {
  d <- filter(comp, set == which_set)
  named <- setdiff(unique(as.character(d$lineage_shown)), "other lineages")
  named <- named[order(match(named, unique(top$lin[top$set == which_set])))]
  colours <- setNames(c(HUES[seq_along(named)], OTHER), c(named, "other lineages"))
  d <- mutate(d, lineage_shown = factor(as.character(lineage_shown),
                                        levels = c(named, "other lineages")))
  tops <- d %>% group_by(band, cell) %>%
    summarise(total = sum(share), k = sum(k), .groups = "drop")

  ggplot(d, aes(cell, share, fill = lineage_shown)) +
    geom_col(width = 0.66, colour = SURFACE, linewidth = 0.5) +
    geom_text(data = tops, aes(x = cell, y = total + 0.012,
                               label = sprintf("%.1f%%  (%d)", 100 * total, k)),
              inherit.aes = FALSE, size = 2.5, colour = INK2) +
    facet_grid(. ~ band) +
    scale_fill_manual(values = colours, name = which_set) +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0, 0.14))) +
    scale_x_discrete(labels = function(x) sub(", ", ",
", sub(" and ", "
and ", x))) +
    labs(x = NULL, y = "share of all Gags in that length band") +
    theme_minimal(base_size = 10) +
    theme(text = element_text(colour = INK),
          strip.text.x = element_text(face = "bold", size = 8.8, colour = INK2),
          axis.title.y = element_text(size = 8.4, colour = INK2),
          axis.text.x = element_text(size = 7.4, colour = INK2, lineheight = 0.95),
          axis.text.y = element_text(size = 8, colour = INK2),
          panel.grid = element_blank(),
          axis.line.y = element_line(colour = GUIDE, linewidth = 0.4),
          axis.ticks.y = element_line(colour = GUIDE, linewidth = 0.4),
          axis.ticks.length = unit(2.5, "pt"),
          panel.spacing.x = unit(8, "pt"),
          legend.position = "right",
          legend.title = element_text(face = "bold", size = 8.8, colour = INK2),
          legend.text = element_text(size = 8.4),
          legend.key.size = unit(10, "pt"),
          plot.background = element_rect(fill = SURFACE, colour = NA),
          plot.margin = margin(6, 10, 4, 10))
}

fig <- panel_for("REXdb  Ty1/copia") / panel_for("REXdb  Ty3/gypsy")

ggsave(file.path(OUT, "fig_knuckle_composition.png"), fig, width = 10.5,
       height = 5.6, dpi = 400, bg = SURFACE)
ggsave(file.path(OUT, "fig_knuckle_composition.pdf"), fig, width = 10.5,
       height = 5.6, bg = SURFACE)
cat("wrote fig_knuckle_composition\n\n")

cat("lineage share of the domain-bearing cells:\n")
print(as.data.frame(
  comp %>% group_by(set, cell, lineage_shown) %>%
    summarise(elements = sum(k), .groups = "drop_last") %>%
    mutate(pct_of_cell = round(100 * elements / sum(elements), 1)) %>%
    ungroup() %>% arrange(set, cell, desc(elements))))
