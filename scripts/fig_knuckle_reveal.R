# Two slides of the same figure: the bars, then the same bars split by lineage.
#
#   fig_reveal_1_plain.{png,pdf}     bars as in fig_knuckle_matched_db
#   fig_reveal_2_lineage.{png,pdf}   identical bars, each subdivided by lineage
#
# The two are drawn by the same function with one flag, so every bar, tick,
# label, star and margin lands in the same place. Advancing the slide changes
# only the fill, which is what makes the reveal read as one picture gaining
# detail rather than two charts being compared.
#
# The bar heights are unchanged between slides: each bar is still the share of
# that group carrying a zinc knuckle, within that length window. The stack
# divides that same height by lineage -- the segments sum to the bar, they do
# not rescale it.
#
# Each superfamily gets its own palette and legend. A shared legend would need
# eight distinguishable hues for eight lineages, and eight cannot be told apart
# reliably; the panels are separate objects so each can name its own.
#
# Run from the repo root, in the conda env:
#   Rscript scripts/fig_knuckle_reveal.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

OUT <- "data/gag_plaac/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"; GUIDE <- "#e6e4df"
WITH_PRLD <- "#7B2D9E"; NO_PRLD <- "#8a8a85"
HUES <- c("#D2551B", "#0E8A5F", "#7B2D9E", "#2D7DD2")
OTHER <- "#a9a8a3"

MIN_CELL <- 20
BREAKS <- c(0, 300, 400, 550, 750, Inf)
BAND_LABELS <- c("<300", "300-400", "400-550", "550-750", ">750")
TOP_LIN <- 4
OFFSET <- 0.18      # half the gap between the two bars of a pair
BAR_W <- 0.33

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

# a band is comparable only where both groups are populated -- same rule as
# fig_knuckle_matched.R, so the two figures show the same windows
usable <- rex %>%
  count(set, band, prd) %>%
  group_by(set, band) %>% filter(n() == 2, min(n) >= MIN_CELL) %>%
  ungroup() %>% distinct(set, band)

dat <- semi_join(rex, usable, by = c("set", "band"))

groups <- dat %>%
  group_by(set, band, prd) %>%
  summarise(n = n(), k = sum(knuckle), .groups = "drop") %>%
  mutate(p = k / n,
         lo = pmax(0, p - 1.96 * sqrt(p * (1 - p) / n)),
         hi = pmin(1, p + 1.96 * sqrt(p * (1 - p) / n)))

tests <- groups %>%
  select(set, band, prd, n, k) %>%
  pivot_wider(names_from = prd, values_from = c(n, k)) %>%
  rowwise() %>%
  mutate(p_value = fisher.test(matrix(c(k_TRUE, n_TRUE - k_TRUE,
                                        k_FALSE, n_FALSE - k_FALSE),
                                      nrow = 2))$p.value) %>%
  ungroup() %>%
  mutate(p_holm = p.adjust(p_value, "holm"),
         label = cut(p_holm, c(-Inf, 0.001, 0.01, 0.05, Inf),
                     labels = c("***", "**", "*", "ns")))

# Named lineages are the most ABUNDANT ones in the plotted bands, not the ones
# contributing most knuckles. Ranking by knuckle count hides exactly the lineage
# this analysis is about: Ale is the largest copia lineage but only ~9%
# knuckle-positive, so it ranks below Ivana, Angela, Tork and SIRE and would
# disappear into "other".
top_lin <- dat %>% count(set, lin, sort = TRUE) %>%
  group_by(set) %>% slice_head(n = TOP_LIN) %>% ungroup()

# the knuckle-bearing members of each group, split by lineage. share is computed
# against the same denominator as the plain bar, so the segments sum to it.
stack <- dat %>%
  filter(knuckle) %>%
  mutate(lineage_shown = ifelse(paste(set, lin) %in% paste(top_lin$set, top_lin$lin),
                                lin, "other lineages")) %>%
  count(set, band, prd, lineage_shown, name = "seg") %>%
  left_join(select(groups, set, band, prd, n), by = c("set", "band", "prd")) %>%
  mutate(share = seg / n)


panel <- function(which_set, split) {
  present <- BAND_LABELS[BAND_LABELS %in%
                           unique(as.character(groups$band[groups$set == which_set]))]
  at <- function(b) match(as.character(b), present)
  shift <- function(prd) ifelse(prd, -OFFSET, OFFSET)

  g <- filter(groups, set == which_set) %>% mutate(x = at(band) + shift(prd))
  t <- filter(tests, set == which_set) %>% mutate(x = at(band))
  s <- filter(stack, set == which_set) %>% mutate(x = at(band) + shift(prd))

  named <- unique(top_lin$lin[top_lin$set == which_set])
  named <- named[named %in% s$lineage_shown]
  colours <- setNames(c(HUES[seq_along(named)], OTHER), c(named, "other lineages"))
  s <- mutate(s, lineage_shown = factor(as.character(lineage_shown),
                                        levels = rev(c(named, "other lineages"))))

  status_fill <- c("carries a prion-like domain" = WITH_PRLD,
                   "no prion-like domain" = NO_PRLD)
  g <- mutate(g, status = factor(ifelse(prd, "carries a prion-like domain",
                                        "no prion-like domain"),
                                 levels = names(status_fill)))

  p <- ggplot()
  if (split) {
    p <- p +
      geom_col(data = s, aes(x, share, fill = lineage_shown), width = BAR_W,
               colour = SURFACE, linewidth = 0.4) +
      scale_fill_manual(values = colours, name = which_set,
                        breaks = c(named, "other lineages"))
  } else {
    p <- p +
      geom_col(data = g, aes(x, p, fill = status), width = BAR_W) +
      scale_fill_manual(values = status_fill, name = which_set)
  }

  # which bar is which has to survive the reveal: on the second slide the fill
  # becomes lineage, so the group is carried by an outline that is present on
  # both slides and unchanged between them
  p <- p +
    geom_col(data = g, aes(x, p, colour = status), fill = NA, width = BAR_W,
             linewidth = 0.75) +
    scale_colour_manual(values = status_fill, name = NULL,
                        guide = guide_legend(order = 2, ncol = 1))

  p +
    geom_errorbar(data = g, aes(x = x, ymin = lo, ymax = hi), width = 0.09,
                  linewidth = 0.35, colour = INK2) +
    geom_text(data = g, aes(x = x, y = hi + 0.018, label = n), size = 2.3,
              colour = INK2) +
    geom_text(data = left_join(t, g %>% group_by(band) %>%
                                 summarise(top = max(hi), .groups = "drop"),
                               by = "band"),
              aes(x = x, y = pmin(top + 0.085, 1.06), label = label),
              size = ifelse(t$p_holm < 0.05, 4.4, 2.7),
              vjust = ifelse(t$p_holm < 0.05, 0.65, 0),
              colour = ifelse(t$p_holm < 0.05, INK, INK2),
              fontface = ifelse(t$p_holm < 0.05, "bold", "italic")) +
    scale_x_continuous(breaks = seq_along(present), labels = present,
                       limits = c(0.4, length(present) + 0.6)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.12),
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "Gag length (aa)", y = "share carrying a zinc knuckle") +
    theme_minimal(base_size = 10) +
    theme(text = element_text(colour = INK),
          axis.title = element_text(size = 8.6, colour = INK2),
          axis.text = element_text(size = 8, colour = INK2),
          panel.grid = element_blank(),
          axis.line.y = element_line(colour = GUIDE, linewidth = 0.4),
          axis.ticks.y = element_line(colour = GUIDE, linewidth = 0.4),
          axis.ticks.length = unit(2.5, "pt"),
          legend.position = "right",
          legend.title = element_text(face = "bold", size = 8.8, colour = INK2),
          legend.text = element_text(size = 8.4),
          legend.key.size = unit(10, "pt"),
          plot.background = element_rect(fill = SURFACE, colour = NA),
          plot.margin = margin(8, 10, 4, 10))
}


build <- function(split, file) {
  fig <- panel("REXdb  Ty1/copia", split) | panel("REXdb  Ty3/gypsy", split)
  ggsave(file.path(OUT, paste0(file, ".png")), fig, width = 11, height = 4.6,
         dpi = 400, bg = SURFACE)
  ggsave(file.path(OUT, paste0(file, ".pdf")), fig, width = 11, height = 4.6,
         bg = SURFACE)
  cat("wrote", file, "\n")
}

build(FALSE, "fig_reveal_1_plain")
build(TRUE, "fig_reveal_2_lineage")

cat("\nbar heights are identical between the two slides; the stack subdivides them\n\n")
print(as.data.frame(
  stack %>%
    left_join(select(groups, set, band, prd, bar_height = p),
              by = c("set", "band", "prd")) %>%
    group_by(set, band, prd) %>%
    mutate(check = sum(share)) %>%
    summarise(bar = round(first(bar_height), 4),
              stack_sums_to = round(first(check), 4), .groups = "drop")))
