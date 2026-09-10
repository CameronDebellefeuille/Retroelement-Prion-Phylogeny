# The knuckle / prion-like domain comparison, held at matched Gag length.
#
#   fig_knuckle_matched_db.{png,pdf}       by database
#   fig_knuckle_matched_lineage.{png,pdf}  by superfamily and lineage
#
# Why this exists: the raw share of knuckle-bearing Gags is confounded, and in
# the direction that hides the effect. Gags with a prion-like domain are much
# longer (Athila 651 aa against 371), and longer Gags carry a knuckle far more
# often for the trivial reason that there is more sequence to carry one in
# (Athila 0% to 24% across length quartiles, CRM 3% to 95%). Comparing raw shares
# therefore gives a knuckle-rate bonus to exactly the group expected to have
# fewer, and in four of eight lineages it reverses the sign of the result.
#
# Fixing that needs nothing more elaborate than comparing like with like: within
# a fixed window of Gag length, what share of each group carries a knuckle. The
# windows are absolute, not per-group quantiles, so a band means the same thing
# in every panel and across both databases.
#
# A band is drawn only when both groups have at least MIN_CELL elements in it.
# Empty space in these panels is honest: it is where one group has too few
# elements of that length to compare.
#
# Each band carries a test. Not a t-test: the measurement per element is binary
# -- knuckle or no knuckle -- so there are no means to compare and no variance to
# pool. The corresponding test for two proportions is Fisher's exact on the 2x2
# of knuckle presence against domain presence, which is exact at any cell count
# and so does not need the large-sample assumption a chi-squared would. Holm
# adjusted p-values across the bands of each figure are printed to the console;
# the figure itself shows the unadjusted value.
#
# Titles are deliberately absent -- these go onto slides that carry their own.
#
# Knuckle calls: REXdb by Pfam E<1e-3 or PROSITE PS50158 (92% sensitivity, 99%
# specificity against InterProScan); GyDB by InterProScan itself.
#
# Run from the repo root, in the conda env:
#   Rscript scripts/fig_knuckle_matched.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

OUT <- "data/gag_plaac/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"; GUIDE <- "#e6e4df"
WITH_PRLD <- "#7B2D9E"
NO_PRLD <- "#8a8a85"

MIN_CELL <- 20
BREAKS <- c(0, 300, 400, 550, 750, Inf)
BAND_LABELS <- c("<300", "300-400", "400-550", "550-750", ">750")

cls <- strsplit(readLines("data/REXdb/Viridiplantae_v4.0.classification"), "\t")
cls <- cls[vapply(cls, function(x) length(x) >= 5 && x[3] == "LTR", logical(1))]
lineage <- tibble(id = vapply(cls, `[`, "", 1),
                  lin = vapply(cls, function(x) x[length(x)], ""))

rex <- read_tsv("data/gag_full/gag_full_plaac.tsv", show_col_types = FALSE) %>%
  select(id = rexdb_id, superfamily, gag_len, has_prd) %>%
  inner_join(read_tsv("data/gag_full/nc_hits.tsv", show_col_types = FALSE) %>%
               select(id = rexdb_id, nc_combined), by = "id") %>%
  left_join(lineage, by = "id") %>%
  mutate(len = as.numeric(gag_len), prd = has_prd == 1, knuckle = nc_combined == 1)

TRUE_KNUCKLE <- c("PF00098", "SM00343", "PS50158", "SSF57756")
ipr <- lapply(list.files("data/gag_plaac/interpro", full.names = TRUE), function(f)
  read_tsv(f, col_names = FALSE, show_col_types = FALSE,
           col_types = readr::cols(.default = "c"))) %>% bind_rows()
ipr_k <- unique(ipr$X1[ipr$X5 %in% TRUE_KNUCKLE])

gydb <- read_tsv("data/gag_plaac/traits_272.tsv", show_col_types = FALSE) %>%
  transmute(id = tip, superfamily, len = as.numeric(gag_len),
            prd = has_prd == 1, knuckle = tip %in% ipr_k)


fisher_by_band <- function(counts) {
  counts %>%
    select(grp, band, prd, n, k) %>%
    pivot_wider(names_from = prd, values_from = c(n, k)) %>%
    rowwise() %>%
    mutate(p_value = fisher.test(matrix(c(k_TRUE, n_TRUE - k_TRUE,
                                          k_FALSE, n_FALSE - k_FALSE),
                                        nrow = 2))$p.value) %>%
    ungroup() %>%
    mutate(p_holm = p.adjust(p_value, "holm"),
           # stars are read off the Holm-adjusted value, not the raw one: these
           # figures make 6 and 11 comparisons, and marking raw significance
           # would promote bands that do not survive the correction
           label = cut(p_holm, c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("***", "**", "*", "ns")))
}


banded <- function(df, group_col) {
  df %>%
    mutate(band = cut(len, BREAKS, labels = BAND_LABELS, right = TRUE)) %>%
    filter(!is.na(band)) %>%
    group_by(grp = .data[[group_col]], band, prd) %>%
    summarise(n = n(), k = sum(knuckle), .groups = "drop") %>%
    mutate(p = k / n,
           lo = pmax(0, p - 1.96 * sqrt(p * (1 - p) / n)),
           hi = pmin(1, p + 1.96 * sqrt(p * (1 - p) / n))) %>%
    group_by(grp, band) %>%
    filter(n() == 2, min(n) >= MIN_CELL) %>%   # both groups present and populated
    ungroup() %>%
    mutate(status = factor(ifelse(prd, "carries a prion-like domain",
                                  "no prion-like domain"),
                           levels = c("carries a prion-like domain",
                                      "no prion-like domain")))
}


draw <- function(data, title, subtitle, ncol) {
  # only the bands that survived the n filter are drawn, in their real order --
  # ggplot will otherwise take whatever order the summarise happened to produce
  present <- BAND_LABELS[BAND_LABELS %in% unique(as.character(data$band))]
  # x is a factor here, so the note sits on a band rather than at a coordinate
  flat <- data %>% group_by(grp) %>% summarise(all_zero = sum(k) == 0,
                                               .groups = "drop") %>%
    filter(all_zero) %>%
    mutate(band = factor(present[ceiling(length(present) / 2)],
                         levels = BAND_LABELS), p = 0.5)
  tests <- fisher_by_band(data) %>%
    left_join(data %>% group_by(grp, band) %>%
                summarise(top = max(hi), .groups = "drop"),
              by = c("grp", "band"))
  ggplot(data, aes(band, p, fill = status)) +
    {if (nrow(flat)) geom_text(data = flat,
                               aes(x = band, y = p,
                                   label = "no zinc knuckle
anywhere in this lineage"),
                               inherit.aes = FALSE, size = 3, colour = INK2,
                               lineheight = 0.95, fontface = "italic")} +
    geom_col(position = position_dodge(width = 0.72), width = 0.64) +
    geom_errorbar(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.72), width = 0.12,
                  linewidth = 0.35, colour = INK2) +
    geom_text(aes(y = hi + 0.018, label = n),
              position = position_dodge(width = 0.72), size = 2.3, colour = INK2) +
    geom_text(data = tests,
              aes(x = band, y = pmin(top + 0.085, 1.06), label = label),
              inherit.aes = FALSE,
              size = ifelse(tests$p_holm < 0.05, 4.4, 2.7),
              vjust = ifelse(tests$p_holm < 0.05, 0.65, 0),
              colour = ifelse(tests$p_holm < 0.05, INK, INK2),
              fontface = ifelse(tests$p_holm < 0.05, "bold", "italic")) +
    facet_wrap(~ grp, ncol = ncol) +
    scale_fill_manual(values = c("carries a prion-like domain" = WITH_PRLD,
                                 "no prion-like domain" = NO_PRLD), name = NULL) +
    scale_x_discrete(limits = present) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.12),
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "Gag length (aa)", y = "share carrying a zinc knuckle",
         caption = paste("Fisher's exact test per band, Holm-adjusted:",
                         "*** p < 0.001, ** p < 0.01, * p < 0.05, ns not significant")) +
    theme_minimal(base_size = 10) +
    theme(text = element_text(colour = INK),
          strip.text = element_text(face = "bold", size = 9, colour = INK2),
          axis.title = element_text(size = 8.6, colour = INK2),
          axis.text = element_text(size = 8, colour = INK2),
          panel.grid = element_blank(),
          axis.line.y = element_line(colour = GUIDE, linewidth = 0.4),
          axis.ticks.y = element_line(colour = GUIDE, linewidth = 0.4),
          axis.ticks.length = unit(2.5, "pt"),
          panel.spacing = unit(10, "pt"),
          legend.position = "bottom", legend.text = element_text(size = 8.8),
          plot.caption = element_text(colour = INK2, size = 7.6, hjust = 0,
                                      margin = margin(t = 8)),
          plot.background = element_rect(fill = SURFACE, colour = NA),
          plot.margin = margin(10, 12, 6, 10))
}


tally <- function(data, what) {
  out <- fisher_by_band(data) %>%
    transmute(grp, band,
              PrLD_pos = sprintf("%d/%d", k_TRUE, n_TRUE),
              pct_pos = round(100 * k_TRUE / n_TRUE, 1),
              PrLD_neg = sprintf("%d/%d", k_FALSE, n_FALSE),
              pct_neg = round(100 * k_FALSE / n_FALSE, 1),
              fisher_p = signif(p_value, 3), holm_p = signif(p_holm, 3),
              direction = ifelse(k_TRUE / n_TRUE < k_FALSE / n_FALSE,
                                 "fewer in PrLD+", "more in PrLD+")) %>%
    arrange(grp, band)
  cat(sprintf(paste("%s: %d comparable bands; knuckle share lower in the PrLD+",
                    "group in %d; significant after Holm in %d\n"),
              what, nrow(out), sum(out$pct_pos < out$pct_neg),
              sum(out$holm_p < 0.05)))
  print(as.data.frame(out))
  cat("\n")
}

# ---- figure 1: by database ---------------------------------------------------

db <- bind_rows(
  rex %>% filter(superfamily == "Ty1/copia") %>% mutate(set = "REXdb  Ty1/copia"),
  rex %>% filter(superfamily == "Ty3/gypsy") %>% mutate(set = "REXdb  Ty3/gypsy"),
  gydb %>% mutate(set = "GyDB  all superfamilies")) %>%
  mutate(set = factor(set, levels = c("REXdb  Ty1/copia", "REXdb  Ty3/gypsy",
                                      "GyDB  all superfamilies")))

db_bands <- banded(db, "set")
fig1 <- draw(db_bands,
             "Zinc knuckle and prion-like domain, compared at matched Gag length",
             paste("share of Gags carrying a zinc knuckle; 95% intervals,",
                   "n above each bar.
GyDB is absent: 25 domains in total",
                   "cannot fill a band at n >= 20 in both groups."),
             ncol = n_distinct(db_bands$grp))
ggsave(file.path(OUT, "fig_knuckle_matched_db.png"), fig1, width = 7.8,
       height = 5, dpi = 400, bg = SURFACE)
ggsave(file.path(OUT, "fig_knuckle_matched_db.pdf"), fig1, width = 7.8,
       height = 5, bg = SURFACE)

# ---- figure 2: by superfamily and lineage ------------------------------------

lin <- bind_rows(
  rex %>% filter(!is.na(lin), lin != "") %>%
    mutate(set = paste0(ifelse(superfamily == "Ty1/copia", "copia  ", "gypsy  "), lin)),
  gydb %>% mutate(set = paste0("GyDB  ", superfamily)))

lin_bands <- banded(lin, "set")
keep <- lin_bands %>% count(grp) %>% filter(n >= 4) %>% pull(grp)  # >=2 bands
lin_bands <- filter(lin_bands, grp %in% keep)

fig2 <- draw(lin_bands,
             "The same comparison inside each lineage",
             paste("every lineage with two or more comparable length windows;",
                   "empty space is where one group has too few Gags of that length"),
             ncol = 4)
ggsave(file.path(OUT, "fig_knuckle_matched_lineage.png"), fig2, width = 11.5,
       height = 6.6, dpi = 400, bg = SURFACE)
ggsave(file.path(OUT, "fig_knuckle_matched_lineage.pdf"), fig2, width = 11.5,
       height = 6.6, bg = SURFACE)

cat("wrote fig_knuckle_matched_db and fig_knuckle_matched_lineage\n\n")
tally(db_bands, "by database")
tally(lin_bands, "by lineage")
