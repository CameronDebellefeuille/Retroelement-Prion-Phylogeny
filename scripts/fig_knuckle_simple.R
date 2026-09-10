# The anti-correlation stated as plainly as it can be stated.
#
#   fig_knuckle_simple.{png,pdf}
#
# Panel A: of the Gags that carry a prion-like domain, what share carry a zinc
# knuckle -- against the same share among Gags with no domain. Two numbers per
# dataset, with 95% intervals. This is the simplest honest form of the claim.
#
# Panel B: the same pair of numbers computed inside each lineage separately, one
# point per lineage, joined. A boxplot needs a distribution, and two percentages
# do not have one -- but across lineages they do, and that distribution is the
# more useful picture: it shows whether the gap is a property of Gag or an
# accident of which clades were sampled. Each lineage contributes one paired
# observation, which is also a crude guard against the pseudoreplication that
# affects every element-level test here.
#
# What panel A cannot do is control for length, and length predicts a domain far
# more strongly than the knuckle does. Gags with a domain are longer, and longer
# Gags are likelier to carry anything at all -- so the raw gap in panel A is a
# floor on the effect, not a measurement of it. fig_knuckle_prld.R has the
# length-matched version.
#
# Run from the repo root, in the conda env:
#   Rscript scripts/fig_knuckle_simple.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

OUT <- "data/gag_plaac/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"; GUIDE <- "#e6e4df"
WITH_PRLD <- "#7B2D9E"   # Gags that carry a prion-like domain
NO_PRLD <- "#8a8a85"     # Gags that do not

MIN_GROUP <- 15          # a lineage needs this many in both groups to be plotted

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
  transmute(id = tip, superfamily, lin = superfamily, len = as.numeric(gag_len),
            prd = has_prd == 1, knuckle = tip %in% ipr_k)

share <- function(df, label) {
  df %>% group_by(prd) %>%
    summarise(n = n(), k = sum(knuckle), .groups = "drop") %>%
    mutate(set = label, p = k / n,
           lo = pmax(0, p - 1.96 * sqrt(p * (1 - p) / n)),
           hi = pmin(1, p + 1.96 * sqrt(p * (1 - p) / n)))
}

bars <- bind_rows(
  share(filter(rex, superfamily == "Ty1/copia"), "REXdb\nTy1/copia"),
  share(filter(rex, superfamily == "Ty3/gypsy"), "REXdb\nTy3/gypsy"),
  share(gydb, "GyDB\nall")) %>%
  mutate(group = ifelse(prd, "carries a prion-like domain", "no prion-like domain"),
         group = factor(group, levels = c("carries a prion-like domain",
                                          "no prion-like domain")),
         set = factor(set, levels = c("REXdb\nTy1/copia", "REXdb\nTy3/gypsy",
                                      "GyDB\nall")))

panel_a <- ggplot(bars, aes(set, p, fill = group)) +
  geom_col(position = position_dodge(width = 0.68), width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.68),
                width = 0.13, linewidth = 0.4, colour = INK2) +
  geom_text(aes(y = hi + 0.035, label = sprintf("%.0f%%\n(%d/%d)", 100 * p, k, n)),
            position = position_dodge(width = 0.68), size = 2.9, colour = INK2,
            lineheight = 0.95) +
  scale_fill_manual(values = c("carries a prion-like domain" = WITH_PRLD,
                               "no prion-like domain" = NO_PRLD), name = NULL) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "share carrying a zinc knuckle",
       title = "A   Gags with a prion-like domain carry a knuckle less often",
       subtitle = "raw shares, not adjusted for Gag length") +
  theme_minimal(base_size = 10) +
  theme(text = element_text(colour = INK),
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(colour = INK2, size = 8.8, margin = margin(b = 8)),
        axis.text.x = element_text(size = 9, colour = INK, lineheight = 0.95),
        axis.text.y = element_text(size = 8, colour = INK2),
        axis.title.y = element_text(size = 8.4, colour = INK2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = GUIDE, linewidth = 0.3),
        legend.position = "bottom", legend.text = element_text(size = 8.6),
        plot.background = element_rect(fill = SURFACE, colour = NA),
        plot.margin = margin(8, 10, 4, 8))

# ---- panel B: one paired observation per lineage ------------------------------

per_lineage <- bind_rows(
  rex %>% filter(!is.na(lin), lin != "") %>% mutate(src = "REXdb"),
  gydb %>% mutate(src = "GyDB")) %>%
  group_by(src, lin, prd) %>%
  summarise(n = n(), k = sum(knuckle), .groups = "drop") %>%
  filter(n >= MIN_GROUP) %>%
  mutate(p = k / n) %>%
  select(src, lin, prd, p, n) %>%
  pivot_wider(names_from = prd, values_from = c(p, n),
              names_prefix = "prd") %>%
  filter(!is.na(p_prdTRUE), !is.na(p_prdFALSE))

paired <- per_lineage %>%
  pivot_longer(c(p_prdTRUE, p_prdFALSE), names_to = "grp", values_to = "p") %>%
  mutate(group = ifelse(grp == "p_prdTRUE", "carries a prion-like domain",
                        "no prion-like domain"),
         group = factor(group, levels = c("carries a prion-like domain",
                                          "no prion-like domain")))

wilcox <- suppressWarnings(
  wilcox.test(per_lineage$p_prdTRUE, per_lineage$p_prdFALSE, paired = TRUE))
lower <- sum(per_lineage$p_prdTRUE < per_lineage$p_prdFALSE)

panel_b <- ggplot(paired, aes(group, p)) +
  geom_boxplot(aes(fill = group), width = 0.45, outlier.shape = NA,
               colour = INK2, linewidth = 0.4, alpha = 0.55) +
  geom_line(aes(group = paste(src, lin)), colour = "#b9b8b3", linewidth = 0.4) +
  geom_point(aes(colour = group), size = 2.2) +
  geom_text(data = filter(paired, group == "no prion-like domain"),
            aes(label = lin), hjust = -0.22, size = 2.7, colour = INK2) +
  scale_fill_manual(values = c("carries a prion-like domain" = WITH_PRLD,
                               "no prion-like domain" = NO_PRLD), guide = "none") +
  scale_colour_manual(values = c("carries a prion-like domain" = WITH_PRLD,
                                 "no prion-like domain" = NO_PRLD), guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(expand = expansion(add = c(0.6, 0.95))) +
  labs(x = NULL, y = "share carrying a zinc knuckle",
       title = "B   The same comparison made inside each lineage",
       subtitle = sprintf(paste("one point per lineage with at least %d Gags in",
                                "both groups; lower in %d of %d, paired",
                                "Wilcoxon p = %.3f"),
                          MIN_GROUP, lower, nrow(per_lineage), wilcox$p.value)) +
  theme_minimal(base_size = 10) +
  theme(text = element_text(colour = INK),
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(colour = INK2, size = 8.8, margin = margin(b = 8)),
        axis.text.x = element_text(size = 9, colour = INK),
        axis.text.y = element_text(size = 8, colour = INK2),
        axis.title.y = element_text(size = 8.4, colour = INK2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = GUIDE, linewidth = 0.3),
        plot.background = element_rect(fill = SURFACE, colour = NA),
        plot.margin = margin(8, 10, 4, 8))

figure <- panel_a | panel_b

ggsave(file.path(OUT, "fig_knuckle_simple.png"), figure, width = 10, height = 5.2,
       dpi = 400, bg = SURFACE)
ggsave(file.path(OUT, "fig_knuckle_simple.pdf"), figure, width = 10, height = 5.2,
       bg = SURFACE)
cat("wrote fig_knuckle_simple\n\n")
print(as.data.frame(bars %>% select(set, group, k, n, p)), digits = 3)
cat("\nper lineage:\n")
print(as.data.frame(per_lineage), digits = 3)
