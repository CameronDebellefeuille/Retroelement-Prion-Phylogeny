# Zinc knuckle presence, split by whether the Gag carries a prion-like domain.
#
#   fig_knuckle_superfamily.{png,pdf}
#
# The companion to the GyDB RT tree slides: the tree shows WHERE the domains
# are, this shows what travels with them. Within each superfamily, the share of
# Gags carrying a knuckle among those with a PrLD against those without.
#
# Three superfamilies, not four. Bel/Pao has ZERO PrLD-bearing Gags (0/23), so
# it has no left-hand bar to draw -- it is named in the caption instead of being
# given an empty slot that would read as 0%.
#
# Knuckle calls are InterProScan domain-level signatures of IPR001878/IPR036875
# -- PF00098, SM00343, PS50158, SSF57756 -- the definition GYDB_FINDINGS 3.5
# revision 2 settled on. The regex over-calls (72% specific) and the Pfam
# gathering threshold under-calls (61% sensitive); revision 1 mistook the
# resulting attenuation for absence of effect, which is why the definition is
# named on the figure itself.
#
# 95% Wilson intervals are drawn because the PrLD+ groups are n = 3, 17 and 5.
# A bare bar at 0% for three copia elements would imply a precision that does
# not exist; the interval reaching 56% is the honest statement.
#
# Run from the repo root, in the conda env:
#   Rscript scripts/fig_knuckle_superfamily.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

OUT <- "data/gag_plaac/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"; GUIDE <- "#e6e4df"
WITH_PRLD <- "#7B2D9E"
NO_PRLD   <- "#8a8a85"

TRUE_KNUCKLE <- c("PF00098", "SM00343", "PS50158", "SSF57756")
KEEP <- c("Ty1/Copia", "Ty3/Gypsy", "Retroviridae")

# ---- InterProScan calls -----------------------------------------------------
ipr <- lapply(list.files("data/gag_plaac/interpro", pattern = "^part",
                         full.names = TRUE),
              read_tsv, col_names = FALSE, show_col_types = FALSE,
              progress = FALSE) %>% bind_rows()

scanned <- unique(ipr$X1)
knuckled <- unique(ipr$X1[ipr$X5 %in% TRUE_KNUCKLE])

traits <- read_tsv("data/gag_plaac/traits_272.tsv", show_col_types = FALSE) %>%
  filter(element %in% scanned, superfamily %in% KEEP) %>%
  mutate(group = ifelse(has_prd == 1, "PrLD present", "PrLD absent"),
         knuckle = element %in% knuckled)

wilson <- function(k, n, z = 1.96) {
  p <- k / n; d <- 1 + z^2 / n
  centre <- p + z^2 / (2 * n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(lo = (centre - half) / d, hi = (centre + half) / d)
}

bars <- traits %>%
  group_by(superfamily, group) %>%
  summarise(n = n(), k = sum(knuckle), .groups = "drop") %>%
  mutate(pct = k / n,
         lo = mapply(function(k, n) wilson(k, n)["lo"], k, n),
         hi = mapply(function(k, n) wilson(k, n)["hi"], k, n),
         superfamily = factor(superfamily, levels = KEEP),
         group = factor(group, levels = c("PrLD present", "PrLD absent")))

# ---- significance -----------------------------------------------------------
# Fisher's exact, not a t-test: the outcome is binary (knuckle yes/no) and the
# PrLD+ cells are n = 3, 17 and 5, where a test on means is not defined. This is
# the test GYDB_FINDINGS 3.5 uses, so the p-values here are comparable to it.
#
# Stars are RAW p. Holm-adjusted values are printed in the caption because three
# tests are run on one question, and NEITHER nominal result survives correction
# (both go to 0.077). Showing only the stars would overstate the evidence;
# showing only Holm would hide a consistent direction across three groups.
tests <- bars %>%
  select(superfamily, group, k, n) %>%
  mutate(no_k = n - k) %>%
  select(-n) %>%
  pivot_wider(names_from = group, values_from = c(k, no_k)) %>%
  rowwise() %>%
  mutate(p = fisher.test(matrix(c(`k_PrLD present`, `no_k_PrLD present`,
                                  `k_PrLD absent`, `no_k_PrLD absent`),
                                nrow = 2))$p.value) %>%
  ungroup() %>%
  mutate(p_holm = p.adjust(p, method = "holm"),
         star = ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "ns"))))

# The bracket has to clear the PERCENT LABELS, not the error bars -- the labels
# sit at hi + 0.045, so anchoring on hi alone puts the bracket through them.
PCT_OFFSET <- 0.045

brackets <- bars %>% group_by(superfamily) %>%
  summarise(top = max(hi), .groups = "drop") %>%
  left_join(tests, by = "superfamily") %>%
  mutate(x = as.integer(superfamily), y = top + PCT_OFFSET + 0.055)

dodge <- position_dodge(width = 0.68)

plot <- ggplot(bars, aes(superfamily, pct, fill = group)) +
  geom_col(position = dodge, width = 0.60, colour = SURFACE, linewidth = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = dodge, width = 0.14,
                linewidth = 0.4, colour = INK2) +
  # the rate on the bar, the denominator under it -- n = 3 and n = 5 are the
  # whole story here and must not be something the reader has to hunt for
  geom_text(aes(y = hi + PCT_OFFSET, label = sprintf("%.0f%%", 100 * pct)),
            position = dodge, size = 3.2, colour = INK, fontface = "bold") +
  geom_text(aes(y = -0.055, label = sprintf("%d/%d", k, n)),
            position = dodge, size = 2.9, colour = INK2) +
  geom_segment(data = brackets, aes(x = x - 0.17, xend = x + 0.17,
                                    y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.35, colour = INK2) +
  geom_segment(data = brackets, aes(x = x - 0.17, xend = x - 0.17,
                                    y = y, yend = y - 0.022),
               inherit.aes = FALSE, linewidth = 0.35, colour = INK2) +
  geom_segment(data = brackets, aes(x = x + 0.17, xend = x + 0.17,
                                    y = y, yend = y - 0.022),
               inherit.aes = FALSE, linewidth = 0.35, colour = INK2) +
  geom_text(data = brackets,
            aes(x = x, y = y + 0.012, label = star),
            inherit.aes = FALSE, size = 4.4, colour = INK, fontface = "bold",
            vjust = 0) +
  scale_fill_manual(values = c("PrLD present" = WITH_PRLD,
                               "PrLD absent" = NO_PRLD), name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, 0.25), limits = c(-0.09, 1.20),
                     expand = c(0, 0)) +
  labs(y = "Gags carrying a zinc knuckle",
       x = NULL,
       caption = paste0(
         "GyDB Gag cores with an InterProScan result. Knuckle = domain-level ",
         "IPR001878 / IPR036875
",
         "(PF00098, SM00343, PS50158, SSF57756). Bars are 95% Wilson intervals.
",
         "Bel/Pao is not shown: 0 of its 23 Gags carry a PrLD, so it has no ",
         "comparison to make.
",
         "Stars are Fisher's exact on raw p. Holm-adjusted across the three ",
         "tests, neither nominal result survives (both p = 0.077).")) +
  theme_minimal(base_size = 11) +
  theme(
    text             = element_text(colour = INK),
    axis.text.x      = element_text(size = 11, face = "bold", colour = INK),
    axis.text.y      = element_text(size = 9, colour = INK2),
    axis.title.y     = element_text(size = 9.4, colour = INK2,
                                    margin = margin(r = 6)),
    panel.grid       = element_blank(),
    axis.line.y      = element_line(colour = GUIDE, linewidth = 0.4),
    axis.ticks.y     = element_line(colour = GUIDE, linewidth = 0.4),
    axis.ticks.length.y = unit(3, "pt"),
    legend.position  = "top",
    legend.text      = element_text(size = 9.6),
    legend.margin    = margin(b = -4),
    plot.caption     = element_text(colour = INK2, size = 7.4, hjust = 0,
                                    margin = margin(t = 10)),
    plot.background  = element_rect(fill = SURFACE, colour = NA),
    plot.margin      = margin(10, 14, 8, 10))

ggsave(file.path(OUT, "fig_knuckle_superfamily.png"), plot, width = 7.2,
       height = 5, dpi = 400, bg = SURFACE)
ggsave(file.path(OUT, "fig_knuckle_superfamily.pdf"), plot, width = 7.2,
       height = 5, bg = SURFACE)

cat("wrote fig_knuckle_superfamily\n")
bars %>% mutate(across(c(pct, lo, hi), ~ sprintf("%.0f%%", 100 * .x))) %>%
  select(superfamily, group, k, n, pct, lo, hi) %>% print(n = 20)
tests %>% select(superfamily, p, p_holm, star) %>%
  mutate(across(c(p, p_holm), ~ round(.x, 4))) %>% print()
