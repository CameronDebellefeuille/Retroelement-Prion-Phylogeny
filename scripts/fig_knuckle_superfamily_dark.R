# fig_knuckle_superfamily_dark.{png,pdf} -- knuckle rate by PrLD status, for a
# dark slide.
#
# The slide version of scripts/fig_knuckle_superfamily.R: same data, same
# tests, same geometry, repainted light-on-dark. Unlike the other dark slide
# figures this one is drawn on SOLID BLACK rather than transparent, by request.
# The bar borders are painted in the background colour too -- they exist to cut
# a gap between the dodged bars, so they have to be whatever is behind them.
#
# The caption is dropped, matching how the figure is cropped in the deck. Read
# the warning under "significance" below before showing this without it.
#
# Run from the repo root, in the conda env:
#   Rscript scripts/fig_knuckle_superfamily_dark.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

OUT <- "figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

BG <- "#000000"
INK <- "#ffffff"        # percent labels, class names, stars
INK2 <- "#c9c7d6"       # tick labels, error bars, denominators, brackets
GUIDE <- "#6f6c82"      # axis line and ticks

# Purple kept, but lifted well off #7B2D9E, which goes to mud on black.
# Deliberately not the #C58BF5 the other slides use for Retroviridae -- that
# hue already means a superfamily in this deck, and Retroviridae is one of the
# groups on this very axis.
WITH_PRLD <- "#b45cf0"
NO_PRLD   <- "#9b98a8"

TRUE_KNUCKLE <- c("PF00098", "SM00343", "PS50158", "SSF57756")
KEEP <- c("Ty1/Copia", "Ty3/Gypsy", "Retroviridae")

ipr <- lapply(list.files("data/processed/gydb/InterProScan", pattern = "^part",
                         full.names = TRUE),
              read_tsv, col_names = FALSE, show_col_types = FALSE,
              progress = FALSE) %>% bind_rows()

scanned <- unique(ipr$X1)
knuckled <- unique(ipr$X1[ipr$X5 %in% TRUE_KNUCKLE])

traits <- read_tsv("data/processed/gydb/traits_272.tsv", show_col_types = FALSE) %>%
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
# WARNING, and the reason the source figure carries a caption: the stars are
# RAW Fisher's exact p. Three tests are run on one question, and once Holm
# correction is applied NEITHER nominal result survives -- both go to 0.077.
# This slide version has no caption to say so, so the caveat has to be carried
# by whoever is presenting it.
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

# The bracket has to clear the PERCENT LABELS, not the error bars.
PCT_OFFSET <- 0.045

brackets <- bars %>% group_by(superfamily) %>%
  summarise(top = max(hi), .groups = "drop") %>%
  left_join(tests, by = "superfamily") %>%
  mutate(x = as.integer(superfamily), y = top + PCT_OFFSET + 0.065)

dodge <- position_dodge(width = 0.68)

plot <- ggplot(bars, aes(superfamily, pct, fill = group)) +
  # No bar border. The paper figure strokes each bar in the background colour
  # to cut a gap between the dodged pair, but the dodge already leaves one, and
  # a border painted in the backing colour cannot survive being made
  # transparent. Dropping it lets one spec serve both backings.
  geom_col(position = dodge, width = 0.60, colour = NA) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = dodge, width = 0.14,
                linewidth = 0.5, colour = INK2) +
  geom_text(aes(y = hi + PCT_OFFSET, label = sprintf("%.0f%%", 100 * pct)),
            position = dodge, size = 4.4, colour = INK, fontface = "bold") +
  geom_text(aes(y = -0.055, label = sprintf("%d/%d", k, n)),
            position = dodge, size = 3.6, colour = INK2) +
  geom_segment(data = brackets, aes(x = x - 0.17, xend = x + 0.17,
                                    y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.45, colour = INK2) +
  geom_segment(data = brackets, aes(x = x - 0.17, xend = x - 0.17,
                                    y = y, yend = y - 0.022),
               inherit.aes = FALSE, linewidth = 0.45, colour = INK2) +
  geom_segment(data = brackets, aes(x = x + 0.17, xend = x + 0.17,
                                    y = y, yend = y - 0.022),
               inherit.aes = FALSE, linewidth = 0.45, colour = INK2) +
  geom_text(data = brackets, aes(x = x, y = y + 0.014, label = star),
            inherit.aes = FALSE, size = 5.6, colour = INK, fontface = "bold",
            vjust = 0) +
  scale_fill_manual(values = c("PrLD present" = WITH_PRLD,
                               "PrLD absent" = NO_PRLD), name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, 0.25), limits = c(-0.09, 1.20),
                     expand = c(0, 0)) +
  labs(y = "Gags carrying a zinc knuckle", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    text                = element_text(colour = INK),
    axis.text.x         = element_text(size = 13, face = "bold", colour = INK),
    axis.text.y         = element_text(size = 11, colour = INK2),
    axis.title.y        = element_text(size = 12, colour = INK,
                                       margin = margin(r = 8)),
    panel.grid          = element_blank(),
    axis.line.y         = element_line(colour = GUIDE, linewidth = 0.5),
    axis.ticks.y        = element_line(colour = GUIDE, linewidth = 0.5),
    axis.ticks.length.y = unit(3.5, "pt"),
    legend.position     = "top",
    legend.text         = element_text(size = 12, colour = INK),
    legend.key.size     = unit(13, "pt"),
    legend.margin       = margin(b = -2),
    plot.margin         = margin(12, 16, 8, 10))

# Two backings off the one spec: solid black, and the same figure with every
# pane knocked out for slides that would rather not carry a black rectangle.
# Only the backing differs -- identical data, geometry, type and colour.
#
# type = "cairo" for the anti-aliasing, as with the other dark figures; the
# Windows GDI default would quantise this to a handful of flat colours and,
# on the transparent build, throw away the alpha channel entirely.
render <- function(name, bg) {
  pane <- if (identical(bg, "transparent")) NA else bg
  p <- plot + theme(
    plot.background   = element_rect(fill = pane, colour = NA),
    panel.background  = element_rect(fill = pane, colour = NA),
    legend.background = element_rect(fill = pane, colour = NA),
    legend.key        = element_rect(fill = pane, colour = NA))
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = 8.4, height = 5.4,
         dpi = 400, bg = bg, type = "cairo")
  cat("wrote", name, "\n")
}

render("fig_knuckle_superfamily_dark_transparent", "transparent")

bars %>% mutate(across(c(pct, lo, hi), ~ sprintf("%.0f%%", 100 * .x))) %>%
  select(superfamily, group, k, n, pct, lo, hi) %>% print(n = 20)
tests %>% select(superfamily, p, p_holm, star) %>%
  mutate(across(c(p, p_holm), ~ round(.x, 4))) %>% print()
