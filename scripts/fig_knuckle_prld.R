# The zinc knuckle / prion-like domain anti-correlation, shown honestly.
#
#   fig_knuckle_prld.{png,pdf}
#
# The difficulty with this result is that Gag length predicts a prion-like domain
# far more strongly than the knuckle does -- roughly a hundredfold per e-fold of
# length -- so a bare 2x2 shows something real for the wrong reason. Both panels
# are built to make that visible rather than to hide it.
#
# Panel A plots the domain rate against length, separately for Gags with and
# without a knuckle. The steep rise along x IS the length effect; the vertical
# gap between the two curves at the same x is the knuckle effect, which is the
# claim. If the curves overlapped, there would be nothing here.
#
# Panel B is the same comparison as an odds ratio adjusted for log(length), in
# each stratum that can carry the test on its own: two REXdb lineages with enough
# domains, the REXdb superfamilies, and GyDB on its InterProScan calls. Strata
# are independent samples of the same question, so consistency across them is the
# argument -- not any single p-value.
#
# Knuckle calls: REXdb uses `nc_combined` (Pfam E<1e-3 or PROSITE PS50158, 92%
# sensitive and 99% specific against InterProScan); GyDB uses InterProScan itself.
#
# Neither panel corrects for phylogeny. Elements within a lineage are related, so
# the intervals are narrower than the evidence warrants. That is stated on the
# figure rather than left to the caption.
#
# Run from the repo root, in the conda env:
#   Rscript scripts/fig_knuckle_prld.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

OUT <- "data/gag_plaac/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"; GUIDE <- "#e6e4df"
NO_K <- "#D2551B"    # no knuckle
HAS_K <- "#2D7DD2"   # knuckle present

# ---- REXdb ------------------------------------------------------------------

# the classification has a ragged number of columns -- the lineage is whatever
# the last field on the row is, so it is read line by line rather than as a table
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

# ---- GyDB, on the InterProScan calls -----------------------------------------

TRUE_KNUCKLE <- c("PF00098", "SM00343", "PS50158", "SSF57756")
ipr <- lapply(list.files("data/gag_plaac/interpro", full.names = TRUE), function(f)
  read_tsv(f, col_names = FALSE, show_col_types = FALSE,
           col_types = readr::cols(.default = "c"))) %>%
  bind_rows()
ipr_k <- unique(ipr$X1[ipr$X5 %in% TRUE_KNUCKLE])

gydb <- read_tsv("data/gag_plaac/traits_272.tsv", show_col_types = FALSE) %>%
  transmute(id = tip, superfamily, len = as.numeric(gag_len),
            prd = has_prd == 1, knuckle = tip %in% ipr_k)

# ---- panel A: rate against length --------------------------------------------

strata <- bind_rows(
  rex %>% filter(superfamily == "Ty1/copia") %>% mutate(panel = "REXdb  Ty1/copia"),
  rex %>% filter(lin == "Ale") %>% mutate(panel = "REXdb  Ale lineage"),
  rex %>% filter(lin == "Athila") %>% mutate(panel = "REXdb  Athila lineage"),
  gydb %>% mutate(panel = "GyDB  all superfamilies")) %>%
  mutate(panel = factor(panel, levels = c("REXdb  Ty1/copia", "REXdb  Ale lineage",
                                          "REXdb  Athila lineage",
                                          "GyDB  all superfamilies")))

binned <- strata %>%
  group_by(panel) %>%
  mutate(bin = cut(len, breaks = unique(quantile(len, seq(0, 1, 0.2), na.rm = TRUE)),
                   include.lowest = TRUE)) %>%
  group_by(panel, bin, knuckle) %>%
  summarise(n = n(), hits = sum(prd), mid = median(len), .groups = "drop") %>%
  filter(n >= 25) %>%
  mutate(rate = hits / n,
         lo = pmax(0, rate - 1.96 * sqrt(rate * (1 - rate) / n)),
         hi = pmin(1, rate + 1.96 * sqrt(rate * (1 - rate) / n)),
         status = ifelse(knuckle, "zinc knuckle present", "no knuckle"))

panel_a <- ggplot(binned, aes(mid, rate, colour = status, fill = status)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.13, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(size = n), shape = 21, colour = "white", stroke = 0.7) +
  scale_colour_manual(values = c("no knuckle" = NO_K,
                                 "zinc knuckle present" = HAS_K), name = NULL) +
  scale_fill_manual(values = c("no knuckle" = NO_K,
                               "zinc knuckle present" = HAS_K), name = NULL) +
  scale_size_area(max_size = 5, guide = "none") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = scales::pretty_breaks(4)) +
  facet_wrap(~ panel, scales = "free", nrow = 1) +
  labs(x = "Gag length (aa), binned by quintile", y = "elements with a prion-like domain",
       title = "A   At matched Gag length, knuckle-bearing Gags carry fewer prion-like domains",
       subtitle = paste("the rise along x is the length effect; the gap between the",
                        "curves at the same x is the knuckle effect")) +
  theme_minimal(base_size = 10) +
  theme(text = element_text(colour = INK),
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(colour = INK2, size = 8.8,
                                     margin = margin(b = 8)),
        strip.text = element_text(face = "bold", size = 8.6, colour = INK2),
        axis.title = element_text(size = 8.4, colour = INK2),
        axis.text = element_text(size = 7.8, colour = INK2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = GUIDE, linewidth = 0.3),
        legend.position = "bottom", legend.text = element_text(size = 8.6),
        plot.background = element_rect(fill = SURFACE, colour = NA),
        plot.margin = margin(8, 10, 4, 8))

# ---- panel B: adjusted odds ratios, one per independent stratum ---------------

fit_or <- function(df, label) {
  df <- filter(df, !is.na(len), len > 0)
  if (sum(df$prd) < 5 || length(unique(df$knuckle)) < 2) return(NULL)
  m <- suppressWarnings(glm(prd ~ knuckle + log(len), binomial, df))
  s <- summary(m)$coefficients
  if (!"knuckleTRUE" %in% rownames(s)) return(NULL)
  b <- s["knuckleTRUE", ]
  tibble(label = label, or = exp(b[1]),
         lo = exp(b[1] - 1.96 * b[2]), hi = exp(b[1] + 1.96 * b[2]),
         p = b[4], n = nrow(df), events = sum(df$prd))
}

forest <- bind_rows(
  fit_or(filter(rex, superfamily == "Ty1/copia"), "REXdb  Ty1/copia, all"),
  fit_or(filter(rex, lin == "Ale"), "REXdb  Ale lineage"),
  fit_or(filter(rex, lin == "Athila"), "REXdb  Athila lineage"),
  fit_or(filter(rex, superfamily == "Ty3/gypsy",
                !lin %in% c("Ogre", "TatI", "TatII", "TatIII", "Tatius", "Retand")),
         "REXdb  Ty3/gypsy, no Tat/Retand"),
  fit_or(gydb, "GyDB  all superfamilies"),
  fit_or(filter(gydb, superfamily == "Ty3/Gypsy"), "GyDB  Ty3/Gypsy"),
  fit_or(filter(gydb, superfamily == "Retroviridae"), "GyDB  Retroviridae")) %>%
  mutate(label = factor(label, levels = rev(label)),
         favours = ifelse(or < 1, "fewer domains when a knuckle is present", "the reverse"))

panel_b <- ggplot(forest, aes(or, label)) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.4, linetype = "22") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.7,
                 colour = HAS_K) +
  geom_point(size = 2.6, colour = HAS_K) +
  geom_text(aes(x = 14, label = sprintf("%d domains / %d", events, n)),
            hjust = 0, size = 2.8, colour = INK2) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 0.5, 1, 2, 4),
                labels = c("0.001", "0.01", "0.1", "0.5", "1", "2", "4"),
                limits = c(0.0008, 60)) +
  annotate("text", x = 0.0009, y = length(levels(forest$label)) + 0.45,
           label = "fewer domains when a knuckle is present", hjust = 0,
           size = 2.9, colour = INK2, fontface = "italic") +
  annotate("text", x = 1.35, y = length(levels(forest$label)) + 0.45,
           label = "the reverse", hjust = 0, size = 2.9, colour = INK2,
           fontface = "italic") +
  coord_cartesian(clip = "off") +
  labs(x = "odds ratio for a prion-like domain, adjusted for log(Gag length)",
       y = NULL,
       title = "B   The same effect in each stratum that can be tested on its own",
       subtitle = paste("independent samples of one question; consistency of",
                        "direction is the argument, not any single interval")) +
  theme_minimal(base_size = 10) +
  theme(text = element_text(colour = INK),
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(colour = INK2, size = 8.8,
                                     margin = margin(b = 8)),
        axis.title.x = element_text(size = 8.4, colour = INK2,
                                    margin = margin(t = 6)),
        axis.text.y = element_text(size = 8.6, colour = INK),
        axis.text.x = element_text(size = 7.8, colour = INK2),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = GUIDE, linewidth = 0.3),
        plot.background = element_rect(fill = SURFACE, colour = NA),
        plot.margin = margin(8, 10, 4, 8))

caption <- paste("Knuckle calls: REXdb by Pfam E<1e-3 or PROSITE PS50158",
                 "(92% sensitivity, 99% specificity against InterProScan);",
                 "GyDB by InterProScan.\nNo model corrects for shared ancestry:",
                 "elements within a lineage are related, so the intervals are",
                 "narrower than the evidence warrants.")

figure <- panel_a / panel_b +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(caption = caption,
                  theme = theme(plot.caption = element_text(hjust = 0, size = 7.6,
                                                            colour = INK2,
                                                            margin = margin(t = 8)),
                                plot.background = element_rect(fill = SURFACE,
                                                               colour = NA)))

ggsave(file.path(OUT, "fig_knuckle_prld.png"), figure, width = 11, height = 8.4,
       dpi = 400, bg = SURFACE)
ggsave(file.path(OUT, "fig_knuckle_prld.pdf"), figure, width = 11, height = 8.4,
       bg = SURFACE)
cat("wrote fig_knuckle_prld\n")
print(as.data.frame(forest %>% select(label, or, lo, hi, p, events, n)), digits = 3)
