# Disorder by protein class -- the companion to fig_gydb_classes (PLAAC LLR).
#   fig_gydb_classes_disorder.{png,pdf}
#
# Same layout as the LLR figure so the two can be read side by side. The point of
# the pair is the contrast: PLAAC calls are confined to Gag, disorder is not.
#
# Run from the repo root, in the conda env:
#   Rscript scripts/figures_disorder_classes.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

OUT <- "figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ACCENT <- "#2a78d6"; NEUTRAL <- "#9d9c96"
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"; SURFACE <- "#fcfcfb"
VIOLIN <- "#8fcbbc"; HITDOT <- "#0d366b"

theme_pub <- function() {
  theme_minimal(base_size = 11) +
    theme(
      text             = element_text(colour = INK),
      axis.text        = element_text(colour = INK2),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
      panel.background = element_rect(fill = SURFACE, colour = NA),
      plot.background  = element_rect(fill = SURFACE, colour = NA),
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(colour = INK2, size = 8.8),
      plot.caption     = element_text(colour = INK2, size = 7.6, hjust = 0),
      legend.position  = "top",
      legend.key.size  = unit(0.8, "lines"))
}

d <- read_tsv("data/processed/gydb/disorder_all.tsv", show_col_types = FALSE) %>%
  group_by(protein_class) %>% mutate(n_scored = n()) %>% ungroup() %>%
  filter(n_scored >= 20)

order_by_median <- d %>% group_by(protein_class) %>%
  summarise(m = median(frac_disordered), .groups = "drop") %>%
  arrange(m) %>% pull(protein_class)

d <- mutate(d, protein_class = factor(protein_class, levels = order_by_median))
hits <- filter(d, has_prd == 1)
counts <- count(d, protein_class, name = "n_scored")

p <- ggplot(d, aes(x = protein_class, y = frac_disordered)) +
  geom_violin(fill = VIOLIN, colour = "#3f8a78", alpha = 0.55,
              linewidth = 0.3, scale = "width", width = 0.9) +
  geom_boxplot(width = 0.16, outlier.shape = NA, colour = INK2,
               fill = "#00000018", linewidth = 0.3) +
  geom_jitter(width = 0.22, height = 0, size = 0.45, alpha = 0.32, colour = INK) +
  geom_point(data = hits, size = 1.9, shape = 21, fill = HITDOT, colour = INK,
             stroke = 0.3,
             position = position_jitter(width = 0.22, height = 0, seed = 1)) +
  geom_text(data = counts, aes(y = 1.07, label = paste0("n = ", n_scored)),
            size = 3, colour = INK2, fontface = "italic") +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.02, 0.02))) +
  theme_pub() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "#d5d4cf", fill = NA,
                                    linewidth = 0.4),
        axis.text.x = element_text(size = 9, colour = INK),
        legend.position = "none") +
  labs(x = NULL, y = "Fraction of core in a predicted IDR (metapredict)",
       title = "Disorder is concentrated in Gag; the enzymes are ordered",
       subtitle = paste0(
         "GAG has the highest median of any class (0.29); RT, INT, AP and RNaseH sit at zero. ",
         "Dark points are the ", nrow(hits), " PLAAC calls."))

# The caption was dropped from the figure by request. What it said is kept here,
# because the numbers qualify how the figure should be read:
#
#   Protein classes with at least 20 scoreable cores, one point per core.
#
#   146 of 2,636 cores are excluded -- 145 contain non-standard residues that
#   metapredict rejects, 1 is under 20 aa.
#
#   CHR appears here but not in the LLR figure: metapredict has no 60 aa window
#   requirement, so short cores can be scored.
#
#   The spikes at 1.00 are a short-core artifact. The 95 cores scored fully
#   disordered have median length 69 aa against 226 for the rest, and 28 of them
#   are CHR. Fraction disordered is unstable below ~100 aa -- read the medians,
#   not the extremes.

ggsave(file.path(OUT, "fig_gydb_classes_disorder.png"), p, width = 8.8,
       height = 6.0, dpi = 300, bg = SURFACE)
cat("wrote fig_gydb_classes_disorder\n")

d %>% group_by(protein_class) %>%
  summarise(n = n(), median = round(median(frac_disordered), 3),
            .groups = "drop") %>% arrange(desc(median)) %>% print(n = 20)
