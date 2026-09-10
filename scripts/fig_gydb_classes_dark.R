# fig_gydb_classes_dark.{png,pdf} -- LLR by protein class, for a dark slide.
#
# The slide version of fig_gydb_classes from figures_272.R (now in backup): same
# data, same geometry, repainted light-on-dark on a transparent background.
#
# Title, subtitle and caption are gone on purpose. The slide carries its own
# headline, and on the deck the whole block was being cropped away anyway.
# What is left is the panel, the y axis, the class names and the n per class.
#
# Palette matches scripts/fig_gydb_rt_layers.R so the PrLD dots are the same
# blue here as on the tree slides -- the marker means one thing in the deck.
#
# Run from the repo root:
#   Rscript scripts/fig_gydb_classes_dark.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

OUT <- "figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Pitched against a pure-white slide headline, so the chrome runs close to
# white: at deck size the earlier mid-greys read as a figure sitting behind the
# title rather than beside it. The point cloud is the one thing deliberately
# held back -- 2,194 dots at full white would out-shout the 26 that matter.
INK <- "#ffffff"        # class names, axis title
INK2 <- "#f7f6fa"       # tick labels, box outlines
FAINT <- "#ffffff"      # the n = labels
CLOUD <- "#cfcddb"      # the non-hit cores
HIT <- "#6cbaff"        # PrLD called -- same blue as the tree slides
VIOLIN <- "#4fbfa4"
VIOLIN_EDGE <- "#7fe8cf"
RULE <- "#b9b7c6"       # the zero line
BORDER <- "#eae7f2"

# Read PLAAC's own output (alpha = 0.5) rather than a reshaped intermediate.
# Protein class is the header prefix GyDB uses to name its cores; PR_ULP1,
# PR_OTU and ORF1_Nterdomain are two-token domain names and keep both tokens,
# otherwise they collapse into a meaningless "PR" / "ORF1".
TWO_TOKEN <- c("PR_ULP1", "PR_OTU", "ORF1_Nterdomain")

protein_class_of <- function(seqid) {
  two <- sub("^([^_]+_[^_]+).*$", "\\1", seqid)
  ifelse(two %in% TWO_TOKEN, two, sub("_.*$", "", seqid))
}

# Same filtering as the paper figure, and for the same reason: a core shorter
# than PLAAC's 60 aa window scores NaN and was never measured, so it is dropped
# rather than counted as a zero. That is what takes CHR out -- 154 cores, 135
# of them too short, leaving 19, under the 20 needed to draw a distribution.
classes <- read_tsv("data/processed/gydb/plaac_data/all_cores_a05.tsv", show_col_types = FALSE) %>%
  transmute(seqid = trimws(SEQid),
            protein_class = protein_class_of(trimws(SEQid)),
            llr = LLR,
            has_prd = as.integer(PRDlen > 0)) %>%
  filter(!is.na(llr)) %>%
  group_by(protein_class) %>% mutate(n_scored = n()) %>% ungroup() %>%
  filter(n_scored >= 20)

order_by_median <- classes %>% group_by(protein_class) %>%
  summarise(m = median(llr), .groups = "drop") %>% arrange(m) %>%
  pull(protein_class)
classes <- mutate(classes,
                  protein_class = factor(protein_class, levels = order_by_median))

# One point per protein. The paper figure jitters ALL cores grey and then
# re-jitters the 26 hits blue on a different seed, so every hit is drawn twice
# at two different x offsets and the cloud holds 2,246 marks for 2,220
# proteins. Splitting the layers is the whole fix.
hits <- filter(classes, has_prd == 1)
rest <- filter(classes, has_prd == 0)
counts <- count(classes, protein_class, name = "n_scored")

p <- ggplot(classes, aes(x = protein_class, y = llr)) +
  geom_hline(yintercept = 0, colour = RULE, linewidth = 0.35) +
  geom_violin(fill = VIOLIN, colour = VIOLIN_EDGE, alpha = 0.26,
              linewidth = 0.4, scale = "width", width = 0.9) +
  geom_boxplot(width = 0.16, outlier.shape = NA, colour = INK2,
               fill = "#ffffff1c", linewidth = 0.35) +
  geom_jitter(data = rest, width = 0.22, height = 0, size = 0.5,
              alpha = 0.38, colour = CLOUD) +
  geom_point(data = hits, size = 2.4, shape = 21, fill = HIT, colour = INK,
             stroke = 0.4,
             position = position_jitter(width = 0.22, height = 0, seed = 1)) +
  geom_text(data = counts, aes(y = max(classes$llr) * 1.13,
                               label = paste0("n = ", n_scored)),
            inherit.aes = TRUE, size = 3.4, colour = FAINT,
            fontface = "italic") +
  scale_y_continuous(breaks = seq(-60, 60, 20),
                     expand = expansion(mult = c(0.03, 0.10))) +
  labs(x = NULL, y = "PLAAC log-likelihood ratio (alpha = 0.5, core 60)") +
  theme_minimal(base_size = 11) +
  theme(
    text             = element_text(colour = INK),
    # The axis title is already pure white, so the only lever left for it is
    # size -- at 11 pt on a 8.8 in figure it read grey purely from thinness.
    axis.title.y     = element_text(colour = INK, size = 12.5),
    axis.text.y      = element_text(colour = INK2, size = 10.5),
    axis.text.x      = element_text(colour = INK, size = 11.5, face = "bold"),
    panel.grid       = element_blank(),
    panel.border     = element_rect(colour = BORDER, fill = NA, linewidth = 0.7),
    # Every pane NA, not a dark fill: the transparency has to go all the way
    # through or the figure lands on the slide inside its own box.
    panel.background = element_rect(fill = NA, colour = NA),
    plot.background  = element_rect(fill = NA, colour = NA),
    legend.position  = "none",
    plot.margin      = margin(6, 8, 4, 4))

# type = "cairo" is not optional: R on Windows defaults the png device to GDI,
# which has no partial alpha and would quantise this to a handful of colours
# with on/off transparency, throwing away every anti-aliased edge.
ggsave(file.path(OUT, "fig_gydb_classes_dark.png"), p, width = 8.8,
       height = 6.0, dpi = 400, bg = "transparent", type = "cairo")

cat(sprintf("wrote fig_gydb_classes_dark   %d cores, %d classes, %d PrLD (%s)\n",
            nrow(classes), length(unique(classes$protein_class)), nrow(hits),
            paste(sort(unique(as.character(hits$protein_class))),
                  collapse = ", ")))
print(as.data.frame(counts))
