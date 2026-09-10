# Where the Brassicaceae elements sit on the REXdb RT trees.
#
#   fig_brass_copia.{png,pdf}   688 tips, 112 Brassicaceae
#   fig_brass_gypsy.{png,pdf}   433 tips,  46 Brassicaceae
#
# Stripped back from plot_tree_full.R: no lineage colour, no trait ring, no PrLD
# marker. Every branch is grey and the only ink on the tree is a tip point for a
# Brassicaceae element, coloured by species. The question here is where those
# elements are, not what they score.
#
# Six hues, validated all-pairs rather than adjacent-only -- tips of one species
# are scattered round the circle, so any two species can end up side by side and
# the adjacent-pair test would not cover it. Worst pair: CVD (protan/deutan)
# dE 8.6, tritan 8.2, normal vision 15.0, all six clearing 3:1 on the surface.
# The four-hue limit noted in plot_tree_full.R applies to that script's branch
# colouring, where hue must survive being drawn as a hairline; a filled point
# with a dark ring carries a hue at smaller separation than a 0.28pt branch does.
#
# Colours are assigned on the sorted species name, so a species is the same hue
# in both panels.
#
# Drawn as IQ-TREE left it, NOT midpoint-rooted. plot_tree_full.R midpoint-roots
# for display, but phangorn is no longer in the conda env (nor is ggrepel), so
# that script does not currently run. Rooting is open anyway (L-3) and on a
# cladogram it only decides where the circle starts, so nothing here needs it and
# nothing here asserts an outgroup.
#
# Run from the repo root, in the conda env:
#   Rscript scripts/fig_brassicaceae.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(ape); library(ggtree)
})

OUT <- "data/tree_full/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"
BRANCH <- "#c9c8c3"

# validated all-pairs on SURFACE -- see the header
HUES <- c("#8A4B1A", "#1B4FD8", "#C0397F", "#3E93C9", "#4E8A16", "#7A3A9E")

GENERA <- c("Arabidopsis", "Brassica", "Capsella", "Eutrema", "Boechera",
            "Camelina", "Raphanus", "Eutrema", "Thlaspi", "Sisymbrium")

OPEN_ANGLE <- 16


build <- function(short) {
  tips <- read_tsv(sprintf("data/cluster/tips_%s.tsv", short),
                   show_col_types = FALSE) %>%
    select(rexdb_id, species)

  tree <- read.tree(sprintf("data/tree_full/%s.treefile", short))

  data <- tips %>%
    filter(rexdb_id %in% tree$tip.label) %>%
    mutate(brass = grepl(paste(GENERA, collapse = "|"), species))

  species <- sort(unique(data$species[data$brass]))
  palette <- setNames(HUES[seq_along(species)], species)

  p <- ggtree(tree, layout = "circular", open.angle = OPEN_ANGLE,
              branch.length = "none", colour = BRANCH, linewidth = 0.26)

  radius <- max(p$data$x, na.rm = TRUE)

  marks <- p$data %>% filter(isTip) %>%
    left_join(data, by = c("label" = "rexdb_id")) %>%
    filter(brass)

  counts <- marks %>% count(species)
  legend_labels <- setNames(
    sprintf("%s  (n = %d)", counts$species, counts$n), counts$species)

  plot <- p +
    geom_point(data = marks, aes(x = radius, y = y, fill = species),
               inherit.aes = FALSE, shape = 21, size = 2.4, stroke = 0.3,
               colour = INK) +
    scale_fill_manual(values = palette, labels = legend_labels, name = NULL) +
    xlim(NA, radius * 1.06) +
    guides(fill = guide_legend(override.aes = list(size = 3.2))) +
    labs(subtitle = sprintf(
      "REXdb %s RT, %d cluster representatives   |   %d Brassicaceae tips in %d species",
      short, nrow(data), nrow(marks), length(species))) +
    theme_void(base_size = 11) +
    theme(
      text            = element_text(colour = INK),
      plot.subtitle   = element_text(colour = INK2, size = 9.5, hjust = 0.5,
                                     margin = margin(t = 4, b = 2)),
      plot.background = element_rect(fill = SURFACE, colour = NA),
      legend.position = "right",
      legend.text     = element_text(size = 9, face = "italic"),
      plot.margin     = margin(2, 2, 2, 2))

  file <- sprintf("fig_brass_%s", short)
  ggsave(file.path(OUT, paste0(file, ".png")), plot, width = 10, height = 8.2,
         dpi = 400, bg = SURFACE)
  ggsave(file.path(OUT, paste0(file, ".pdf")), plot, width = 10, height = 8.2,
         bg = SURFACE)
  cat(sprintf("wrote %s -- %d of %d tips are Brassicaceae, %d species\n",
              file, nrow(marks), nrow(data), length(species)))
}


for (short in c("copia", "gypsy")) build(short)
