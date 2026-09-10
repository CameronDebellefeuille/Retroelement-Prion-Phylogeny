# The REXdb RT cladograms with the full-length Gag trait ringed around them.
#
#   fig_copia_nterm.{png,pdf}   N-terminal disorder   fig_gypsy_nterm.{png,pdf}
#   fig_copia_all.{png,pdf}     whole-Gag disorder    fig_gypsy_all.{png,pdf}
#
# Same construction as the GyDB figure: cladogram so every tip sits at one
# radius, disorder as bar LENGTH rather than colour (the values bunch in a range
# no one reads reliably as lightness), and a called PrLD as the only tip marker,
# in a navy reserved for it so a marker can never be read as a lineage.
#
# Branch colour is lineage, but only for the four largest -- four hues is what
# survives the colourblind check when any two clades can end up side by side on
# a circle. Every other lineage is grey and named on the rim instead, so identity
# is carried by the label rather than by a fifth hue nobody can separate.
#
# 688 copia tips against GyDB's 272, so tip names are dropped; the rim carries
# lineage names with the median trait and n.
#
# Midpoint-rooted for display only -- a circle has to start somewhere, and
# rooting is still an open decision. Nothing here asserts an outgroup.
#
# Run from the repo root, in the conda env:
#   Rscript scripts/plot_tree_full.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(ape); library(ggtree)
})

OUT <- "data/tree_full/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"
HIT <- "#0d366b"; GUIDE <- "#d5d4cf"; OTHER <- "#8a8a85"

# Four hues, validated all-pairs for protan/deutan separation on this surface.
HUES <- c("#D2551B", "#0E8A5F", "#7B2D9E", "#2D7DD2")

FOCUS <- list(copia = c("Ale", "Ivana", "Tork", "SIRE"),
              gypsy = c("Reina", "CRM", "Retand", "Athila"))

LABEL_MIN <- 15      # rim gets a lineage name at this many tips or more
BAR_BASE <- 1.04     # bars start here, as a multiple of the tree radius
BAR_SPAN <- 0.26     # a trait value of 1.00 reaches this far
LAB_GAP <- 0.06
OPEN_ANGLE <- 16      # the wedge left open at 3 o'clock, in degrees


spread <- function(y, gap, lo, hi) {
  # Push labels apart, then pull them back down if the run overshoots -- the
  # band ends short of the first and last tip because the open wedge sits
  # between them and belongs to the bar scale. Two passes: a forward-only pass
  # parks the overflow on top of that scale.
  for (i in seq_along(y)[-1]) y[i] <- max(y[i], y[i - 1] + gap)
  for (i in rev(seq_len(length(y) - 1))) y[i] <- min(y[i], y[i + 1] - gap)
  if (max(y) > hi) y <- y - (max(y) - hi)
  pmax(y, lo)
}


longest_run <- function(ys) {
  # a lineage split into blocks around the circle gets its label on the biggest
  # block, not on the mean, which would land inside someone else's clade
  ys <- sort(ys)
  breaks <- c(0, which(diff(ys) > 1), length(ys))
  best <- which.max(diff(breaks))
  mean(range(ys[(breaks[best] + 1):breaks[best + 1]]))
}


build <- function(short, column, regions, subtitle, file) {
  tips <- read_tsv(sprintf("data/cluster/tips_%s.tsv", short),
                   show_col_types = FALSE) %>%
    select(rexdb_id, lineage = lineage_leaf)
  traits <- read_tsv("data/gag_full/traits_full.tsv", show_col_types = FALSE) %>%
    select(rexdb_id, has_prd, prd_region, disorder_scored, frac_disordered,
           disorder_nterm)

  tree <- phangorn::midpoint(read.tree(sprintf("data/tree_full/%s.treefile", short)))
  data <- tips %>% inner_join(traits, by = "rexdb_id") %>%
    filter(rexdb_id %in% tree$tip.label) %>%
    mutate(trait = .data[[column]],
           focus = ifelse(lineage %in% FOCUS[[short]], lineage, "other"))

  palette <- c(setNames(HUES, FOCUS[[short]]), other = OTHER)
  groups <- split(data$rexdb_id, data$focus)

  p <- ggtree(groupOTU(tree, groups), layout = "circular", open.angle = OPEN_ANGLE,
              branch.length = "none", aes(colour = group), linewidth = 0.28) +
    scale_colour_manual(values = palette, guide = "none")

  radius <- max(p$data$x, na.rm = TRUE)
  base <- radius * BAR_BASE
  span <- radius * BAR_SPAN

  tipdata <- p$data %>% filter(isTip) %>%
    left_join(data, by = c("label" = "rexdb_id")) %>%
    # a tip whose sequence carries an X has no disorder score; it keeps its
    # branch colour and gets no bar, rather than a bar of zero
    mutate(bar_end = base + ifelse(is.na(trait), 0, trait) * span)

  # the open wedge, in tip units: no label may enter it, the bar scale lives there
  wedge <- nrow(tipdata) * OPEN_ANGLE / (360 - OPEN_ANGLE)

  labels <- tipdata %>% group_by(lineage) %>%
    summarise(y = longest_run(y), n = n(),
              med = median(trait, na.rm = TRUE),
              focus = first(focus), .groups = "drop") %>%
    filter(n >= LABEL_MIN) %>%
    arrange(y) %>%
    # two small clades a few tips apart get labels on top of each other, so the
    # angles are spread. The anchor is still the clade; only the angle moves.
    mutate(y = spread(y, nrow(tipdata) * 0.045, 1 + wedge * 0.5,
                      nrow(tipdata) - wedge * 0.5),
           text = sprintf("%s\n%.2f (n = %d)", lineage, med, n))

  # a called domain is marked only if it sits in the region this figure is
  # about: 53 of the 54 copia calls are C-terminal, so an unfiltered marker on
  # the N-terminal figure would point at a domain the ring does not measure
  marks <- filter(tipdata, has_prd == 1, prd_region %in% regions)
  rings <- tibble(x = base + c(0.25, 0.5, 0.75, 1) * span,
                  lab = c("", "0.5", "", "1.0"))
  # the scale sits in the middle of the open wedge, the only arc with no tips
  label_x <- base + span + radius * LAB_GAP

  plot <- p +
    geom_vline(xintercept = rings$x, colour = GUIDE, linewidth = 0.25) +
    geom_vline(xintercept = base, colour = "#a9a8a3", linewidth = 0.3) +
    geom_segment(data = filter(tipdata, !is.na(trait)),
                 aes(x = base, xend = bar_end, y = y, yend = y, colour = focus),
                 inherit.aes = FALSE, linewidth = 0.42, lineend = "butt") +
    geom_point(data = marks, aes(x = radius, y = y, shape = "PrLD called"),
               inherit.aes = FALSE, size = 1.9, stroke = 0.25, colour = INK,
               fill = HIT) +
    geom_text(data = labels, aes(x = label_x, y = y, label = text,
                                 colour = focus),
              inherit.aes = FALSE, size = 3.6, lineheight = 1.0,
              fontface = "bold", show.legend = FALSE) +
    annotate("text", x = rings$x, y = nrow(tipdata) + wedge / 2,
             label = rings$lab, size = 2.6,
             colour = INK2, hjust = 0.5) +
    scale_shape_manual(values = c("PrLD called" = 21), name = NULL,
                       drop = FALSE, limits = "PrLD called") +
    xlim(NA, label_x + radius * 0.16) +
    guides(shape = guide_legend(override.aes = list(size = 2.6, fill = HIT,
                                                    colour = INK))) +
    labs(subtitle = sprintf("%s   |   marked: %d of %d tips with a PrLD called %s",
                            subtitle, nrow(marks), nrow(tipdata),
                            if (length(regions) == 1) "in the N-terminal region"
                            else "anywhere in Gag")) +
    theme_void(base_size = 11) +
    theme(
      text             = element_text(colour = INK),
      plot.subtitle    = element_text(colour = INK2, size = 9.5, hjust = 0.5,
                                      margin = margin(t = 4, b = 2)),
      plot.background  = element_rect(fill = SURFACE, colour = NA),
      legend.position  = "right",
      legend.text      = element_text(size = 9),
      legend.margin    = margin(l = 0, r = 2),
      plot.margin      = margin(2, 2, 2, 2))

  ggsave(file.path(OUT, paste0(file, ".png")), plot, width = 10.5,
         height = 9.4, dpi = 400, bg = SURFACE)
  ggsave(file.path(OUT, paste0(file, ".pdf")), plot, width = 10.5,
         height = 9.4, bg = SURFACE)
  cat(sprintf("wrote %s -- %d tips, %d PrLD, median %s %.3f\n", file,
              nrow(tipdata), nrow(marks), column,
              median(tipdata$trait, na.rm = TRUE)))
}


for (short in c("copia", "gypsy")) {
  build(short, "disorder_nterm", "N-term",
        "N-terminal disorder (mean metapredict score, upstream of the capsid core)",
        sprintf("fig_%s_nterm", short))
  build(short, "frac_disordered", c("N-term", "middle", "C-term"),
        "Whole-Gag disorder (fraction of residues scoring >= 0.5)",
        sprintf("fig_%s_all", short))
}
