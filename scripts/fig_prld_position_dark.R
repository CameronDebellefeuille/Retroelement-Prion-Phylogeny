# fig_prld_position_dark.{png,pdf} -- where the PrLD sits inside Gag, for a
# dark slide.
#
# The slide version of panel A from scripts/fig_prld_position.R: every called
# GyDB domain drawn on its own Gag to scale. Same data, same geometry,
# repainted light-on-dark on a transparent background.
#
# Stripped for the deck: no panel tag, no title, no subtitle, no zinc knuckles
# and so no legend. The slide carries its own headline, and the knuckles are a
# separate claim that this slide is not making.
#
# Superfamily colours are the brightened dark-slide set from
# scripts/fig_gydb_rt_layers.R -- the paper hues (#D2551B, #0E8A5F, #7B2D9E)
# go muddy against near-black.
#
# Run from the repo root:
#   Rscript scripts/fig_prld_position_dark.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(grid)
})

OUT <- "figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

INK <- "#ffffff"        # element names
INK2 <- "#e4e2ee"       # axis text
BACKBONE <- "#6f6c82"   # the ungoverned length of Gag
GRIDLINE <- "#403e56"

FAMILY <- c("Ty3/Gypsy"    = "#FF8A45",
            "Ty1/Copia"    = "#37D6A0",
            "Retroviridae" = "#C58BF5",
            "Bel/Pao"      = "#FF6BA6")

traits <- read_tsv("data/processed/gydb/traits_272.tsv", show_col_types = FALSE)

called <- traits %>%
  filter(has_prd == 1) %>%
  mutate(across(c(gag_len, prd_start, prd_end, prd_mid_pct), as.numeric),
         superfamily = factor(superfamily, levels = names(FAMILY))) %>%
  arrange(superfamily, prd_mid_pct) %>%
  mutate(tip = factor(tip, levels = tip))

p <- ggplot(called) +
  geom_segment(aes(x = 1, xend = gag_len, y = tip, yend = tip),
               colour = BACKBONE, linewidth = 1.9, lineend = "round") +
  geom_segment(aes(x = prd_start, xend = prd_end, y = tip, yend = tip,
                   colour = superfamily),
               linewidth = 5.2, lineend = "butt") +
  scale_colour_manual(values = FAMILY, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03)),
                     breaks = seq(0, 900, 200)) +
  facet_grid(superfamily ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  labs(x = "residue position in Gag", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    text               = element_text(colour = INK),
    axis.text.y        = element_text(size = 10, colour = INK),
    axis.text.x        = element_text(size = 10, colour = INK2),
    axis.title.x       = element_text(size = 11, colour = INK,
                                      margin = margin(t = 8)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = GRIDLINE, linewidth = 0.35),
    panel.spacing      = unit(7, "pt"),
    strip.placement    = "outside",
    strip.text.y.left  = element_text(angle = 0, hjust = 1, face = "bold",
                                      size = 12),
    legend.position    = "none",
    # Every pane NA rather than a dark fill, or the figure lands on the slide
    # inside its own box.
    panel.background   = element_rect(fill = NA, colour = NA),
    plot.background    = element_rect(fill = NA, colour = NA),
    strip.background   = element_rect(fill = NA, colour = NA),
    plot.margin        = margin(10, 12, 4, 8))

# Each facet strip is rendered as its own grob, so a theme colour applies to
# all of them or none -- there is no vectorised strip.text. Build the gtable
# and repaint the strips individually instead: they come out in facet order,
# which is the order the superfamily factor was levelled in.
paint_text <- function(g, col) {
  if (inherits(g, "text")) {
    g$gp <- modifyList(if (is.null(g$gp)) gpar() else g$gp, list(col = col))
    return(g)
  }
  if (!is.null(g$children)) g$children <- lapply(g$children, paint_text, col = col)
  if (!is.null(g$grobs)) g$grobs <- lapply(g$grobs, paint_text, col = col)
  g
}

gt <- ggplotGrob(p)
strips <- grep("^strip-l", gt$layout$name)
strips <- strips[order(gt$layout$t[strips])]
present <- levels(droplevels(called$superfamily))
stopifnot(length(strips) == length(present))
for (i in seq_along(strips)) {
  gt$grobs[[strips[i]]] <- paint_text(gt$grobs[[strips[i]]],
                                      unname(FAMILY[present[i]]))
}

# type = "cairo": R on Windows defaults the png device to GDI, which has no
# partial alpha and would quantise this to a few colours with on/off
# transparency, throwing away every anti-aliased edge.
ggsave(file.path(OUT, "fig_prld_position_dark.png"), gt, width = 9.2,
       height = 6.0, dpi = 400, bg = "transparent", type = "cairo")

cat(sprintf("wrote fig_prld_position_dark -- %d domains; strips %s\n",
            nrow(called),
            paste(sprintf("%s=%s", present, FAMILY[present]), collapse = " ")))
