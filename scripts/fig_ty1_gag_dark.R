# fig_ty1_gag_dark.{png,pdf} -- the Ty1 Gag domain schematic, for a dark slide.
#
# Redrawn from the version in the deck: same 440 aa bar, same boundaries, same
# two colours, but on a transparent background with white type and white rules
# instead of black. The dashed leader lines are not reproduced -- they were
# slide furniture, not part of the diagram.
#
# Boundaries are drawn at the residue numbers the diagram itself labels. In the
# original, six of the seven ticks sit exactly where a linear 1..440 axis puts
# them and the 355 tick is drawn about nine residues right of its own label;
# that is followed here as a drawing slip rather than a second opinion about
# where the domain ends.
#
# Which span each name belongs to was read off the original by measuring where
# the words are centred, not assumed: PrLD, CA and NAC each sit at the exact
# midpoint of 66-136, 159-355 and 355-401. So NAC labels 355-401, and the last
# stretch to 440 is deliberately unnamed -- as it is in the source diagram.
#
# Run from the repo root:
#   Rscript scripts/fig_ty1_gag_dark.R

suppressPackageStartupMessages({
  library(ggplot2)
})

OUT <- "data/gag_plaac/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ORANGE <- "#df912f"     # sampled from the original
GREEN <- "#3a673e"      # sampled from the original
INK <- "#ffffff"        # was black: type and every rule on the bar

GAG_LEN <- 440
PRLD <- c(66, 136)
BREAKS <- c(1, 66, 136, 159, 355, 401, 440)
CUTS <- c(159, 355, 401)          # internal dividers
BOT <- 0                          # bar geometry, in plot units
TOP <- 1

spans <- data.frame(
  label = c("PrLD", "CA", "NAC"),
  start = c(PRLD[1], 159, 355),
  end   = c(PRLD[2], 355, 401))
spans$mid <- (spans$start + spans$end) / 2

p <- ggplot() +
  # the bar, then the PrLD block over it, then the rules on top of both
  annotate("rect", xmin = 1, xmax = GAG_LEN, ymin = BOT, ymax = TOP,
           fill = ORANGE, colour = NA) +
  annotate("rect", xmin = PRLD[1], xmax = PRLD[2], ymin = BOT, ymax = TOP,
           fill = GREEN, colour = NA) +
  annotate("segment", x = CUTS, xend = CUTS, y = BOT, yend = TOP,
           colour = INK, linewidth = 0.6) +
  annotate("segment", x = PRLD, xend = PRLD, y = BOT, yend = TOP,
           colour = INK, linewidth = 0.6) +
  annotate("rect", xmin = 1, xmax = GAG_LEN, ymin = BOT, ymax = TOP,
           fill = NA, colour = INK, linewidth = 0.7) +
  # ticks and residue numbers above
  annotate("segment", x = BREAKS, xend = BREAKS, y = TOP + 0.10,
           yend = TOP + 0.30, colour = INK, linewidth = 0.5) +
  # Ticks stay exact; only the numerals move. 136 and 159 are 23 residues
  # apart and the two labels need about 22 at this size, so that one pair is
  # eased apart rather than shrinking every number on the figure to suit it.
  annotate("text", x = BREAKS + c(0, 0, -7, 7, 0, 0, 0), y = TOP + 0.42,
           label = BREAKS, colour = INK, size = 5.2, vjust = 0) +
  # domain names below, each centred on its own span
  geom_text(data = spans, aes(x = mid, y = BOT - 0.22, label = label),
            colour = INK, size = 6, vjust = 1) +
  # the row label, outside the bar on the left
  annotate("text", x = -10, y = (BOT + TOP) / 2, label = "Ty1 Gag",
           colour = INK, size = 6, hjust = 1) +
  scale_x_continuous(limits = c(-115, 460), expand = c(0, 0)) +
  scale_y_continuous(limits = c(BOT - 0.75, TOP + 0.95), expand = c(0, 0)) +
  theme_void() +
  theme(
    # Every pane NA rather than a dark fill, or the schematic lands on the
    # slide inside its own box.
    plot.background  = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA),
    plot.margin      = margin(4, 6, 4, 6))

# type = "cairo": R on Windows defaults the png device to GDI, which has no
# partial alpha and would quantise this to a few colours with on/off
# transparency, throwing away every anti-aliased edge.
ggsave(file.path(OUT, "fig_ty1_gag_dark.png"), p, width = 9.0, height = 2.1,
       dpi = 400, bg = "transparent", type = "cairo")
ggsave(file.path(OUT, "fig_ty1_gag_dark.pdf"), p, width = 9.0, height = 2.1,
       bg = "transparent", device = cairo_pdf)

cat(sprintf("wrote fig_ty1_gag_dark -- %d aa, PrLD %d-%d, spans %s\n",
            GAG_LEN, PRLD[1], PRLD[2],
            paste(sprintf("%s %d-%d", spans$label, spans$start, spans$end),
                  collapse = "  ")))
