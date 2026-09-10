# One figure: the cladogram with disorder as a bar ring and PrLDs as markers.
#   fig_gydb_clad_combined.{png,pdf}
#
# Disorder is continuous, so it is encoded as bar LENGTH, not colour -- most of
# the values sit between 0.17 and 0.42 and nobody reads a 0.25 lightness
# difference reliably. Bars take the superfamily colour, so the ring needs no
# legend of its own and reads as four blocks. PrLD calls stay as markers, which
# puts the co-occurrence (hits median 0.50 disordered vs 0.29) on one figure.
#
# No title -- it goes on the slide. Legend on the right.
#
# Run from the repo root, in the conda env:
#   Rscript scripts/figures_combined.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(ape); library(ggtree)
})

OUT <- "figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"
HIT <- "#0d366b"; GUIDE <- "#d5d4cf"

FAMILY <- c("Ty3/Gypsy"    = "#D2551B",
            "Ty1/Copia"    = "#0E8A5F",
            "Retroviridae" = "#7B2D9E",
            "Bel/Pao"      = "#C0166B")

# Radial budget, as multiples of the tree radius.
BAR_BASE <- 1.05     # bars start here
BAR_SPAN <- 0.30     # a fraction of 1.00 reaches this far
LAB_GAP  <- 0.05     # gap between the longest bar and the tip names
CHAR_W   <- 0.0140   # width of one tip-label character, in tree-radius units

traits <- read_tsv("data/processed/gydb/traits_272.tsv", show_col_types = FALSE) %>%
  left_join(read_tsv("data/processed/gydb/disorder_272.tsv", show_col_types = FALSE),
            by = "tip")
tree <- phangorn::midpoint(read.tree("data/processed/gydb/tree_272/rt.treefile"))
groups <- split(traits$tip, traits$superfamily)

# Display-only shortening: CopiaSL_monotypic-Chr01_1s215 -> Chr01_1s215. These
# are 29 characters and alone set the radius the lineage labels have to clear;
# the branch colour already says which superfamily they belong to.
short_name <- function(x) sub("^(Copia|Gypsy)SL_monotypic-", "", x)


longest_run <- function(ys) {
  ys <- sort(ys)
  breaks <- c(0, which(diff(ys) > 1), length(ys))
  best <- which.max(diff(breaks))
  mean(range(ys[(breaks[best] + 1):breaks[best + 1]]))
}


p <- ggtree(groupOTU(tree, groups), layout = "circular", open.angle = 16,
            branch.length = "none", aes(colour = group), linewidth = 0.5) +
  scale_colour_manual(values = FAMILY, guide = "none")

radius <- max(p$data$x, na.rm = TRUE)
base <- radius * BAR_BASE
span <- radius * BAR_SPAN

tipdata <- p$data %>% filter(isTip) %>%
  left_join(traits, by = c("label" = "tip")) %>%
  mutate(display = short_name(label),
         bar_end = base + frac_disordered * span)

label_start <- base + span + radius * LAB_GAP

# Each lineage label clears its OWN longest tip name, so short-named clades keep
# their label close in. The text is CENTRED on its anchor, so it also has to
# clear half its own width -- that half-width was what used to land on the tip
# names in the Retroviridae sector.
LAB_SIZE <- 6.4
lab_char <- CHAR_W * (LAB_SIZE / 1.35) * radius   # one character, in x units

labels <- tipdata %>% group_by(superfamily) %>%
  summarise(y = longest_run(y), n = n(), med = median(frac_disordered),
            reach = max(nchar(display)), .groups = "drop") %>%
  mutate(text = sprintf("%s\n%.2f (n = %d)", superfamily, med, n),
         widest = pmax(nchar(superfamily), nchar(sprintf("%.2f (n = %d)", med, n))),
         x = label_start + reach * CHAR_W * radius +
             widest * lab_char / 2 + radius * 0.03)

marks <- filter(tipdata, has_prd == 1)
rings <- tibble(x = base + c(0.25, 0.5, 0.75, 1) * span,
                lab = c("", "0.5", "", "1.0"))

plot <- p +
  geom_vline(xintercept = rings$x, colour = GUIDE, linewidth = 0.25) +
  geom_vline(xintercept = base, colour = "#a9a8a3", linewidth = 0.3) +
  geom_segment(data = tipdata,
               aes(x = base, xend = bar_end, y = y, yend = y,
                   colour = superfamily),
               inherit.aes = FALSE, linewidth = 0.85, lineend = "butt") +
  # on the tip itself, at the foot of its own bar
  geom_point(data = marks, aes(x = radius, y = y, shape = "PrLD called"),
             inherit.aes = FALSE, size = 2.3, stroke = 0.3, colour = INK,
             fill = HIT) +
  # names: hits bold and in the marker colour, everything else grey
  geom_tiplab(data = function(d)
                mutate(d, label = ifelse(label %in% marks$label, NA,
                                         short_name(label))),
              size = 1.35, offset = label_start - radius, colour = "#54534f",
              na.rm = TRUE) +
  geom_tiplab(data = function(d)
                mutate(d, label = ifelse(label %in% marks$label,
                                         short_name(label), NA)),
              size = 1.6, offset = label_start - radius, colour = HIT,
              fontface = "bold", na.rm = TRUE) +
  geom_text(data = labels, aes(x = x, y = y, label = text,
                               colour = superfamily),
            inherit.aes = FALSE, size = LAB_SIZE, lineheight = 1.0,
            fontface = "bold", show.legend = FALSE) +
  # bar scale, printed once in the open wedge
  annotate("text", x = rings$x, y = 0.6, label = rings$lab, size = 2.6,
           colour = INK2, hjust = 0.5) +
  scale_shape_manual(values = c("PrLD called" = 21), name = NULL) +
  xlim(NA, max(labels$x) + radius * 0.05) +
  guides(shape = guide_legend(override.aes = list(size = 2.8, fill = HIT,
                                                  colour = INK))) +
  theme_void(base_size = 11) +
  theme(
    text             = element_text(colour = INK),
    plot.background  = element_rect(fill = SURFACE, colour = NA),
    legend.position  = "right",
    legend.direction = "vertical",
    legend.text      = element_text(size = 9),
    legend.margin    = margin(l = 0, r = 2),
    plot.margin      = margin(2, 2, 2, 2))

ggsave(file.path(OUT, "fig_gydb_clad_combined.png"), plot, width = 11.5,
       height = 10.2, dpi = 400, bg = SURFACE)
cat("wrote fig_gydb_clad_combined\n")
cat(sprintf("lineage label radii: %s\n",
            paste(sprintf("%s %.2fR", labels$superfamily, labels$x / radius),
                  collapse = "  ")))
