# Where the prion-like domain sits inside Gag.
#
#   fig_prld_position.{png,pdf}
#
# Panel A is the whole GyDB result at element resolution: every one of the 25
# called domains, drawn on its own Gag to scale, with the Pfam zinc knuckles
# beside them. 25 elements is few enough to name, which is the reason to lead
# with GyDB rather than REXdb's 1,131 calls.
#
# Panel B is the same question asked of both databases at once, as the share of
# domains falling in each third of Gag. It carries the cross-database
# comparison: GyDB copia is N-terminal (2 of 3 calls, both yeast Pseudovirus),
# REXdb plant copia is overwhelmingly C-terminal.
#
# Position is a fraction of each element's own Gag, so panel B is comparable
# across proteins of different length; panel A keeps absolute residues so the
# length differences stay visible.
#
# Colour: superfamily in panel A, from the palette used by the other GyDB
# figures. Panel B encodes position along the protein, which is ordered, so it
# takes one hue light-to-dark rather than three unrelated ones.
#
# Run from the repo root, in the conda env:
#   Rscript scripts/fig_prld_position.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork)
})

OUT <- "data/gag_plaac/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"
BACKBONE <- "#d9d8d3"; HIT <- "#0d366b"

FAMILY <- c("Ty3/Gypsy"    = "#D2551B",
            "Ty1/Copia"    = "#0E8A5F",
            "Retroviridae" = "#7B2D9E",
            "Bel/Pao"      = "#C0166B")

# one hue, light to dark, because a third of a protein is an ordered quantity
ZONE <- c("N-terminal" = "#c6d9ee", "middle" = "#5d8ec4", "C-terminal" = "#123f6b")

traits <- read_tsv("data/gag_plaac/traits_272.tsv", show_col_types = FALSE)
knuckle <- read_tsv("data/gag_plaac/gydb_nc_hits.tsv", show_col_types = FALSE) %>%
  filter(nc_strict == 1) %>%
  transmute(tip = rexdb_id, k_start = best_start, k_end = best_end)

called <- traits %>%
  filter(has_prd == 1) %>%
  mutate(across(c(gag_len, prd_start, prd_end, prd_mid_pct), as.numeric),
         superfamily = factor(superfamily, levels = names(FAMILY))) %>%
  arrange(superfamily, prd_mid_pct) %>%
  mutate(tip = factor(tip, levels = tip))

marks <- knuckle %>%
  inner_join(select(called, tip, superfamily), by = "tip") %>%
  mutate(mid = (k_start + k_end) / 2)

panel_a <- ggplot(called) +
  geom_segment(aes(x = 1, xend = gag_len, y = tip, yend = tip),
               colour = BACKBONE, linewidth = 1.6, lineend = "round") +
  geom_segment(aes(x = prd_start, xend = prd_end, y = tip, yend = tip,
                   colour = superfamily),
               linewidth = 4.4, lineend = "butt") +
  geom_point(data = marks, aes(x = mid, y = tip, shape = "zinc knuckle"),
             fill = HIT, colour = INK, size = 2.1, stroke = 0.3) +
  scale_colour_manual(values = FAMILY, guide = "none") +
  scale_shape_manual(values = c("zinc knuckle" = 23), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03)),
                     breaks = seq(0, 900, 200)) +
  facet_grid(superfamily ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "residue position in Gag", y = NULL,
       title = "A   Prion-like domains in the GyDB Gag cores",
       subtitle = "every called domain, on its own Gag to scale") +
  guides(shape = guide_legend(override.aes = list(size = 2.6))) +
  theme_minimal(base_size = 10) +
  theme(
    text             = element_text(colour = INK),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(colour = INK2, size = 9, margin = margin(b = 8)),
    axis.text.y      = element_text(size = 7.6, colour = INK),
    axis.text.x      = element_text(size = 8, colour = INK2),
    axis.title.x     = element_text(size = 8.6, colour = INK2,
                                    margin = margin(t = 6)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "#eceae5", linewidth = 0.3),
    panel.spacing    = unit(6, "pt"),
    strip.placement  = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold",
                                     size = 8.4, colour = INK2),
    legend.position  = "bottom",
    legend.text      = element_text(size = 8.4),
    legend.key.height = unit(9, "pt"),
    plot.background  = element_rect(fill = SURFACE, colour = NA),
    plot.margin      = margin(8, 10, 4, 8))

# --- panel B: both databases, share of domains per third of Gag ---------------

zones_of <- function(df, label) {
  df %>%
    filter(has_prd == 1, !is.na(prd_start)) %>%
    mutate(frac = ((as.numeric(prd_start) + as.numeric(prd_end)) / 2) /
             as.numeric(gag_len),
           zone = case_when(frac < 1/3 ~ "N-terminal",
                            frac <= 2/3 ~ "middle",
                            TRUE ~ "C-terminal")) %>%
    count(zone) %>%
    mutate(set = label, share = n / sum(n), total = sum(n))
}

rexdb <- read_tsv("data/gag_full/gag_full_plaac.tsv", show_col_types = FALSE)

zones <- bind_rows(
  zones_of(filter(traits, superfamily == "Ty1/Copia"), "GyDB  Ty1/Copia"),
  zones_of(filter(traits, superfamily == "Ty3/Gypsy"), "GyDB  Ty3/Gypsy"),
  zones_of(filter(traits, superfamily == "Retroviridae"), "GyDB  Retroviridae"),
  zones_of(filter(rexdb, superfamily == "Ty1/copia"), "REXdb  Ty1/copia"),
  zones_of(filter(rexdb, superfamily == "Ty3/gypsy"), "REXdb  Ty3/gypsy")) %>%
  mutate(zone = factor(zone, levels = names(ZONE)),
         set = factor(set, levels = rev(c(
           "GyDB  Ty1/Copia", "GyDB  Ty3/Gypsy", "GyDB  Retroviridae",
           "REXdb  Ty1/copia", "REXdb  Ty3/gypsy"))))

labels <- distinct(zones, set, total)

panel_b <- ggplot(zones, aes(share, set, fill = zone)) +
  geom_col(width = 0.62, colour = SURFACE, linewidth = 0.6,
           position = position_stack(reverse = TRUE)) +
  geom_text(data = labels, aes(x = 1.02, y = set, label = sprintf("n = %d", total)),
            inherit.aes = FALSE, hjust = 0, size = 2.9, colour = INK2) +
  scale_fill_manual(values = ZONE, name = NULL) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.13),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "share of called domains", y = NULL,
       title = "B   Which third of Gag the domain falls in",
       subtitle = "position as a fraction of each element's own Gag") +
  theme_minimal(base_size = 10) +
  theme(
    text             = element_text(colour = INK),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(colour = INK2, size = 9, margin = margin(b = 8)),
    axis.text.y      = element_text(size = 8.4, colour = INK),
    axis.text.x      = element_text(size = 8, colour = INK2),
    axis.title.x     = element_text(size = 8.6, colour = INK2,
                                    margin = margin(t = 6)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "#eceae5", linewidth = 0.3),
    legend.position  = "bottom",
    legend.text      = element_text(size = 8.4),
    legend.key.size  = unit(10, "pt"),
    plot.background  = element_rect(fill = SURFACE, colour = NA),
    plot.margin      = margin(8, 10, 4, 8))

figure <- panel_a / panel_b + plot_layout(heights = c(2.5, 1))

ggsave(file.path(OUT, "fig_prld_position.png"), figure, width = 7.6,
       height = 9.2, dpi = 400, bg = SURFACE)
ggsave(file.path(OUT, "fig_prld_position.pdf"), figure, width = 7.6,
       height = 9.2, bg = SURFACE)
cat(sprintf("wrote fig_prld_position -- %d GyDB domains, %d knuckles marked\n",
            nrow(called), nrow(marks)))
