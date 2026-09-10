# The GyDB RT tree with Gag PrLDs mapped onto the tips.
#
#   fig_gydb_rt_full_dark.png   every called domain, PrLD rate per superfamily
#
# 272 tips. The tree is built from RT; the PrLD is a Gag trait joined on element
# name. Tips are named on the rim, PrLD-bearing ones marked and bolded, and each
# superfamily block carries its hit rate as hits/n (pct).
#
# Written light-on-dark on a transparent background, so it drops onto a dark
# slide with no white box around it. build() still takes dark = FALSE for a
# paper-surface version; only the dark one is built.
#
# No disorder ring here -- that is figures_combined.R, which draws the same tree
# with disorder as a bar ring. Both are midpoint-rooted so tip order is
# identical between them. Rooting is display-only (L-3).
#
# This script previously also built a bare tree, an N-terminal-only tree and a
# text-free registered overlay set, as a slide reveal. Those are retired; the
# untrimmed version is in the backup as fig_gydb_rt_layers.R_untrimmed. build()
# still accepts layers=, zones=, frame=, fill= and legend=, so the reveal can be
# rebuilt by adding calls back.
#
# Run from the repo root:
#   Rscript scripts/fig_gydb_rt_layers.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(ape); library(ggtree)
})

OUT <- "figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

INK <- "#0b0b0b"; INK2 <- "#52514e"; SURFACE <- "#fcfcfb"
HIT <- "#0d366b"
GREY <- "#8d8c87"
TIPCOL <- "#54534f"

FAMILY <- c("Ty3/Gypsy"    = "#D2551B",
            "Ty1/Copia"    = "#0E8A5F",
            "Retroviridae" = "#7B2D9E",
            "Bel/Pao"      = "#C0166B")

# Dark-slide palette. Not the paper palette filtered: each ink is re-picked so
# it holds its hue against near-black. The branches ARE the figure on the bare
# slide, so they go near-white and heavier -- a 0.5pt grey hairline that reads
# fine on paper vanishes when projected onto a dark ground.
DARK_INK <- "#f4f3f8"; DARK_GREY <- "#dedbe9"; DARK_TIPCOL <- "#a8a5bb"
DARK_HIT <- "#6cbaff"

DARK_FAMILY <- c("Ty3/Gypsy"    = "#FF8A45",
                 "Ty1/Copia"    = "#37D6A0",
                 "Retroviridae" = "#C58BF5",
                 "Bel/Pao"      = "#FF6BA6")

ALL_ZONES <- c("N-terminal", "middle", "C-terminal")

# Radial budget. No disorder ring, so the names come in close: the annulus that
# held the bars in figures_combined.R is gone, not merely left empty.
TIP_GAP    <- 0.03     # tips -> marker ring
NAME_GAP   <- 0.06     # marker ring -> tip names
CHAR_W     <- 0.0140   # width of one tip-label character, in tree-radius units
LAB_SIZE   <- 6.4
OPEN_ANGLE <- 16
BRANCH_WIDTH <- 0.5
DARK_BRANCH_WIDTH <- 0.85

short_name <- function(x) sub("^(Copia|Gypsy)SL_monotypic-", "", x)

# Push a sorted set of label positions apart until no two are closer than
# `sep`, staying inside [lo, hi]. Both members of an overlapping pair move, so
# a cluster opens around its own centre rather than drifting off one end, and
# the original order is never crossed -- a label can only be read against the
# wrong tip if it overtakes a neighbour, and this cannot.
spread_labels <- function(y, sep, lo, hi) {
  for (i in seq_len(500)) {
    d <- diff(y)
    bad <- which(d < sep)
    if (!length(bad)) break
    push <- (sep - d[bad]) / 2
    y[bad] <- y[bad] - push
    y[bad + 1L] <- y[bad + 1L] + push
    y <- pmin(pmax(y, lo), hi)
  }
  y
}

longest_run <- function(ys) {
  ys <- sort(ys)
  breaks <- c(0, which(diff(ys) > 1), length(ys))
  best <- which.max(diff(breaks))
  mean(range(ys[(breaks[best] + 1):breaks[best + 1]]))
}

traits <- read_tsv("data/processed/gydb/traits_272.tsv", show_col_types = FALSE)
tree <- phangorn::midpoint(read.tree("data/processed/gydb/tree_272/rt.treefile"))
groups <- split(traits$tip, traits$superfamily)


build <- function(file, layers = character(0), zones = ALL_ZONES,
                  fill = FALSE, dark = FALSE, frame = NULL, legend = TRUE,
                  hit_size = 1.6, mark_size = 2.3, hit_sep = 5) {

  # One palette swap up front, so the rest of the function never asks again.
  ink    <- if (dark) DARK_INK else INK
  grey   <- if (dark) DARK_GREY else GREY
  tipcol <- if (dark) DARK_TIPCOL else TIPCOL
  hit    <- if (dark) DARK_HIT else HIT
  family <- if (dark) DARK_FAMILY else FAMILY
  stroke <- if (dark) DARK_BRANCH_WIDTH else BRANCH_WIDTH
  paper  <- if (dark) NA else SURFACE

  # Hit names are the only text on the overlay slides, so they carry the
  # reading and the marker carries the colour coding: near-white beats the
  # marker blue at projector distance, and the dot is already blue.
  hitlab <- if (dark) ink else HIT

  shown <- traits %>% filter(has_prd == 1, prd_zone %in% zones)

  p <- ggtree(groupOTU(tree, groups), layout = "circular",
              open.angle = OPEN_ANGLE, branch.length = "none",
              aes(colour = group), linewidth = stroke) +
    scale_colour_manual(
      values = if ("colour" %in% layers) family
               else setNames(rep(grey, length(family)), names(family)),
      guide = "none")

  radius <- max(p$data$x, na.rm = TRUE)
  mark_x <- radius * (1 + TIP_GAP)
  label_start <- mark_x + radius * NAME_GAP

  tipdata <- p$data %>% filter(isTip) %>%
    left_join(traits, by = c("label" = "tip")) %>%
    mutate(display = short_name(label))

  lab_char <- CHAR_W * (LAB_SIZE / 1.35) * radius

  # Rate is computed on the zone-filtered hits, so each slide reports its own.
  labels <- tipdata %>% group_by(superfamily) %>%
    summarise(y = longest_run(y), n = n(),
              reach = max(nchar(display)), .groups = "drop") %>%
    left_join(count(shown, superfamily, name = "hits"), by = "superfamily") %>%
    mutate(hits = ifelse(is.na(hits), 0L, hits),
           text = sprintf("%s\n%d/%d (%.0f%%)", superfamily, hits, n,
                          100 * hits / n),
           widest = pmax(nchar(superfamily),
                         nchar(sprintf("%d/%d (%.0f%%)", hits, n,
                                       100 * hits / n))),
           x = label_start + reach * CHAR_W * radius +
               widest * lab_char / 2 + radius * 0.03)

  marks <- filter(tipdata, label %in% shown$tip)

  # Where the frame stops, decided before anything is painted so the layers
  # cannot drift apart. `frame` pins it to a fixed multiple of the tree radius;
  # the default derives it from the superfamily blocks, which is fine for a
  # standalone figure but NOT safe to register against -- see FRAME below.
  outer <- if (fill) radius * 1.02
           else if (!is.null(frame)) radius * frame
           else max(labels$x) + radius * 0.05

  plot <- p

  if ("prld" %in% layers) {
    plot <- plot +
      geom_point(data = marks, aes(x = mark_x, y = y, shape = "PrLD called"),
                 inherit.aes = FALSE, size = mark_size, stroke = 0.3,
                 colour = ink, fill = hit) +
      scale_shape_manual(values = c("PrLD called" = 21), name = NULL) +
      guides(shape = guide_legend(override.aes = list(size = 2.8, fill = hit,
                                                      colour = ink)))
  }

  if ("tips" %in% layers) {
    plot <- plot +
      geom_tiplab(data = function(d)
                    mutate(d, label = ifelse(label %in% shown$tip, NA,
                                             short_name(label))),
                  size = 1.35, offset = label_start - radius,
                  colour = tipcol, na.rm = TRUE) +
      geom_tiplab(data = function(d)
                    mutate(d, label = ifelse(label %in% shown$tip,
                                             short_name(label), NA)),
                  size = 1.6, offset = label_start - radius, colour = hit,
                  fontface = "bold", na.rm = TRUE)
  }

  # "hits" is "tips" with the other 247 names deleted rather than greyed out.
  # Nothing else competes for the annulus, so the names that are left can run
  # several times the size they do on the labelled figures.
  #
  # 25 names on a 272-tip circle is not a density problem on average -- it is a
  # clustering problem. Nomad/HMS-Beagle/Yoyo/Burdock are four CONSECUTIVE
  # tips, and at a size worth reading their names sit on top of each other.
  # Radial lanes are the obvious fix and the wrong one here: four lanes of ten
  # characters would eat more radius than the tree. So the labels slide along
  # the ring instead, and a connector puts each one back on its tip.
  if ("hits" %in% layers) {
    tp <- filter(p$data, isTip)
    slope <- diff(range(tp$angle)) / diff(range(tp$y))
    base <- min(tp$angle) - slope * min(tp$y)

    lab <- marks %>% arrange(y) %>%
      mutate(y_lab = spread_labels(y, hit_sep, min(tp$y), max(tp$y)),
             ang = base + slope * y_lab,
             # Left half reads upside down unless it is turned and anchored
             # from the other end, exactly as geom_tiplab does it.
             flip = ang > 90 & ang < 270,
             tang = ifelse(flip, ang + 180, ang),
             hj = ifelse(flip, 1, 0))

    link <- bind_rows(
      transmute(lab, id = label, x = mark_x, y = y),
      transmute(lab, id = label, x = label_start - radius * 0.02, y = y_lab))

    plot <- plot +
      geom_path(data = link, aes(x = x, y = y, group = id),
                inherit.aes = FALSE, colour = hit, linewidth = 0.35,
                alpha = 0.8) +
      geom_text(data = lab, aes(x = label_start, y = y_lab, label = display,
                                angle = tang, hjust = hj),
                inherit.aes = FALSE, size = hit_size, colour = hitlab,
                fontface = "bold")
  }

  if ("blocks" %in% layers) {
    plot <- plot +
      geom_text(data = labels, aes(x = x, y = y, label = text,
                                   colour = superfamily),
                inherit.aes = FALSE, size = LAB_SIZE, lineheight = 1.0,
                fontface = "bold", show.legend = FALSE)
  }

  # On the dark build every one of these panes is NA, not a dark fill: the hole
  # has to go all the way through, or the legend key punches a box in the slide.
  plot <- plot +
    xlim(NA, outer) +
    theme_void(base_size = 11) +
    theme(
      text              = element_text(colour = ink),
      plot.background   = element_rect(fill = paper, colour = NA),
      panel.background  = element_rect(fill = paper, colour = NA),
      legend.background = element_rect(fill = paper, colour = NA),
      legend.key        = element_rect(fill = paper, colour = NA),
      legend.position   = if (legend && "prld" %in% layers) "right" else "none",
      legend.direction  = "vertical",
      legend.text       = element_text(size = 9, colour = ink),
      legend.margin     = margin(l = 0, r = 2),
      plot.margin       = margin(2, 2, 2, 2))

  # type = "cairo" is not optional. R on Windows defaults the png device to
  # GDI, which cannot do partial alpha: it quantises to a 6-colour palette with
  # on/off transparency, so every branch comes out aliased and the soft edges
  # that make a hairline readable are thrown away. Cairo writes real RGBA.
  bg <- if (dark) "transparent" else SURFACE
  ggsave(file.path(OUT, paste0(file, ".png")), plot, width = 11.5,
         height = 10.2, dpi = 400, bg = bg, type = "cairo")
  cat(sprintf("wrote %-27s %3d marks   %s\n", file, nrow(marks),
              paste(sprintf("%s %d/%d", labels$superfamily, labels$hits,
                            labels$n), collapse = "  ")))
}


# One tree: every called domain, tips named, superfamily blocks and hit counts
# on the rim. Light and dark are the same figure repainted.
#
# The bare and N-terminal layers, and the text-free registered overlay set that
# went with them, were a presentation reveal and are no longer built here. The
# untrimmed script is in the backup as fig_gydb_rt_layers.R_untrimmed if that
# sequence is ever wanted again; build() still takes the layers=, zones=,
# frame=, fill= and legend= arguments they used.
build("fig_gydb_rt_full_dark", dark = TRUE,
      layers = c("colour", "prld", "tips", "blocks"), zones = ALL_ZONES)
