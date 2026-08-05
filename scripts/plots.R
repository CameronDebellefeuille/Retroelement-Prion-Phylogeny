# Publication figures from data/traits.tsv and the RT tree.
# Run: Rscript scripts/plots.R      Writes PNG + PDF to data/figures/

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(ape); library(ggtree)
})

dir.create("data/figures", showWarnings = FALSE, recursive = TRUE)

# Reference categorical palette, slots 1-3, unchanged and in order.
COPIA <- "#2a78d6"; GYPSY <- "#eb6834"; TY1 <- "#1baf7a"
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"; SURFACE <- "#fcfcfb"

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
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.key.size  = unit(0.8, "lines")
    )
}

traits <- read_tsv("data/traits.tsv", show_col_types = FALSE) |>
  mutate(group = case_when(is_reference == 1 ~ "Ty1 (yeast)",
                           superfamily == "Ty1/copia" ~ "Ty1/copia (plant)",
                           TRUE ~ "Ty3/gypsy (plant)"),
         group = factor(group, levels = c("Ty1/copia (plant)",
                                          "Ty3/gypsy (plant)", "Ty1 (yeast)")))
pal <- c("Ty1/copia (plant)" = COPIA, "Ty3/gypsy (plant)" = GYPSY,
         "Ty1 (yeast)" = TY1)

# n = 5 for the yeast controls, so they are drawn as points; a box would be a
# meaningless flat line next to boxes summarising thousands of elements.
plants <- filter(traits, group != "Ty1 (yeast)")
yeast  <- filter(traits, group == "Ty1 (yeast)")

# --- 1. architecture: how long is the N-terminal extension --------------------
f1 <- ggplot(mapping = aes(group, gag_upstream)) +
  geom_boxplot(data = plants, aes(fill = group), width = 0.5,
               outlier.size = 0.3, outlier.alpha = 0.2,
               linewidth = 0.3, colour = INK2) +
  geom_point(data = yeast, colour = TY1, size = 1.8, alpha = 0.9) +
  geom_hline(yintercept = 180, linetype = "dashed", colour = INK2, linewidth = 0.3) +
  annotate("text", x = 0.55, y = 196, label = "180 aa", hjust = 0,
           size = 2.9, colour = INK2) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = c("copia\n(plant)", "gypsy\n(plant)", "Ty1\n(yeast)")) +
  coord_cartesian(ylim = c(0, 420)) +
  labs(title = "N-terminal extension length",
       subtitle = "Residues upstream of the capsid core. Yeast controls shown as points (n = 5).",
       x = NULL, y = "N-terminal residues") +
  theme_pub()

# --- 2. mechanism: prion-like domains track that length -----------------------
bins <- plants |>
  mutate(bin = cut(gag_upstream, c(0, 75, 120, 180, 250, 1100),
                   labels = c("0-75", "75-120", "120-180", "180-250", "250+"),
                   right = FALSE)) |>
  filter(!is.na(bin)) |>
  group_by(group, bin) |>
  summarise(n = n(), prd = sum(has_prd == 1), .groups = "drop") |>
  mutate(rate = 100 * prd / n) |>
  filter(n >= 20)

f2 <- ggplot(bins, aes(bin, rate, colour = group, group = group)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.3) +
  scale_colour_manual(values = pal) +
  labs(title = "Prion-like domains track length",
       subtitle = "PLAAC domains per length bin (bins with n >= 20)",
       x = "N-terminal residues available", y = "elements with a PRD (%)") +
  theme_pub()

# --- 3. disorder is a property of the N-terminus, not the core ----------------
longer <- function(d) {
  d |> select(group, `N-terminal region` = disorder_nterm,
              `capsid core` = disorder_core) |>
    pivot_longer(-group, names_to = "region", values_to = "disorder") |>
    filter(!is.na(disorder)) |>
    mutate(region = factor(region, levels = c("N-terminal region", "capsid core")))
}
dis_p <- longer(filter(plants, disorder_scored == 1))
dis_y <- longer(filter(yeast,  disorder_scored == 1))

f3 <- ggplot(mapping = aes(region, disorder)) +
  geom_boxplot(data = dis_p, aes(fill = group), width = 0.55,
               outlier.size = 0.25, outlier.alpha = 0.15,
               linewidth = 0.3, colour = INK2, position = position_dodge(0.65)) +
  geom_point(data = dis_y, colour = TY1, size = 1.6, alpha = 0.9,
             position = position_nudge(x = 0.30)) +
  scale_fill_manual(values = pal) +
  labs(title = "Disorder is N-terminal, not core",
       subtitle = "metapredict mean per-residue disorder. Green points: yeast controls.",
       x = NULL, y = "mean disorder") +
  theme_pub()

# --- 4. PLAAC scores, with the controls in frame ------------------------------
f4 <- ggplot(traits, aes(llr, fill = group, colour = group)) +
  geom_density(alpha = 0.25, linewidth = 0.5, adjust = 1.2) +
  geom_vline(xintercept = 0, colour = INK2, linewidth = 0.3) +
  annotate("text", x = 2, y = 0.088, label = "LLR = 0", hjust = 0,
           size = 2.9, colour = INK2) +
  scale_fill_manual(values = pal) + scale_colour_manual(values = pal) +
  coord_cartesian(xlim = c(-60, 45)) +
  labs(title = "Prion-like score vs the yeast control",
       subtitle = "PLAAC log-likelihood ratio, alpha = 1.0, core length 60",
       x = "PLAAC LLR", y = "density") +
  theme_pub()

save2 <- function(name, plot, w, h) {
  ggsave(paste0("data/figures/", name, ".png"), plot, width = w, height = h, dpi = 300)
  ggsave(paste0("data/figures/", name, ".pdf"), plot, width = w, height = h)
}
save2("fig1_architecture", f1, 5.2, 4)
save2("fig2_mechanism",    f2, 5.6, 4)
save2("fig3_disorder",     f3, 5.6, 4)
save2("fig4_plaac",        f4, 5.6, 4)
save2("panel", (f1 | f2) / (f3 | f4) + plot_annotation(tag_levels = "A"), 11, 8)

# --- 5. the tree, with the trait on it ---------------------------------------
tree <- read.tree("data/tree/rt.treefile")
tip <- traits |> filter(rexdb_id %in% tree$tip.label) |>
  transmute(label = rexdb_id, gag_upstream,
            prd = ifelse(has_prd == 1, "PRD", NA_character_))

f5 <- ggtree(tree, layout = "fan", open.angle = 8, linewidth = 0.15,
             colour = INK2) %<+% tip +
  geom_tippoint(aes(colour = gag_upstream), size = 1.1, alpha = 0.85) +
  scale_colour_gradient(low = "#cfe0f5", high = COPIA, name = "N-term aa",
                        na.value = GRID) +
  geom_tippoint(aes(subset = !is.na(prd)), colour = GYPSY, size = 1.7, shape = 18) +
  labs(title = "Copia RT phylogeny, N-terminal extension mapped on tips",
       subtitle = sprintf(paste("%d tips, LG+G4, SH-aLRT. Orange diamonds mark a PLAAC domain;",
                                "they are clustered, not scattered."),
                          length(tree$tip.label))) +
  theme(plot.title = element_text(face = "bold", size = 11, colour = INK),
        plot.subtitle = element_text(colour = INK2, size = 8.8),
        plot.background = element_rect(fill = SURFACE, colour = NA),
        legend.position = c(0.97, 0.5),
        plot.margin = margin(4, 4, 4, 4))

save2("fig5_tree", f5, 7, 6.4)
cat("wrote data/figures/\n")
