#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
  library(purrr)
  library(scales)
})
# ---------- CLI ----------
opt <- OptionParser() |>
  add_option(c("-p", "--path"),    type = "character", help = "Path to merged_summary CSV", metavar = "FILE") |>
  add_option(c("-o", "--outfile"), type = "character", help = "Optional output image path (png/pdf/svg)") |>
  add_option(c("-m", "--othermethod"), type = "character", default="baseMemoir", help = "used for naming the other method") |>
  add_option(c("--bins"),          type = "integer",   default = 150, help = "Shared number of bins [default %default]") |>
  add_option(c("--threshold"),     type = "double",    default = 0.70, help = "Vertical threshold line x [default %default]") |>
  parse_args()

if (is.null(opt$path)) {
  stop("Please provide --path to the CSV.")
}
# ---------- Load & prepare data ----------
merged_geno <- read_csv(opt$path, show_col_types = FALSE)

# path = "/Users/gc3045/git/laml2-experiments/real_data/baseMemoir/plots/colony2_merged_summary.csv"
# merged_geno <- read_csv(path)

col_uned_agree <- "#228833"  # green (Agree, unedited)
col_uned_dis   <- "#EE6677"  # coral (Disagree, unedited)
col_ed_agree   <- "#4477AA"  # blue (Agree, edited)
col_ed_dis_u   <- "#CCBB44"  # mustard (Disagree: LP unedited)
col_ed_dis_d   <- "#AA3377"  # wine (Disagree: LP different edit)

# -------------------------------------------------------------------
# Expect merged_geno to have at least: pub_geno (0=unedited, !=0 edited),
# lp_state (0=unedited, !=0 edited, -1 missing), and pub_pmax in [0,1].
# -------------------------------------------------------------------

# Coerce needed cols & drop missings for LP
df <- merged_geno %>%
  mutate(
    pub_geno  = as.integer(pub_geno),
    lp_state = as.integer(lp_state)
  ) %>%
  filter(!is.na(pub_pmax), lp_state != -1L)

# =======================
# A) Left "heatmap" panel
# =======================

uu <- sum(df$pub_geno == 0L & df$lp_state == 0L)
ue <- sum(df$pub_geno == 0L & df$lp_state != 0L)
eu <- sum(df$pub_geno != 0L & df$lp_state == 0L)
ee <- sum(df$pub_geno != 0L & df$lp_state != 0L)

same_edited <- sum(df$pub_geno != 0L & df$lp_state != 0L & df$pub_geno == df$lp_state)
diff_edited <- ee - same_edited

heat_df <- tibble::tribble(
  ~pub,       ~lp,        ~n,            ~label,         ~col,
  "Unedited","Unedited", uu, sprintf("Agree:\n%s", comma(uu)), "#228833",  # orange-ish
  "Unedited","Edited",   ue, sprintf("Disagree:\n%s", comma(ue)), "#EE6677",# light blue
  "Edited",  "Unedited", eu, sprintf("Disagree:\n%s", comma(eu)), "#CCBB44",
  "Edited",  "Edited",   ee, "",                                     "#4477AA"
) %>%
  mutate(
    pub = factor(pub, levels = c("Edited","Unedited")), # so Edited is lower row (like your mock)
    lp = factor(lp, levels = c("Unedited","Edited"))
  )

# assume heat_df, same_edited, diff_edited are already defined as before

# subset for the three non-diagonal cells; keep per-row color in 'col'
label_df <- heat_df %>%
  filter(!(pub == "Edited" & lp == "Edited"))

# numeric coordinates for Edited–Edited cell
ee_x <- which(levels(heat_df$lp) == "Edited")
ee_y <- which(levels(heat_df$pub) == "Edited")

sz_cell_text <- 6     # "Agree/Disagree" in non-EE cells
sz_ee_text   <- 6     # the two labels inside Edited–Edited
sz_axis_text <- 14
sz_axis_lab  <- 15

p_heat <- ggplot(heat_df, aes(x = lp, y = pub)) +
  geom_tile(fill = "white", color = "black", linewidth = 0.8, width = 0.98, height = 0.98) +
  # diagonal split in Edited–Edited tile
  geom_segment(
    data = subset(heat_df, pub == "Edited" & lp == "Edited"),
    aes(x = as.numeric(lp) - 0.49, xend = as.numeric(lp) + 0.49,
        y = as.numeric(pub) - 0.49, yend = as.numeric(pub) + 0.49),
    inherit.aes = FALSE, linewidth = 1, color = "black"
  ) +
  # non-diagonal cell labels (bigger)
  geom_text(
    data = label_df,
    aes(label = label, colour = col),
    size = sz_cell_text, lineheight = 0.95, show.legend = FALSE
  ) +
  scale_colour_identity() +
  # Edited–Edited split labels (bigger)
  annotate("text", x = ee_x - 0.35, y = ee_y + 0.2,
           label = sprintf("Disagree:\n%s", comma(diff_edited)),
           hjust = 0, vjust = 1, color = "#AA3377",
           size = sz_ee_text, lineheight = 0.95) +
  annotate("text", x = ee_x + 0.35, y = ee_y - 0.2,
           label = sprintf("Agree:\n%s", comma(same_edited)),
           hjust = 1, vjust = 0, color = "#4477AA",
           size = sz_ee_text, lineheight = 0.95) +
  labs(x = "LAML-Pro genotype call", y = paste0(opt$othermethod, " genotype call")) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.title  = element_text(color = "black"),
    axis.text.x = element_text(color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 13, angle = 90, vjust = 0.5)
  )

p_heat +
  theme(
    # shrink space between Y *tick labels* and panel
    axis.text.y  = element_text(
      color = "black", size = 13, 
      angle = 90, vjust = 0.5,
      margin = margin(r = 0)   # ↓ smaller = closer; try 0–3
    ),
    # shrink space between Y *axis title* and tick labels
    axis.title.y = element_text(
      color = "black", size = 14,
      margin = margin(r = 0)   # default ~10; reduce a lot
    ),
    # shorter tick marks to save a bit more room
    axis.ticks.length = unit(1.5, "pt"),
    # if you still need to pull the whole panel left/right
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )
p_heat

# ============================================
# B) Right: two mirrored hist panels (top/bot)
# ============================================

# Unedited (pub==0): Agree if lp==0 else Disagree
unedited_df <- df %>%
  filter(pub_geno == 0L) %>%
  mutate(agreement = if_else(lp_state == 0L, "Agree", "Disagree"))

# Edited (pub!=0): three buckets
edited_df <- df %>%
  filter(pub_geno != 0L) %>%
  mutate(cmp = case_when(
    lp_state == 0L ~ "Disagree: LP unedited",
    lp_state != pub_geno ~ "Disagree: LP different edit",
    TRUE ~ "Agree"
  ),
  cmp = factor(cmp, levels = c("Agree","Disagree: LP unedited","Disagree: LP different edit")))

# Shared bins
bins  <- 150L
edges <- seq(0, 1, length.out = bins + 1L)
lefts <- head(edges, -1)
width <- diff(edges)

# Helper to bin and turn into mirrored bars
normalize_by <- "group"

bin_panel <- function(dat, x_col, group_col, levels, colors,
                      normalize_by = c("panel","group","global"),
                      global_n = NULL) {
  normalize_by <- match.arg(normalize_by)
  x <- dat[[x_col]]
  g <- factor(dat[[group_col]], levels = levels)
  
  # bin
  idx <- cut(x, breaks = edges, include.lowest = TRUE, right = TRUE)
  tab <- as.data.frame(table(idx, g), stringsAsFactors = FALSE)
  
  # bin geometry
  tab$bin_left  <- lefts[match(tab$idx, levels(idx))]
  tab$bin_width <- width[match(tab$idx, levels(idx))]
  tab$col       <- colors[match(tab$g, levels)]
  
  # denominators
  if (normalize_by == "panel") {
    denom_map <- setNames(rep(max(length(x), 1), length(levels)), levels)
  } else if (normalize_by == "group") {
    # count rows per group = total sites in that group within this panel
    n_by_g <- tapply(x, g, length)
    denom_map <- setNames(pmax(as.numeric(n_by_g[levels]), 1), levels)
  } else { # "global"
    stopifnot(!is.null(global_n))
    denom_map <- setNames(rep(max(global_n, 1), length(levels)), levels)
  }
  # proportions per row (bin, group), mirrored for all but first level
  tab$denom <- denom_map[tab$g]
  tab$y <- as.numeric(tab$Freq) / tab$denom
  tab$y[tab$g != levels[1]] <- -tab$y[tab$g != levels[1]]
  tab
}

uned_tab <- bin_panel(
  unedited_df, "pub_pmax", "agreement",
  levels = c("Agree","Disagree"),
  colors = c(col_uned_agree, col_uned_dis),
  normalize_by = normalize_by
)

edit_tab <- bin_panel(
  edited_df, "pub_pmax", "cmp",
  levels = c("Agree","Disagree: LP unedited","Disagree: LP different edit"),
  colors = c(col_ed_agree, col_ed_dis_u, col_ed_dis_d),
  normalize_by = normalize_by
)

# Build the two ggplots sharing x
p_top <- ggplot(uned_tab, aes(x = bin_left, y = y, width = bin_width)) +
  # draw negative first (no border)
  geom_col(data = dplyr::filter(uned_tab, g == "Disagree"),
           position = "identity",
           fill = alpha(col_uned_dis, 0.),
           color = col_uned_dis, linewidth = 0.7) +
  # then positive with a border
  geom_col(data = dplyr::filter(uned_tab, g == "Agree"),
           position = "identity",
           fill = alpha(col_uned_agree, 0.20),
           color = col_uned_agree, linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  scale_y_continuous(labels = function(v) sprintf("%g", abs(v))) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(y = ylab, x = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank())
p_bot <- ggplot(edit_tab, aes(x = bin_left, y = y, width = bin_width)) +
  geom_col(data = dplyr::filter(edit_tab, g != "Agree"),
           position = "identity",
           fill = alpha(edit_tab$col[edit_tab$g != "Agree"], 0.20),
           color = edit_tab$col[edit_tab$g != "Agree"], linewidth = 0.6) +
  geom_col(data = dplyr::filter(edit_tab, g == "Agree"),
           position = "identity",
           fill = alpha(col_ed_agree, 0.20),
           color = col_ed_agree, linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  scale_y_continuous(labels = function(v) sprintf("%g", abs(v))) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(y = ylab, x = paste0(opt$othermethod, " posterior probability")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")


ylab <- if (normalize_by == "group") "prop. within group" else if (normalize_by == "panel") "prop. of sites" else "prop. of all sites"
p_top <- p_top + labs(y = ylab)
p_bot <- p_bot + labs(y = ylab)

# Make symmetric y-limits across both panels
ymax <- max(abs(ggplot_build(p_top)$layout$panel_scales_y[[1]]$range$range),
            abs(ggplot_build(p_bot)$layout$panel_scales_y[[1]]$range$range))
p_top <- p_top + coord_cartesian(ylim = c(-ymax, ymax))
p_bot <- p_bot + coord_cartesian(ylim = c(-ymax, ymax))
p_top <- p_top +
  geom_vline(xintercept = 0.70, linetype = "dashed", linewidth = 0.5, colour = "black")
p_bot <- p_bot +
  geom_vline(xintercept = 0.70, linetype = "dashed", linewidth = 0.5, colour = "black")

big_axis_text <- 14
big_axis_title <- 16

p_top <- p_top +
  theme(
    axis.title.y = element_text(size = big_axis_title, face = "bold"),
    axis.title.x = element_text(size = big_axis_title, face = "bold"),
    axis.text.x  = element_text(size = big_axis_text, face = "bold", color = "black"),
    axis.text.y  = element_text(size = big_axis_text, face = "bold", color = "black")
  )
p_bot <- p_bot +
  theme(
    axis.title.y = element_text(size = big_axis_title, face = "bold"),
    axis.title.x = element_text(size = big_axis_title, face = "bold"),
    axis.text.x  = element_text(size = big_axis_text, face = "bold", color = "black"),
    axis.text.y  = element_text(size = big_axis_text, face = "bold", color = "black")
  )

# =================
# Compose the figure
# =================
right_stack <- p_top / p_bot
final_plot  <- p_heat | right_stack + plot_layout(widths = c(1.1, 1.3))

# final_plot

if (!is.null(opt$outfile)) {
  ggsave(opt$outfile, final_plot, width = 12, height = 4.8, dpi = 300, limitsize = FALSE)
  message("✅ Saved: ", opt$outfile)
} else {
  print(final_plot)
}