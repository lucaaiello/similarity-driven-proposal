# Shared ggplot theme and method palette for manuscript figures.
theme_manuscript <- function(base_size = 11, diagonal = FALSE) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(plot.title = element_blank(), plot.subtitle = element_blank(),
          plot.tag = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "#E4E4E4", linewidth = 0.3),
          axis.text = element_text(colour = "#333333"),
          axis.text.x = element_text(angle = if (diagonal) 40 else 0,
                                    hjust = if (diagonal) 1 else 0.5),
          axis.title = element_text(size = base_size),
          legend.position = "bottom", legend.box = "vertical", legend.title = element_blank(),
          legend.text = element_text(size = base_size - 2),
          legend.key.height = grid::unit(0.35, "cm"),
          plot.margin = margin(5, 8, 5, 5))
}
palette <- c(uniform = "#666666", f_fixed = "#56B4E9", f_online = "#0072B2",
             bic_fixed = "#E69F00", bic_online = "#D55E00", locally_balanced = "#009E73",
             bic_locally_balanced = "#CC79A7", f_flip_only = "#56B4E9",
             f_uniform_swap = "#0072B2", f_informed_swap = "#003F6B",
             bic_flip_only = "#E69F00", bic_uniform_swap = "#D55E00",
             bic_informed_swap = "#8C2D04", uniform_flip = "#666666",
             optimized_similarity_fixed = "#E69F00", online_similarity_flip = "#0072B2",
             online_similarity_local = "#009E73", gagnon_laplace_rj = "#CC79A7",
             f = "#0072B2", bic = "#D55E00", marginal = "#009E73")
