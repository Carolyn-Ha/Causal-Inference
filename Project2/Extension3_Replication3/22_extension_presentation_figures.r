# ============================================================
# 22_extension_presentation_figures.r
#
# Goal:
# Create the presentation figures used for the extension intro
# / methodology slide:
#   1. Example preferential-attachment networks
#   2. A three-panel network comparison figure
#   3. A simple table summarizing the simulation design
#
# This script is intentionally standalone and does NOT run any
# simulations. It only creates illustrative visuals.
# ============================================================


# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

library(igraph)
library(tidyverse)
library(gridExtra)
library(grid)

set.seed(6165)

output_dir <- "results/presentation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Settings for illustrative networks
# ------------------------------------------------------------

n_nodes <- 50
m_edges <- 2
attachment_levels_plot <- c(0, 1, 2)

# These are the attachment levels used in the actual extension
attachment_levels_all <- c(0, 0.5, 1, 1.5, 2)


# ------------------------------------------------------------
# 3. Generate example preferential-attachment networks
# ------------------------------------------------------------

example_networks <- lapply(
  attachment_levels_plot,
  function(a) {
    sample_pa(
      n = n_nodes,
      power = a,
      m = m_edges,
      directed = FALSE,
      algorithm = "psumtree"
    )
  }
)

names(example_networks) <- paste0("a_", attachment_levels_plot)


# ------------------------------------------------------------
# 4. Helper function to plot one network
# ------------------------------------------------------------

plot_network_example <- function(g, a_value, file_path = NULL) {

  # A fixed layout makes each panel visually stable
  lay <- layout_with_fr(g)

  if (!is.null(file_path)) {
    png(
      filename = file_path,
      width = 1800,
      height = 1800,
      res = 250
    )
  }

  par(mar = c(0.5, 0.5, 2, 0.5))

  plot(
    g,
    layout = lay,
    vertex.size = 6,
    vertex.label = NA,
    edge.arrow.mode = 0,
    main = paste0("Preferential attachment: a = ", a_value)
  )

  if (!is.null(file_path)) {
    dev.off()
  }
}


# ------------------------------------------------------------
# 5. Save each example network separately
# ------------------------------------------------------------

for (i in seq_along(example_networks)) {

  a_val <- attachment_levels_plot[i]
  g <- example_networks[[i]]

  file_name <- paste0(
    "slide1_network_a",
    gsub("\\.", "_", as.character(a_val)),
    ".png"
  )

  plot_network_example(
    g = g,
    a_value = a_val,
    file_path = file.path(output_dir, file_name)
  )
}


# ------------------------------------------------------------
# 6. Save combined three-panel network figure
# ------------------------------------------------------------

png(
  filename = file.path(output_dir, "slide1_network_panel.png"),
  width = 4200,
  height = 1500,
  res = 250
)

par(mfrow = c(1, 3), mar = c(0.5, 0.5, 2, 0.5))

for (i in seq_along(example_networks)) {
  g <- example_networks[[i]]
  a_val <- attachment_levels_plot[i]

  lay <- layout_with_fr(g)

  plot(
    g,
    layout = lay,
    vertex.size = 6,
    vertex.label = NA,
    edge.arrow.mode = 0,
    main = paste0("a = ", a_val)
  )
}

dev.off()


# ------------------------------------------------------------
# 7. Create a simulation-design summary table
# ------------------------------------------------------------

design_table <- tibble(
  Component = c(
    "Network model",
    "Attachment levels",
    "Networks per attachment level",
    "Nodes per network",
    "Edges per network",
    "Seeding strategies compared",
    "Assignment probability",
    "Seeds per network",
    "Adoption periods",
    "Adoption parameters",
    "Monte Carlo repetitions"
  ),
  Value = c(
    "Preferential attachment",
    "{0, 0.5, 1, 1.5, 2}",
    "150",
    "50",
    "Approximately 97",
    "One-hop vs random",
    "0.5 / 0.5",
    "2",
    "3",
    expression(alpha %in% c(-2, -1) * ", " *
                 beta %in% c(0, 2, 4) * ", " *
                 gamma %in% c(0, 0.1)),
    "2,000 per setting"
  )
)

# Convert expression row to character for table display
design_table$Value <- c(
  "Preferential attachment",
  "{0, 0.5, 1, 1.5, 2}",
  "150",
  "50",
  "Approximately 97",
  "One-hop vs random",
  "0.5 / 0.5",
  "2",
  "3",
  "alpha in {-2,-1}, beta in {0,2,4}, gamma in {0,0.1}",
  "2,000 per setting"
)

table_grob <- tableGrob(
  design_table,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 12,
    core = list(fg_params = list(hjust = 0, x = 0.02)),
    colhead = list(fg_params = list(fontface = "bold"))
  )
)

png(
  filename = file.path(output_dir, "slide1_design_table.png"),
  width = 2600,
  height = 1400,
  res = 220
)

grid.newpage()
grid.draw(table_grob)

dev.off()


# ------------------------------------------------------------
# 8. Optional console output
# ------------------------------------------------------------

cat("\nPresentation figures saved to:\n", output_dir, "\n")
cat("Created files:\n")
cat(" - slide1_network_a0.png\n")
cat(" - slide1_network_a1.png\n")
cat(" - slide1_network_a2.png\n")
cat(" - slide1_network_panel.png\n")
cat(" - slide1_design_table.png\n")
