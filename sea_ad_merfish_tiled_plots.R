# Load packages
library(anndata)
library(tidyverse)

# Load spatial coordinates and metadata
hfivead <- read_h5ad("./data/raw/SEAAD_MTG_MERFISH_all-nuclei.2023-05-08 (1).h5ad")
spatial_coordinates <- hfivead$obsm
hfivead_metadata <- hfivead$obs

# Create data frame for coordinates and metadata
hfivead_metadata_obsm <- cbind(hfivead_metadata, spatial_coordinates)

# Remove unused cells
hfivead_metadata_obsm_filtered <- hfivead_metadata_obsm %>%
  filter(`Used in analysis` == 1)

# Plot all tissue samples in a tiled plot
ggplot(data = hfivead_metadata_obsm_filtered, aes(x = X_selected_cell_spatial_tiled.1, y = X_selected_cell_spatial_tiled.2, color = Subclass)) +
  geom_point(size = 0.00000001, show.legend = FALSE) +
  theme_minimal() +
  ggtitle("X_selected_cell_spatial_tiled")