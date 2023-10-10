# load libraries ----
library(dplyr)
library(ggplot2)
library(Seurat)
library(SeuratObject)
library(Giotto) # remotes::install_github("RubD/Giotto") 

# load data ----
merged_cortex = readRDS("data/Hypoxia_Cortex_Final.rds")

dim(merged_cortex)
unique(merged_cortex$CellType)
unique(merged_cortex$CellSubType)
unique(merged_cortex$orig.ident)

Idents(merged_cortex) = "CellType"
DimPlot(merged_cortex, reduction = "umap", label = TRUE, raster = FALSE)


ImageDimPlot(merged_cortex, group.by = 'orig.ident', axes = TRUE)

p1 <- ImageFeaturePlot(merged_cortex, features = "Cux2")
p2 <- ImageDimPlot(merged_cortex, molecules = "Cux2", nmols = 1000, alpha = 1, mols.cols = "red")
p1 + p2


get_sample_data = function(data, sample) {
  
  meta = data@meta.data[data@meta.data[["orig.ident"]] == sample, c("orig.ident", "CellType", 'CellSubType')]
  image = c(paste0(sample, "_ALL"), paste0(sample, "_ALL.1"))
  
  coors = data@images[[image]]@boundaries[["centroids"]]@coords
  coordinates_df = as.data.frame(coors)
  
  cell_names = data@images[[image]]@boundaries[["centroids"]]@cells
  rownames(coordinates_df) = cell_names
  
  meta$cell = rownames(meta)
  coordinates_df$cell = rownames(coordinates_df)
  
  merged = merge(meta, coordinates_df, by = "cell")
  
  return(merged)
}

samps = unique(merged_cortex$orig.ident)
tmp = do.call(rbind, lapply(samps, function(x) get_sample_data(merged_cortex, x)))


get_sample_data = function(data, sample, num_images=2) {
  
  meta = data@meta.data[data@meta.data[["orig.ident"]] == sample, c("orig.ident", "CellType", 'CellSubType')]
  meta$cell = rownames(meta)
  
  image_names = sapply(0:(num_images-1), function(i) ifelse(i==0, paste0(sample, "_ALL"), paste0(sample, "_ALL.", i)))
  
  combined_coordinates_df = do.call(rbind, lapply(image_names, function(image) {
    image_data = data@images[[image]]@boundaries[["centroids"]]
    coordinates_df = as.data.frame(image_data@coords, row.names = image_data@cells)
    coordinates_df$cell = rownames(coordinates_df)
    return(coordinates_df)
  }))
  
  return(merge(meta, combined_coordinates_df, by = "cell"))
}



# Subset Seurat object for a specific sample
specific_sample <- "HX1"  # replace with your desired sample name
subset_data <- subset(merged_cortex, subset = orig.ident == specific_sample)

# Plot the spatial data for the subset
ImageDimPlot(subset_data, axes = TRUE)









