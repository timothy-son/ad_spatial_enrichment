
# load libraries ----
library(dplyr)
library(ggplot2)
library(data.table)
library(Seurat)
library(SeuratObject)

# devtools::install_github("drieslab/Giotto@suite")
library(Giotto) 
# devtools::install_github("RubD/GiottoUtils")
# devtools::install_github("RubD/GiottoClass")
# devtools::install_github("drieslab/GiottoData")
library(GiottoData) 
if(!checkGiottoEnvironment()){
  installGiottoEnvironment() # run once 
}

# spatial enrichment analysis ----

# for each sample 
# - get cell identifiers 
# - get x/y coordinates from image slot 
get_spat_data = function(data, samples) {
  do.call(rbind, lapply(samples, function(samp) {
    img = data@images[[samp]]
    data.frame(x = img@boundaries[["centroids"]]@coords[, 1], 
               y = img@boundaries[["centroids"]]@coords[, 2], 
               cell_ID = img@boundaries[["centroids"]]@cells)
  }))
}

# for each sample 
# - create Giotto object 
# - create Delaunay network 
# - calculate cell proximity enrichments 
get_cp_enrichments = function(samp, sobject, spatial_data) {
  meta = subset(sobject@meta.data, subset = orig.ident == samp)
  meta$cell_ID = rownames(meta)
  
  gobject = createGiottoObject(
    expression = sobject@assays[[sobject@active.assay]]@data[, meta$cell_ID],
    spatial_locs = spatial_data[spatial_data$cell_ID %in% meta$cell_ID,],
    cell_metadata = meta,
    dimension_reduction = list(
      pca = sobject@reductions$pca@cell.embeddings[meta$cell_ID, ],
      umap = sobject@reductions$umap@cell.embeddings[meta$cell_ID, ])
  )
  
  gobject = createSpatialDelaunayNetwork(gobject)
  cpe_scores = cellProximityEnrichment(gobject, cluster_column = "CellType", # change cell resolution here 
                                      number_of_simulations = 1000, adjust_method = "fdr") # change number of simulations 
  
  cpe_scores$raw_sim_table$ID = samp
  cpe_scores$enrichm_res$ID = samp
  cpe_scores$raw_sim_table$Condition = meta$Condition[1]
  cpe_scores$enrichm_res$Condition = meta$Condition[1]
  
  return(cpe_scores)
}

# load data as Seurat object 
sobject = readRDS("data/Hypoxia_Cortex_Final.rds")
spatial_data = get_spat_data(sobject, names(sobject@images))

samples = unique(sobject$orig.ident)
results = lapply(samples, get_cp_enrichments, sobject, spatial_data)
names(results) = samples

# combine results across samples 
comb_raw_sim_table = do.call(rbind, lapply(results, function(x) x$raw_sim_table))
comb_enrichm_res = do.call(rbind, lapply(results, function(x) x$enrichm_res))

fisher_method = function(p_vals){
  p_vals[p_vals == 0] = .Machine$double.eps
  -2 * sum(log(p_vals))
}

# summarize results for each condition (NX, HX)
# - use Fisher method to combine p-values 
# TODO: PI_value is currently averaged, consider recalcualting from combined p-value
combine_results = function(condition){
  enrichm_res = comb_enrichm_res %>%
    filter(Condition == condition) %>%
    group_by(unified_int, type_int, Condition) %>%
    summarise(
      across(c(original, simulations, enrichm, PI_value, int_ranking), mean),
      across(c(p_higher_orig, p_lower_orig, p.adj_higher, p.adj_lower), fisher_method)) %>%
    ungroup() %>% as.data.table()
  
  raw_sim_table = comb_raw_sim_table %>%
    filter(Condition == condition) %>%
    group_by(unified_int, type_int, round, orig, Condition) %>%
    summarise(V1 = mean(V1, na.rm = TRUE)) %>%
    ungroup() %>% as.data.table()
  
  list(enrichm_res = enrichm_res, raw_sim_table = raw_sim_table)
}

# plotting 
for(condition in c("NX", "HX")){
  results = combine_results(condition)
  
  cellProximityBarplot(gobject = NULL, results, p_val = 0.05, show_plot = TRUE, return_plot = FALSE, save_plot = FALSE)
  cellProximityHeatmap(gobject = NULL, results, show_plot = TRUE, return_plot = FALSE, save_plot = FALSE) 
  cellProximityNetwork(gobject = NULL, results, show_plot = TRUE, return_plot = FALSE, save_plot = FALSE)
}



































