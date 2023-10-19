
# load libraries ----
library(dplyr)
library(ggplot2)
library(patchwork)
library(data.table)
library(svglite)
library(Seurat)
library(SeuratObject)
library(future)
options(future.globals.maxSize = 1024 * 1024 * 1024 * 8) # 8 GB

# devtools::install_github("sqjin/CellChat")
library(CellChat)
# devtools::install_github("Wei-BioMath/NeuronChat")
library(NeuronChat)

# devtools::install_github("drieslab/Giotto@suite")
library(Giotto) 
# devtools::install_github("RubD/GiottoUtils")
# devtools::install_github("RubD/GiottoClass")
# devtools::install_github("drieslab/GiottoData")
library(GiottoData) 
if(!checkGiottoEnvironment()){
  installGiottoEnvironment() # run once 
}

# functions ----

# for each sample 
# - get cell identifiers 
# - get x/y coordinates from image slot 
# - spatially translate by some constant to prevent overlap
get_spat_data = function(data, samples) {
  max_x = 10000
  max_y = 10000
  do.call(rbind, lapply(1:length(samples), function(i) {
    samp = samples[i]
    img = data@images[[samp]]
    x_offset = (i - 1) * max_x
    y_offset = (i - 1) * max_y
    data.frame(
      x = img@boundaries[["centroids"]]@coords[, 1] + x_offset,
      y = img@boundaries[["centroids"]]@coords[, 2] + y_offset,
      cell_ID = img@boundaries[["centroids"]]@cells,
      image_ID = samp
    )
  }))
}

run_giotto_analyses = function(cond, sobject, spatial_data, LR_data) {
  # get metadata and add cell ids as column
  meta = subset(sobject@meta.data, subset = Condition == cond)
  meta$cell_ID = rownames(meta)
  # create Giotto object
  gobject = createGiottoObject(
    expression = list(normalized = sobject@assays$integrated@data[, meta$cell_ID]),
    spatial_locs = spatial_data[spatial_data$cell_ID %in% meta$cell_ID,],
    cell_metadata = meta,
    dimension_reduction = list(
      pca = sobject@reductions$pca@cell.embeddings[meta$cell_IDD, ],
      umap = sobject@reductions$umap@cell.embeddings[meta$cell_ID, ])
  )
  # calculate Delaunay network 
  gobject = createSpatialDelaunayNetwork(gobject)
  # spatPlot(gobject, show_network = TRUE, point_size = 0)
  
  # get cell proximity enrichment scores 
  cpe_scores = cellProximityEnrichment(gobject, 
                                       cluster_column = "CellSubType", # change cell resolution here 
                                       number_of_simulations = 1000, 
                                       adjust_method = "fdr") # change number of simulations 
  cpe_scores$raw_sim_table$Condition = cond
  cpe_scores$enrichm_res$Condition = cond
  
  # filter LR data by expression data
  LR_data = LR_data[LR_data$ligand %in% gobject@feat_ID$rna & LR_data$receptor %in% gobject@feat_ID$rna,]
  # get statistical significance of LR expression changes upon cell-cell interaction
  spatial_all_scores = spatCellCellcom(gobject,
                                       cluster_column = 'CellSubType',
                                       random_iter = 100,
                                       feat_set_1 = LR_data$ligand,
                                       feat_set_2 = LR_data$receptor,
                                       adjust_method = 'fdr',
                                       do_parallel = TRUE, cores = 24, 
                                       verbose = "a lot")
  # get statistical significance of LR expression changes based on expression
  expr_only_scores = exprCellCellcom(gobject,
                                     cluster_column = 'CellSubType',
                                     random_iter = 100,
                                     feat_set_1 = LR_data$ligand,
                                     feat_set_2 = LR_data$receptor) 
  # combine spatial and expression based cell-cell communication LR scores 
  comb_comm = combCCcom(spatialCC = spatial_all_scores,
                        exprCC = expr_only_scores)
  # save
  res = list(cpe_scores, spatial_all_scores, expr_only_scores, comb_comm)
  names(res) = c("cpe_scores", "spatial_all_scores", "expr_only_scores", "comb_comm")
  return(res)
}

# script ----

# load data as Seurat object 
sobject = readRDS("data/Hypoxia_Cortex_Final.rds")
spatial_data = get_spat_data(sobject, names(sobject@images))

# get ligand-receptor pairs from CellChat
data(list='interactionDB_mouse'); neuron_chat_db = eval(parse(text = 'interactionDB_mouse'))
neuron_chat_db = purrr::map_dfr(neuron_chat_db, ~expand.grid(
  interaction_name = .x$interaction_name,
  ligand = .x$lig_contributor,
  receptor = .x$receptor_subunit))
# get ligand-receptor pairs from NeuronChat
cell_chat_db = CellChatDB.mouse$interaction[, c("interaction_name", "ligand", "receptor")]
LR_data = dplyr::distinct(rbind(cell_chat_db, neuron_chat_db))
rownames(LR_data) = NULL

# run
conditions = c("NX","HX")
results = lapply(conditions, run_giotto_analyses, sobject, spatial_data, LR_data)
names(results) = conditions
#saveRDS(results, file = "./output/results.rds")

# plotting ----

results = readRDS(file = "./output/results_cortex.rds")

for(cond in c("NX", "HX")){
  cpe_scores = results[[cond]]$cpe_scores
  
  png(filename = paste0("./figures/CPE_heatmap_",cond,".png"), width = 8.5, height = 8, units = "in", res = 300)
  cellProximityHeatmap(gobject = NULL, cpe_scores, show_plot = TRUE, return_plot = TRUE, save_plot = FALSE) 
  dev.off()
  
  cellProximityNetwork(gobject = NULL, cpe_scores, show_plot = TRUE, return_plot = TRUE, save_plot = FALSE)
  ggsave(paste0("./figures/CPE_network_",cond,".png"), device = "png", width = 11, height = 9, units = "in")
  
  spatial_all_scores = results[[cond]]$spatial_all_scores
  
  for(i in 2:3){
    if(i == 2) name = "spatial" else name = "expression"
    LR_scores = results[[cond]][[i]]
    LR_scores = na.omit(LR_scores)
    
    LR_scores_selected = LR_scores[p.adj <= 0.2 & abs(log2fc) > 0.1 & lig_nr >= 2 & rec_nr >= 2]
    data.table::setorder(LR_scores_selected, -PI)
    
    top_LR_ints = unique(LR_scores_selected[order(-abs(PI))]$LR_comb)[1:30]
    top_LR_cell_ints = unique(LR_scores_selected[order(-abs(PI))]$LR_cell_comb)[1:30]
    
    plotCCcomHeatmap(gobject = NULL,
                     LR_scores, 
                     selected_LR = top_LR_ints,
                     selected_cell_LR = top_LR_cell_ints,
                     aggl_method = 'ward.D2',
                     show = 'LR_expr',
                     show_plot = FALSE, return_plot = TRUE, save_plot = FALSE) +
      ggtitle(paste(cond, "cell-cell ligand-receptor by", name)) +
      ggplot2::scale_fill_gradientn(colours = c('white', 'red', 'darkred')) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(face = "bold"))
    ggsave(paste0("./figures/LR_heatmap_",name,"_",cond,".png"), device = "png", width = 6.5, height = 6.5, units = "in")
    
    plotCCcomDotplot(gobject,
                     LR_scores, 
                     selected_LR = top_LR_ints,
                     selected_cell_LR = top_LR_cell_ints,
                     aggl_method = 'ward.D2',
                     cluster_on = 'log2fc',
                     show_plot = FALSE, return_plot = TRUE, save_plot = FALSE) +
      ggtitle(paste(cond, "cell-cell ligand-receptor by", name)) +
      ggplot2::scale_color_gradientn(
        colours = c('darkblue', 'blue', 'white', 'red', 'darkred'), 
        limits = c(-5, 5),  breaks = c(-5, 0, 5)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(face = "bold"))
    ggsave(paste0("./figures/LR_spatial_dotplot_",cond,".png"), device = "png", width = 7.5, height = 7.5, units = "in")
    }
    
  
    comb_comm = results[[cond]]$comb_comm
    comb_comm = na.omit(comb_comm)
    
    # p1 = plotRecovery(gobject = NULL, comb_comm, ground_truth = "spatial", show_plot = FALSE, return_plot = TRUE, save_plot = FALSE)
    # p2 = plotRecovery(gobject = NULL, comb_comm, ground_truth = "expression", show_plot = FALSE, return_plot = TRUE, save_plot = FALSE)
    # p1 + p2
    
    comb_comm_selected = comb_comm[p.adj <= 0.2 & abs(log2fc) > 0.1 & lig_nr >= 2 & rec_nr >= 2]
    data.table::setorder(comb_comm_selected, -PI_spat)
    
    top_LR_ints = unique(comb_comm_selected[order(-abs(PI_spat))]$LR_comb)[1:20]
    top_LR_cell_ints = unique(comb_comm_selected[order(-abs(PI_spat))]$LR_cell_comb)[1:12]

    plotCombineCellCellCommunication(gobject = NULL,
                                     comb_comm,
                                     selected_LR = top_LR_ints,
                                     selected_cell_LR = top_LR_cell_ints,
                                     detail_plot = F,
                                     simple_plot = T,
                                     show_plot = TRUE, return_plot = TRUE, save_plot = FALSE)
    ggsave(paste0("./figures/LR_detailed_comparisons_",cond,".png"), device = "png", width = 12, height = 9, units = "in")
}













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


results = tmp
results = results[[1]]




































