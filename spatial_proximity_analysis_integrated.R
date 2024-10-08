# Load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(glmmTMB)
  library(broom.mixed)
  library(DHARMa)
  library(seriation)
  library(pheatmap)
  library(RColorBrewer)
  library(igraph)
  library(ggraph)
})

# Set up Python ---- 

library(reticulate)
# reticulate::install_python()
virtualenv_create("python_environment")
use_virtualenv("python_environment", required = TRUE)
sapply(c("pandas", "numpy", "scipy", "matplotlib", "seaborn"), function(pkg) 
  if(!py_module_available(pkg)) virtualenv_install("python_environment", pkg))

source_python("./python_functions_modified.py")

# R functions ##################################################################

#' Retrieve Spatial Data from a Seurat Object
#'
#' This function extracts spatial data from a Seurat object, combining 
#' information from images and metadata. It processes each image to extract 
#' coordinates and cell IDs, then merges this data with the metadata. The function 
#' allows for optional filtering by sample IDs and includes cell type data.
#'
#' @param seurat_object A Seurat object containing spatial data and metadata.
#' @param cell_type_column The name of the column in the metadata that contains
#' cell type information. This column name is modified for consistent formatting.
#' @param sample Optional vector of sample IDs for filtering; if provided, only 
#' cells from the specified samples are included in the output.
#'
#' @return A data frame containing spatial coordinates (x, y), cell IDs, 
#' sample IDs, image IDs, and cell types for each cell in the Seurat object. 
#' The data frame includes columns: `x`, `y`, `cell_ID`, `sample_ID`, 
#' `image_ID`, and the specified `cell_type_column`.
#'
#' @examples
#'
#' # Retrieve spatial data from a Seurat object with sample filtering and cell 
#' type information:
#' 
#' get_spatial_data(seurat_object, cell_type_column = "cell_type", 
#' sample = c("Sample1", "Sample2"))

get_spatial_data = function(seurat_object, cell_type_column, sample=NULL) {
  images = names(seurat_object@images)
  meta = seurat_object@meta.data %>% 
    mutate(cell_ID = row.names(.), sample_ID = orig.ident) 
  # mutate(!!cell_type_column := str_replace_all(
  #   !!sym(cell_type_column), " ", "")) 
  if (!is.null(sample)) {
    meta = filter(meta, sample_ID %in% sample)
  }
  spat_data = lapply(images, function(i) {
    data = seurat_object@images[[i]]
    df = data.frame(
      x = data@boundaries[["centroids"]]@coords[, 1],
      y = data@boundaries[["centroids"]]@coords[, 2],
      cell_ID = data@boundaries[["centroids"]]@cells
    ) %>% 
      inner_join(meta, by = "cell_ID") %>% 
      mutate(image_ID = i)}) %>% 
    bind_rows()
  print(spat_data[1:10, 1:5])
  return(spat_data)
}

#' Retrieve Expression Data from a Seurat Object
#'
#' This function extracts expression data from a Seurat object. It accesses the
#' scaled expression data from the integrated assay of the Seurat object and 
#' optionally filters it based on sample IDs. The function ensures that only 
#' expression data corresponding to the cells present in the specified samples 
#' (if any) are returned.
#'
#' @param seurat_object A Seurat object containing expression data within its
#' integrated assay.
#' @param sample Optional vector of sample IDs for filtering; if provided, only 
#' expression data from the specified samples are included in the output.
#'
#' @return A matrix of expression data with genes as rows and cells as columns. 
#' Only cells matching the provided sample IDs (if any) are included. The first 
#' few rows and columns of the matrix are displayed as a preview.
#'
#' @examples
#' # Retrieve expression data from a Seurat object without sample filtering
#' get_expression_data(seurat_object)
#'
#' # Retrieve expression data from a Seurat object with sample filtering
#' get_expression_data(seurat_object, sample = "NX1")

get_expression_data = function(seurat_object, sample=NULL) {
  meta = seurat_object@meta.data %>% 
    mutate(cell_ID = row.names(.), sample_ID = orig.ident)
  if (!is.null(sample)) {
    meta = filter(meta, sample_ID %in% sample)
  }
  expr_data = seurat_object@assays$integrated@scale.data %>%
    t() %>%
    as.data.frame() %>%
    filter(row.names(.) %in% meta$cell_ID)
  print(expr_data[1:10, 1:5])
  return(expr_data)
}

#' Retrieve Ligand-Receptor Pairs from NeuronChat and CellChat Databases
#'
#' This function consolidates ligand-receptor interaction data from two distinct
#' sources: NeuronChat and CellChat. It processes data from the interaction
#' databases of both packages, formats them into a uniform structure, 
#' and combines them into a single data frame for further analysis.
#'
#' @return A data frame containing consolidated ligand-receptor pairs with 
#' the following columns:
#'   - interaction_name: Name of the interaction.
#'   - ligand: Name of the ligand involved in the interaction.
#'   - receptor: Name of the receptor involved in the interaction.
#'   - DB: Source database of the interaction ('neuronchat' or 'cellchat').
#'
#' @examples
#' LR_pairs = get_LR_pairs()

get_LR_pairs = function() {
  # devtools::install_github("sqjin/CellChat")
  library(CellChat)
  # devtools::install_github("Wei-BioMath/NeuronChat")
  library(NeuronChat)
  
  data(list='interactionDB_mouse')
  neuron_chat_db = eval(parse(text = 'interactionDB_mouse'))
  neuron_chat_db = purrr::map_dfr(neuron_chat_db, ~expand.grid(
    interaction_name = .x$interaction_name,
    ligand = .x$lig_contributor,
    receptor = .x$receptor_subunit)) %>%
    mutate(DB = "neuronchat")
  
  cell_chat_db = CellChatDB.mouse$interaction %>%
    dplyr::select(interaction_name_2) %>%
    rename(interaction_name = interaction_name_2) %>%
    separate(interaction_name, into = c("ligand", "receptor"), sep = " - ", remove = F) %>%
    mutate(receptor = str_replace_all(receptor, "[()]", "")) %>%
    separate_rows(receptor, sep = "\\+") %>%
    mutate(DB = "cellchat")
  
  LR_pairs = rbind(cell_chat_db, neuron_chat_db) %>%
    distinct(ligand, receptor, .keep_all = TRUE) %>%
    mutate(across(where(is.character), trimws)) 
  
  rownames(LR_pairs) = NULL
  print(LR_pairs[1:10, 1:4])
  return(LR_pairs)
}

#' Calculate differential spatial statistics for cell types
#'between two conditions.
#'
#' This function evaluates spatial differences between cell types in different
#' conditions using a generalized linear mixed-effects model (GLMM). It filters
#' and prepares data based on specified thresholds, and then applies the GLMM
#' to assess the impact of conditions on the spatial distribution of cell types.
#'
#' @param spat_stats A data frame containing spatial statistics, which is the 
#' output of 'get_spatial_stats'. 
#' @param ref_level A character string specifying the reference level for the condition groups.
#' @param comp_level A character string specifying the comparison level for the condition groups.
#' @param min_nonzero_threshold A numeric value representing the minimum 
#' threshold for the non-zero count ratio for each cell type pair and sample. 
#' Default is 0.10.
#' @param min_data_count_threshold An integer indicating the minimum number 
#' of measures per cell type pair and sample. Default is 10.
#' @param model_diagnostics A logical value indicating whether to perform and display
#' diagnostic tests for the model. Default is TRUE. For interpretation of model diagnostics see:
#' https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html 
#'
#' @return A data frame with computed model fits, including coefficients, adjusted p-values, 
#' and various diagnostic values if model diagnostics are enabled.
#'
#' @examples
#' spat_stats = get_spatial_stats(data)
#' get_spatial_differences(spat_stats, ref_level = "control")
#'
#' @seealso
#' get_spatial_stats for generating the input spatial statistics data frame.
#' plot_diff_heatmap for visualizing the differential spatial statistics.


get_spatial_differences = function(spat_stats,
                                   min_nonzero_threshold = 0.10,
                                   min_data_count_threshold = 10,
                                   model_diagnostics = TRUE) {
  valid_pairs = spat_stats %>%
    group_by(condition, cell_type_pair, sample_ID) %>%
    summarise(non_zero_count = sum(b_ratio != 0, na.rm = TRUE),
              total_count = n()) %>%
    group_by(cell_type_pair) %>%
    filter(all(non_zero_count / total_count >= min_nonzero_threshold,
               total_count > min_data_count_threshold)) %>%
    pull(cell_type_pair) %>%
    unique()
  if (length(valid_pairs) == 0) {
    stop("Error: No cell type pairs have sufficient data based on the specified thresholds.")
  }
  
  spat_stats %>%
    filter(cell_type_pair %in% valid_pairs) %>%
    mutate(condition = as.numeric(condition), sample_ID = factor(sample_ID)) %>%
    group_by(cell_type_a, cell_type_b) %>%
    do({
      model = glmmTMB(b_ratio ~ condition + (1 | sample_ID),
                      family = gaussian(), data = .)
      model_fit = broom.mixed::tidy(model)
      if(model_diagnostics){
        simulation_output = simulateResiduals(fittedModel = model, plot = F)
        print(plot(simulation_output))
        model_fit$uniformity_pval = testUniformity(simulation_output, plot = F)$p.value
        model_fit$outlier_pval = testOutliers(simulation_output, plot = F)$p.value
        model_fit$dispersion_pval = testDispersion(simulation_output, plot = F)$p.value
        model_fit$zeroInflation_pval = testZeroInflation(simulation_output, plot = F)$p.value
      }
      model_fit
    }) %>%
    ungroup() %>%
    filter(term == "condition") %>%
    mutate(padj = p.adjust(p.value, "BH"),
           cell_type_pair = paste(cell_type_a, cell_type_b, sep = " -- "))
}


# Example usage
# spat_stats should be your dataframe with continuous condition values ranging from 1 to 10


#' Heatmaps for spatial statistics data
#'
#' This function creates heatmaps for spatial statistics data, highlighting 
#' the average ratios between cell types under different conditions. Note that
#' the rows are 'cell_type_a' and the columns are query 'cell_type_b'. 
#' There is an option to reorder the data for optimal clustering.
#'
#' @param spat_stats A data frame that is the output of 'get_spatial_stats'.  
#' @param opt_order A logical value indicating whether to apply optimal ordering 
#' (clustering) to the heatmap. Default is FALSE.
#'
#' @return Heatmaps for each condition in the spatial statistics data. 
#' Each heatmap displays the average b_ratio between cell types.
#'
#' @examples
#' plot_stats_heatmap(spat_stats, TRUE)

plot_stats_heatmap = function(spat_stats, opt_order = FALSE) {
  for (c in unique(spat_stats$condition)) {
    mat = spat_stats %>%
      dplyr::select(cell_type_a, cell_type_b, b_ratio, condition) %>%
      filter(condition == c) %>%
      group_by(cell_type_a, cell_type_b) %>%
      summarise(avg_ratio = mean(b_ratio, na.rm = TRUE)) %>%
      pivot_wider(
        names_from = cell_type_b, 
        values_from = avg_ratio, 
        values_fill = NA
      ) %>%
      column_to_rownames("cell_type_a") %>%
      as.matrix()
    mat = mat[, match(rownames(mat), colnames(mat))]
    mat[mat == 0] = NA
    if (opt_order) {
      o1 = seriate(dist(mat), method = "OLO")[[1]]$order
      o2 = seriate(dist(t(mat)), method = "OLO")[[1]]$order
      mat = mat[o1, o2]
    }
    log_breaks = function(n, minval, maxval) {
      breaks = exp(seq(log(minval + 1), log(maxval + 1), length.out = n)) - 1
      return(breaks)
    }
    breaks = log_breaks(n = 51, minval = 0, maxval = 1)
    main = paste("cell type proximity", c, "\nmean_dmin:", 
                 round(mean(unique(spat_stats$d_min))),
                 "; mean_dmax:", round(mean(unique(spat_stats$d_max))))
    pheatmap(mat, 
             main = main, 
             na_col = "lightgrey",
             color = colorRampPalette(c("white","#813e8f","#813e8f", "#813e8f", "#4b0857"))(50), 
             breaks = breaks,
             display_numbers = FALSE, border_color = NA, 
             cluster_rows = FALSE, cluster_cols = FALSE)
  }
}

#' Heatmaps for differential cell proximity analysis
#'
#' This function generates heatmaps to visualize the differences in cell 
#' proximity for various cell types between conditions, highlighting significant 
#' differences based on a specified significance level. Note that the rows are
#' cell_type_a' and the columns are query 'cell_type_b'. 
#'
#' @param spat_diff A data frame that is the output of 'get_spatial_differences'.
#' @param opt_order A logical value indicating whether to apply optimal ordering
#' (clustering) to the heatmap. Default is FALSE.
#' @param sig_level A numeric value specifying the significance level for 
#' highlighting in the heatmap. Default is 0.01.
#'
#' @return Heatmaps representing differential cell proximity, where the color 
#' scale represents estimated effect of condition, with significant differences 
#' marked "*".
#'
#' @examples
#' plot_diff_heatmap(spat_diff, TRUE, 0.05)

plot_diff_heatmap = function(spat_diff, opt_order = FALSE, sig_level = 0.05, main = "") {
  
  all_cell_types = unique(c(as.character(spat_diff$cell_type_a), 
                            as.character(spat_diff$cell_type_b)))
  complete_grid = expand.grid(cell_type_a = all_cell_types, 
                              cell_type_b = all_cell_types)
  
  get_matrix = function(data, col, scale = FALSE){
    if (scale) data[[col]] = as.numeric(scale(data[[col]]))
    mat = data %>%
      dplyr::select(cell_type_a, cell_type_b, col) %>%
      bind_rows(
        complete_grid %>%
          anti_join(spat_diff, by = c("cell_type_a", "cell_type_b")) %>%
          mutate(!!col := NA)
      ) %>%
      pivot_wider(names_from = cell_type_b, values_from = !!col) %>%
      column_to_rownames("cell_type_a") %>%
      as.matrix()
    mat = mat[, match(rownames(mat), colnames(mat))]
  }
  mat_beta = get_matrix(spat_diff, "estimate")
  mat_sig = get_matrix(spat_diff, "padj")
  mat_sig[mat_sig < sig_level] = "*"
  mat_sig[mat_sig != "*"] = ""
  mat_sig[is.na(mat_sig)] = ""
  if (opt_order) {
    mat_tmp = mat_beta
    mat_tmp[is.na(mat_tmp)] = 0
    o1 = seriate(dist(mat_tmp), method = "OLO")[[1]]$order
    o2 = seriate(dist(t(mat_tmp)), method = "OLO")[[1]]$order
    mat_beta = mat_beta[o1, o2]
    mat_sig = mat_sig[o1, o2]
  }
  max_abs = max(abs(mat_beta), na.rm = TRUE)
  pheatmap(mat_beta, 
           display_numbers = mat_sig,
           main = main, 
           na_col = "lightgrey",
           number_color = "white",
           fontsize_number = 20,
           color = colorRampPalette(
             c("#156a2f","#66b66b","#66b66b","white","#813e8f","#813e8f","#4b0857"))(51), 
           breaks = seq(-max_abs, max_abs, length.out = 52),
           border_color = NA, cluster_rows = F, cluster_cols = F)
}

#' Bar plots for differential cell proximity
#'
#' This function creates bar plots to display the top differential spatial
#' statistics, highlighting significant differences and indicating the 
#' direction of change.
#'
#' @param spat_diff A data frame that is output of 'get_spatial_differences'.
#' @param n An integer specifying the number of top differences for positive
#' and negative estimates to display in the bar plot. Default is 20.
#' @param sig_level A numeric value specifying the significance level for 
#' highlighting in the plot. Default is 0.01.
#'
#' @return A ggplot object representing a bar plot of the top n
#' differential statistics, with significant differences marked and 
#' direction of change indicated.
#'
#' @examples
#' plot_diff_barplot(spat_diff, 20, 0.01)

plot_diff_barplot = function(spat_diff, n = 20, sig_level = 0.01) {
  spat_diff %>%
    mutate(
      significance = ifelse(padj < sig_level, "*", ""),
      y_pos = ifelse(estimate >= 0, 
                     estimate + std.error + 0.05, 
                     estimate - std.error - 0.05)
    ) %>%
    group_by(sign = estimate >= 0) %>%
    top_n(n, if_else(sign, estimate, -estimate)) %>%
    ungroup() %>%
    arrange(estimate) %>%
    ggplot(aes(fct_reorder(cell_type_pair, estimate), 
               estimate, fill = estimate > 0)) +
    geom_col() +
    geom_text(aes(label = significance, y = y_pos), 
              vjust = 0.5, hjust = 0.5, color = "black", size = 5) +
    geom_errorbar(aes(ymin = estimate - std.error, 
                      ymax = estimate + std.error), width = 0.2) +
    coord_flip() +
    scale_fill_manual(values = c("#66b66b", "#813e8f"), 
                      labels = c("decrease", "increase")) +
    labs(x = "Cell Type Pair", y = "Differential Cell Density Change", fill = "") +
    theme_minimal() +
    theme(legend.position = "top")
}

#' Network plots for differential cell proximity
#'
#' This function generates network plots to visualize the relationships between 
#' cell types based on their differential proximity statistics. 
#' It highlights significant relationships and allows the inclusion of 
#' insignificant ones as an option.
#'
#' @param spat_diff A data frame that is output of 'get_spatial_differences'.
#' @param sig_level A numeric value specifying the significance level for 
#' highlighting edges in the network plot. Default is 0.05.
#' @param plot_insignificant A logical value indicating whether to include 
#' insignificant relationships in the plot. Default is FALSE.
#'
#' @return A ggraph object representing a network plot, with edges weighted and 
#' colored based on significance and direction of change.
#'
#' @examples
#' plot_diff_network(spat_diff, 0.05, TRUE)

plot_diff_network = function(spat_diff, sig_level = 0.05, plot_insignificant = FALSE, 
                             show_arrows = TRUE) {
  edges_to_plot = if(plot_insignificant) spat_diff 
  else spat_diff[spat_diff$padj < sig_level, ]
  
  igraph_network = graph_from_data_frame(
    d = edges_to_plot[, c('cell_type_a', 'cell_type_b')], directed = TRUE)
  
  E(igraph_network)$weight = abs(edges_to_plot$estimate)
  E(igraph_network)$alpha = scales::rescale(
    edges_to_plot$estimate, to = c(0.1, 1))
  
  E(igraph_network)$color = ifelse(
    edges_to_plot$padj < sig_level, 
    ifelse(edges_to_plot$estimate > 0, 'positive', 'negative'), 
    'insignificant')
  
  node_colors = setNames(
    colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(
      length(unique(c(edges_to_plot$cell_type_a, 
                      edges_to_plot$cell_type_b)))
    ), unique(c(edges_to_plot$cell_type_a, edges_to_plot$cell_type_b))
  )
  
  # For other options see: https://igraph.org/r/doc/layout_.html
  layout = layout_in_circle(igraph_network)
  arrow_type = if(show_arrows) arrow(length = unit(4, 'mm'), type = 'closed', ends = 'first') else NULL
  
  ggraph(igraph_network, layout = layout) +
    geom_node_point(aes(color = name), size = 5, show.legend = FALSE) +
    geom_edge_link(
      aes(edge_width = weight, color = color, alpha = alpha), 
      arrow = arrow_type
    ) +
    geom_edge_loop(
      aes(color = color, edge_width = weight, alpha = alpha, strength = 0.5)
    ) +
    scale_edge_width_continuous(range = c(0.1, 3)) +
    scale_edge_alpha_continuous(range = c(0.1, 1)) + 
    geom_node_text(aes(label = name), repel = TRUE, size = 5, show.legend = FALSE) +
    scale_edge_color_manual(
      values = c('positive' = '#813e8f', 'negative' = '#66b66b', 
                 'insignificant' = if(plot_insignificant) 'gray95' else NA)
    ) +
    theme_void() +
    guides(color = guide_legend())
}

#' Visualize spatial exemplars to highlight differential spatial statistics.
#'
#' Spatial plots indicating the distribution of select center and surround 
#' cell types between conditions. An additional violin plot shows the distribution
#' of the b_ratios of cell type b for every cell of type a in the selected samples.
#' The function will automatically select two samples from each condition that 
#' have the most comparable cell densities.  
#' 
#' This function is not fully tested. It is common to see poorly representative
#' spatial plots given the sparsity of some cell types in MERFISH data.
#' 
#' @param spat_stats A data frame containing spatial statistics, which is the output
#' of 'get_spatial_stats'.
#' @param spat_diff A data frame containing differential spatial statistics, which
#' is the output of 'get_spatial_differences'.
#' @param cell_type_column A character string specifying the column in 'spat_stats'
#' that contains cell type information.
#' @param ref_level A character string specifying the reference level for the condition groups.
#' @param comp_level A character string specifying the comparison level for the condition groups.
#' @param cell_type_a A character string specifying the first cell type of interest.
#' @param cell_type_b A character string specifying the second cell type of interest.
#' @param min_b_count A numeric value representing the minimum count of 'b' cells
#' required for a cell to be considered in analysis. Default is 0.
#' @param min_all_count A numeric value representing the minimum total count of cells
#' required for a cell to be considered in analysis. Default is 0.
#' @param sample_area A numeric value specifying the side length of the square area
#' around the central cell for visualization. Default is 120.
#' @param seed An optional integer value to set the seed for random number generation,
#' ensuring reproducibility. Default is NULL.
#'
#' @return A list of ggplot objects for each selected sample, illustrating the spatial
#' distribution of cell types within the specified area and a violin and boxplot
#' comparison of log-transformed counts between conditions.
#'
#' @examples
#' spat_stats = get_spatial_stats(data)
#' spat_diff = get_spatial_differences(spat_stats, ref_level = "control")
#' plot_spatial_exemplars(spat_stats, spat_diff, "cell_type", "T_cells", "B_cells")
#'
#' @seealso
#' get_spatial_stats and get_spatial_differences for required input data frames.
#' ggplot for the underlying plotting system used.

plot_spatial_exemplars = function(spat_stats, 
                                  spat_diff) {
  
  p3 = spat_stats %>%
    ggplot(aes(x = condition, y = log(b_ratio))) +
    geom_point(aes(color = as.factor(condition)), size = 3) +
    labs(x = "Condition", y = "Log B Ratio (Surround Cell)") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.line = element_line(linewidth = 1)
    )
  p3
} 


# SEA_AD SPATIAL ANALYSIS #####################################################

# Packages
library(anndata)
library(patchwork)
library(ggpubr)


# SEATTLE ALZHEIMER'S BRAIN ATLAS MERFISH DATASET PROCESSING


# Read in spatial coordinates and metadata
hfivead <- read_h5ad("./data/raw/SEAAD_MTG_MERFISH_all-nuclei.2023-05-08 (1).h5ad")
obsm <- hfivead$obsm
hfivead_metadata <- hfivead$obs
hfivead_metadata_obsm <- cbind(hfivead_metadata, obsm)

# Filter unused cells
hfivead_metadata_obsm_filtered <- hfivead_metadata_obsm %>%
  filter(`Used in analysis` == 1)

# Create and format data frame required for local cell type proportion analysis
row_names <- rownames(hfivead_metadata_obsm_filtered)
hfivead_metadata_obsm_filtered$cell_ID <- row_names
hfivead_subset <- hfivead_metadata_obsm_filtered[, c("X_selected_cell_spatial_tiled.1", "X_selected_cell_spatial_tiled.2", "Section", "Class", "Subclass", "Supertype", "Donor ID", "Continuous Pseudo-progression Score", "cell_ID", "Sex")]
names(hfivead_subset) <- c("x", "y", "orig.ident", "Class", "Subclass", "Supertype", "Donor ID", "condition", "cell_ID", "Sex")

samples = unique(hfivead_subset$`orig.ident`)[-59]
samples <- as.character(samples)
hfivead_subset$Subclass <- as.character(hfivead_subset$Subclass)
hfivead_subset$orig.ident <- as.character(hfivead_subset$orig.ident)

# Binarize condition (AD or no_AD) using mean CPS
average_value <- mean(hfivead_subset$condition)
hfivead_subset$binarized_condition <- ifelse(hfivead_subset$condition < average_value, "no_AD", "AD")


# INTEGRATED GLOBAL AND LOCAL CELL TYPE PROPORTION ANALYSIS

# Cell subclass proportions
total_counts <- hfivead_subset %>%
  group_by(orig.ident, condition, Sex) %>%
  summarise(total_count = n(), .groups = 'drop')

subclass_counts <- hfivead_subset %>%
  group_by(orig.ident, Subclass) %>%
  summarise(subclass_count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = Subclass, values_from = subclass_count, values_fill = list(subclass_count = 0))

subclass_proportions <- total_counts %>%
  left_join(subclass_counts, by = "orig.ident")

subclass_individual_proportions <- setdiff(names(subclass_proportions), c("orig.ident", "condition", "Sex", "total_count"))

subclass_proportions <- subclass_proportions %>%
  mutate(across(all_of(subclass_individual_proportions), ~ . / total_count, .names = "proportion_{.col}"))

# Get spatial stats
spat_stats = bind_rows(
  lapply(samples, function(s) {
    spat_data = hfivead_subset %>%
      filter(orig.ident==s)
    spat_data = r_to_py(spat_data)
    # See 'python_functions.py' for 'get_spatial_stats' documentation
    stats = get_spatial_stats(spat_data, proportion_data = subclass_proportions, cell_type_column = "Subclass", 
                              d_min_scale = 0, d_max_scale = 10) %>%
      mutate(cell_type_pair = paste(cell_type_a, cell_type_b, sep = " -- "))
    return(stats)
  })
)

# Format object types of spat_stats
spat_stats$cell_type_a <- as.character(spat_stats$cell_type_a)

# Spatial stats heatmap
plot_stats_heatmap(spat_stats, opt_order = F)

# Get spatial differences
spat_diff = get_spatial_differences(
  spat_stats,
  min_nonzero_threshold = 0.10, min_data_count_threshold = 10, model_diagnostics = F
)

# Spatial differences heatmap
plot_diff_heatmap(spat_diff, opt_order = F, sig_level = 0.1, main = "SEA_AD -- Alzheimer's")

# Spatial differences barplot
plot_diff_barplot(spat_diff, n = 10, sig_level = 0.1)

# Spatial differences network plot
plot_diff_network(spat_diff, sig_level = 0.1, plot_insignificant = T, show_arrows = T)
