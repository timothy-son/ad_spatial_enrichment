
# Load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(seriation)
  library(pheatmap)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(RColorBrewer)
  library(igraph)
  library(ggraph)
})

# Set up Python ---- 

library(reticulate)
# reticulate::install_python()
virtualenv_create("python_environment")
use_virtualenv("python_environment", required = TRUE)
sapply(c("pandas", "numpy", "scipy"), function(pkg) 
  if(!py_module_available(pkg)) virtualenv_install("python_environment", pkg))

source_python("./python_functions.py")

# R functions ##################################################################

#' Retrieve spatial data from a Seurat Object
#'
#' This function extracts spatial data from a Seurat object, combining 
#' information from images and metadata. It processes each image to extract 
#' coordinates and cell IDs, and then merges this data with the metadata.
#'
#' @param seurat_object A Seurat object containing spatial data and metadata, 
#' where subject/sample ID is 'orig.ident' 
#'
#' @return A data frame containing spatial coordinates (x, y), cell IDs, 
#' sample IDs, and image IDs for each cell in the Seurat object.
#'
#' @examples
#' get_spatial_data(seurat_object)

get_spatial_data = function(seurat_object) {
  images = names(seurat_object@images)
  meta = seurat_object@meta.data %>% 
    mutate(cell_ID = row.names(.), sample_ID = orig.ident)
  
  lapply(images, function(i) {
    data = seurat_object@images[[i]]
    df = data.frame(
      x = data@boundaries[["centroids"]]@coords[, 1],
      y = data@boundaries[["centroids"]]@coords[, 2],
      cell_ID = data@boundaries[["centroids"]]@cells
      ) %>% 
      inner_join(meta, by = "cell_ID") %>% 
      mutate(image_ID = i)}) %>% 
    bind_rows()
}

#' Calculate differential spatial statistics for cell types between two 
#' conditions. 
#'
#' This function evaluates spatial differences between cell types using a 
#' generalized linear mixed-effects model (GLMM). 
#'
#' @param spat_stats A data frame containing spatial statistics, which is the 
#' output of 'get_spatial_stats'. 
#' @param min_nonzero_threshold A numeric value representing the minimum 
#' threshold for the non-zero count ratio for each cell type pair and sample. 
#' Default is 0.40.
#' @param min_data_count_threshold An integer indicating the minimum number 
#' of measure per cell type pair and sample. Default is 10.
#' @param ref_level The reference level for the condition groups.
#'
#' @return A data frame with computed model fits, adjusted p-values, and
#' flags for singular fits.
#' 
#' @examples
#' get_spatial_differences(spat_stats, 0.40, 10, "reference_level")

get_spatial_differences = function(spat_stats, 
                                   min_nonzero_threshold = 0.40, 
                                   min_data_count_threshold = 10, 
                                   ref_level) {
  valid_pairs = spat_stats %>%
    group_by(condition, cell_type_pair, sample_ID) %>%
    summarise(non_zero_count = sum(b_ratio != 0, na.rm = TRUE),
              total_count = n()) %>%
    group_by(condition, cell_type_pair) %>%
    filter(all(non_zero_count / total_count >= min_nonzero_threshold,
               total_count > min_data_count_threshold)) %>%
    pull(cell_type_pair) %>% 
    unique()
  
  spat_stats %>%
    filter(cell_type_pair %in% valid_pairs) %>%
    drop_na() %>%
    mutate(condition = factor(condition), sample_ID = factor(sample_ID)) %>%
    mutate(condition = relevel(condition, ref = ref_level)) %>%
    group_by(cell_type_a, cell_type_b) %>%
    do({
      model = glmer(b_ratio + 0.01 ~ condition + (1 | sample_ID), 
                    family = Gamma(link = "log"), data = .)
      model_fit = broom.mixed::tidy(model)
      singular_fit = isSingular(model)
      if (!is.null(model_fit)) {
        model_fit$cell_type_pair = paste(
          unique(.$cell_type_a), unique(.$cell_type_b), sep = "-"
          )
        model_fit$singular_fit = singular_fit
      }
      model_fit
    }) %>%
    ungroup() %>%
    filter(!singular_fit, str_detect(term, "condition")) %>%
    mutate(term = "condition", padj = p.adjust(p.value, "BH"))
}

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
  data = spat_stats[, c("cell_type_a", "cell_type_b", "b_ratio", "condition")]
  for (c in unique(data$condition)) {
    mat = data %>%
      filter(condition == c) %>%
      group_by(cell_type_a, cell_type_b) %>%
      summarise(avg_ratio = mean(b_ratio, na.rm = TRUE)) %>%
      pivot_wider(
        names_from = cell_type_b, 
        values_from = avg_ratio, 
        values_fill = list(avg_ratio = 0)
      ) %>%
      column_to_rownames("cell_type_a") %>%
      as.matrix()
    mat = mat[, match(rownames(mat), colnames(mat))]
    
    if (opt_order) {
      order = seriate(dist(mat), method = "OLO")[[1]]$order
      mat = mat[order, order]
    }
    main = paste("cell type proximity", c, "\nmean_dmin:", 
                 round(mean(unique(spat_stats$d_min))),
                 "; mean_dmax:", round(mean(unique(spat_stats$d_max))))
    pheatmap(mat, 
             main = main, 
             color = colorRampPalette(c("white", "red", "red", "darkred"))(50), 
             breaks = seq(0, 1, length.out = 51),
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

plot_diff_heatmap = function(spat_diff, opt_order = FALSE, sig_level = 0.05) {
  
  all_cell_types = unique(c(spat_diff$cell_type_a, spat_diff$cell_type_b))
  complete_grid = expand.grid(cell_type_a = all_cell_types, 
                              cell_type_b = all_cell_types)
  
  get_matrix = function(data, col, scale = FALSE){
    if (scale) data[[col]] = as.numeric(scale(data[[col]]))
    mat = data %>%
      select(cell_type_a, cell_type_b, col) %>%
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
  mat_beta = get_matrix(spat_diff, "estimate", scale = F)
  mat_sig = get_matrix(spat_diff, "padj")
  mat_sig[mat_sig < sig_level] = "*"
  mat_sig[mat_sig != "*"] = ""
  mat_sig[is.na(mat_sig)] = ""
  
  if(opt_order == TRUE){
    opt_order = seriate(dist(mat_beta), method = "OLO")[[1]]$order
    mat_beta = mat_beta[opt_order, opt_order]
    mat_sig = mat_sig[opt_order, opt_order]
  }
  main = paste("differential cell type proximity\n")
  max_abs = max(abs(mat_beta), na.rm = TRUE)
  pheatmap(mat_beta, 
           display_numbers = mat_sig,
           main = main, 
           na_col = "white",
           number_color = "white",
           fontsize_number = 20,
           color = colorRampPalette(c("darkblue","blue","white","red","darkred"))(51), 
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
    scale_fill_manual(values = c("blue", "red"), 
                      labels = c("decrease", "increase")) +
    labs(x = "Cell Type Pair", y = "Estimate", fill = "Estimate Sign") +
    theme_minimal()
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

plot_diff_network = function(spat_diff, sig_level = 0.05, plot_insignificant = FALSE) {
  
  edges_to_plot = if(plot_insignificant) {
    spat_diff
  } else {
    spat_diff[spat_diff$padj < sig_level, ]
  }
  
  igraph_network = graph_from_data_frame(
    d = edges_to_plot[, c('cell_type_a', 'cell_type_b')], directed = TRUE)
  
  E(igraph_network)$weight = abs(edges_to_plot$estimate)
  E(igraph_network)$alpha = scales::rescale(edges_to_plot$estimate, to = c(0.1, 1))
  E(igraph_network)$color = ifelse(
    edges_to_plot$padj < sig_level,
    ifelse(edges_to_plot$estimate > 0, 'positive', 'negative'),
    'insignificant')
  
  node_colors = setNames(
    colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(
      length(unique(c(spat_diff$cell_type_a, spat_diff$cell_type_b)))
    ),
    unique(c(spat_diff$cell_type_a, spat_diff$cell_type_b))
  )
  # For other options see: https://igraph.org/r/doc/layout_.html
  layout = layout_in_circle(igraph_network)
  
  ggraph(igraph_network, layout = layout) +
    geom_edge_link(aes(edge_width = weight, color = color, alpha = alpha)) +
    geom_edge_loop(aes(color = color, edge_width = weight, alpha = alpha, strength = 0.5)) +
    scale_edge_width_continuous(range = c(0.1, 3)) +
    scale_edge_alpha_continuous(range = c(0.1, 1)) + 
    geom_node_point(aes(color = name), size = 5, show.legend = FALSE) +
    geom_node_text(aes(label = name), repel = TRUE, size = 5, show.legend = FALSE) +
    scale_edge_color_manual(
      values = c('positive' = 'red', 'negative' = 'blue',
                 'insignificant' = if(plot_insignificant) 'gray95' else NA)) +
    theme_void() +
    guides(color = guide_legend())
}

# script #######################################################################

# Load Seurat object
sobject = readRDS("data/Hypoxia_Cortex_Final.rds")
sobject = readRDS("data/SVZ_neuronal_subset_Seurat_FINAL.rds")
sobject = readRDS("data/CC_FINAL_dims10Seurat_Annotated.rds")

# Get spatial coordinates and clean cell type lables
spat_data = get_spatial_data(sobject) %>%
  mutate(across(c("CellType", "CellSubType"), ~str_replace_all(.x, " - ", "_") %>%
                  str_replace_all(" ", "_")))

# For each sample, get cell type proximities 
samples = unique(spat_data$sample_ID)
spat_stats = bind_rows(
  lapply(samples, function(s) {
    data = subset(spat_data, sample_ID == s)
    data = r_to_py(data)
    # See 'python_functions.py' for 'get_spatial_stats' documentation
    # Selection of d_max_scale is important
    stats = get_spatial_stats(data, "CellSubType", d_min_scale = 0, d_max_scale = 10)
    stats$cell_type_pair = paste(stats$cell_type_a, stats$cell_type_b, sep = "-")
    return(stats)
  })
)
# Plot cell type proximity (b_ratio) heatmaps for each condition 
plot_stats_heatmap(spat_stats)

# Get differential proximities (estimate) between conditions
# Singular fits may occur and are removed from the results 
# Adjust thresholds if too many cell types are removed 
spat_diff = get_spatial_differences(spat_stats, min_nonzero_threshold = 0.20, 
                                    min_data_count_threshold = 10, ref_level = "NX")

# Differential cell type proximity plots 
plot_diff_heatmap(spat_diff, opt_order = F, sig_level = 0.05)
plot_diff_barplot(spat_diff, n = 20, sig_level = 0.05)
plot_diff_network(spat_diff, sig_level = 0.05, plot_insignificant = T)

################################################################################





