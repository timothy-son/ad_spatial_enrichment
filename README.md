# BCB330 Global and Local Proportion Spatial Proximity Analysis for Alzheimer's Patients

## Overview
This repository contains the analysis scripts and documentation for a research project focused on the spatial changes associated with Alzheimer's disease (AD). Our study explores the alterations in global cell type proportions, sex differences, and local neighborhood proportions in AD using bioinformatics techniques.

## Project Goals
- Characterize Global Cell Type Proportions: Investigate changes in global cell type proportions in Alzheimer's disease compared to healthy controls.
- Explore Sex Differences: Assess how Alzheimer's disease affects male and female global cell type proportions.
- Local Proportion Changes: Analyze local cell type changes in Alzheimer's disease, identifying local cellular changes, particularly in the context of disease progression.

## Dataset
Raw data was obtained from [SEA_AD](https://portal.brain-map.org/explore/seattle-alzheimers-disease). The SEA AD MERFISH dataset of the middle temporal gyrus was used for the analyses.
- [MERFISH data](https://sea-ad-spatial-transcriptomics.s3.amazonaws.com/index.html#middle-temporal-gyrus/all_donors-h5ad/)
- [Cell pose data](https://sea-ad-spatial-transcriptomics.s3.amazonaws.com/index.html#middle-temporal-gyrus/H21.33.001/1217500590/)

## Files
- **spatial_proximity_analysis**: global and local cell type proportion analysis
- **spatial_proximity_analysis_continuous**: local spatial proximity analysis using continuous scale
- **spatial_proximity_analysis_integrated**: local spatial proximity analysis taking into account global proportions
- **python_functions**: python functions for spatial_proximity_analysis and spatial_proximity_analysis_continuous
- **python_functions_integrated**: python functions for spatial_proximity_analysis_integrated
- **sea_ad_merfish_tiled_plots**: spatial plots for all SEA-AD MERFISH samples
- **cell_pose_data_clustering**: SEA-AD cell pose data clustering by gene expression
