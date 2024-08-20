import pandas as pd
import numpy as np
import gc
from scipy.spatial import KDTree
from scipy.sparse import coo_array
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

"""Calculate spatial statistics for cell types within a spatial dataset.

This function computes the proximity ratio for each cell to all cell types within specified distance scales.
It uses a KDTree for efficient spatial querying and returns a DataFrame with the calculated statistics.
It adjusts for differences in scale between sample images using the median smallest distance between all pairs of cells ('d_scale'). 

Parameters:
- spatial_data (pd.DataFrame): A DataFrame containing the spatial coordinates (x, y), cell types, and other metadata. Must also include a 'condition' column. 
- cell_type_column (str): The column name in `spatial_data` that contains the cell type information.
- d_min_scale (float, optional): The minimum distance scale factor for proximity calculations, defaulting to 0.
  If 'd_min_scale' > 0, returns the proximity ratio within an annulus defined 'd_min_scale' and 'd_max_scale'.
- d_max_scale (float, optional): The maximum distance scale factor for proximity calculations, defaulting to 20.

The function returns a DataFrame with the following columns:
- cell_ID: The identifier of the cell.
- cell_type_a: The cell type of the cell which serves as the centroid. 
- cell_type_b: The other cell type for which the proximity ratio is calculated.
- b_count: The count of cell_type_b cells within the proximity of cell_type_a.
- all_count: The count of all cells of any cell type within the proximity of cell_type_a.
- b_ratio: The ratio of b_count to all_count, representing the proximity ratio.
- d_min: The actual minimum distance used for proximity calculations.
- d_max: The actual maximum distance used for proximity calculations.
- condition: The unique condition from the `spatial_data`.
- sample_ID: The unique identifier for the sample from `spatial_data`.

Usage Example:
--------------
results_df = get_spatial_stats(spatial_data=df, cell_type_column='cell_type', d_min_scale=1, d_max_scale=5)
"""

def get_spatial_stats(spatial_data, proportion_data, cell_type_column, d_min_scale=0, d_max_scale=5, verbose=False):
    tree = KDTree(spatial_data[['x', 'y']])
    d_scale = np.median(tree.query(spatial_data[['x', 'y']], k=2)[0][:, 1])

    d_min = d_min_scale * d_scale
    d_max = d_max_scale * d_scale
    if verbose:
        print(f"d_scale {orig_ident}: {d_scale}")

    pairs = tree.query_pairs(d_max)
    if d_min > 0:
        pairs -= tree.query_pairs(d_min)
    mat = np.array(list(pairs))
    mat = np.concatenate((mat, mat[:, ::-1]))
    sparse_mat = coo_array((np.ones(len(mat), dtype=bool), mat.T),
                           shape=(len(spatial_data), len(spatial_data))).tocsr()

    condition = spatial_data['condition'].unique()[0]
    orig_ident = spatial_data['orig.ident'].unique()[0]

    results = []
    for cell_type_b in spatial_data[cell_type_column].unique():
        cell_type_mask = spatial_data[cell_type_column].values == cell_type_b
        cell_b_count = sparse_mat[:, cell_type_mask].sum(axis=1)
        all_count = sparse_mat.sum(axis=1)
        with np.errstate(invalid='ignore'):
            cell_b_ratio = cell_b_count / all_count - proportion_data.loc[proportion_data['orig.ident'] == orig_ident, f"proportion_{cell_type_b}"].values[0]
        results.append(pd.DataFrame({
            'cell_ID': spatial_data['cell_ID'].values,
            'cell_type_a': spatial_data[cell_type_column].values,
            'cell_type_b': cell_type_b,
            'b_count': cell_b_count,
            'all_count': all_count,
            'b_ratio': cell_b_ratio,
            'd_min': d_min,
            'd_max': d_max,
            'condition': condition,
            'sample_ID': orig_ident
        }))
    results = pd.concat(results).reset_index()
    gc.collect()
    return results


"""
Calculate ligand-receptor (LR) interaction statistics within a spatial dataset.

This function analyzes spatial co-localization of ligand-receptor pairs in cellular data, 
considering spatial proximity and expression levels. It uses a KDTree for efficient spatial querying 
and Pearson's correlation to assess the relationship between ligand and receptor expressions.

Parameters:
- spatial_data (pd.DataFrame): A DataFrame containing spatial coordinates (x, y), cell identifiers, 
  and other metadata. Must include 'Condition' and 'orig.ident' columns.
- expr_data (pd.DataFrame): A DataFrame containing expression data for the ligands and receptors.
- LR_pairs (pd.DataFrame): A DataFrame specifying pairs of ligands and receptors for analysis.
- cell_type_column (str): The column name in `spatial_data` that contains the cell type information.
- d_max_scale (float, optional): The maximum distance scale factor for proximity calculations, defaulting to 5.
- min_pairs (int, optional): The minimum number of cell pairs required for analysis, defaulting to 10.
- min_dup_pct (float, optional): The minimum percentage of expression data that is allowed to be identical, defaulting to 0.90.
- plotting (bool, optional): If True, generates scatter plots of ligand-receptor expression correlations. Not recommended, 
  as this generates 10s of thousands of plots and slows down the function significantly.
- save_plots (str, optional): The path to save plots, if plotting is enabled.

Returns:
- pd.DataFrame: A DataFrame with columns detailing cell types, ligands, receptors, correlation statistics, 
  and other relevant information. Specifically includes 'cell_type_a', 'cell_type_b', 'ligand', 'receptor', 
  'n' (number of pairs), 'cor' (Pearson correlation), 'pvalue' (significance of correlation), 'd_max' 
  (maximum distance for spatial proximity), 'condition', and 'sample_ID'.

Usage Example:
--------------
results_df = get_LR_stats(spatial_data=df_spatial, expr_data=df_expr, 
                          LR_pairs=df_LR, cell_type_column='cell_type', d_max_scale=5)
"""

def get_LR_stats(spatial_data, expr_data, LR_pairs, cell_type_column, d_max_scale=5, min_pairs=10, min_dup_pct=0.90, plotting=False, save_plots=None):

    assert all(spatial_data['cell_ID'] == expr_data.index)
    coords = spatial_data[['x', 'y']]
    d_scale = np.median(KDTree(coords).query(coords, k=2)[0][:, 1])
    d_max = d_max_scale * d_scale

    condition = spatial_data['Condition'].unique()[0]
    orig_ident = spatial_data['orig.ident'].unique()[0]
    unique_cell_types = spatial_data[cell_type_column].unique()
    valid_lr_pairs = set(
        (ligand, receptor)
        for lr in LR_pairs.itertuples()
        for ligand, receptor in [(lr.ligand, lr.receptor), (lr.receptor, lr.ligand)]
        if {ligand, receptor}.issubset(expr_data.columns)
    )
    results = []
    for cell_type_a in unique_cell_types:
        mask_a = spatial_data[cell_type_column] == cell_type_a
        coords_a = coords[mask_a]
        indices_a = np.where(mask_a)[0]
        tree_a = KDTree(coords_a)
        for cell_type_b in unique_cell_types:
            mask_b = spatial_data[cell_type_column] == cell_type_b
            coords_b = coords[mask_b]
            indices_b = np.where(mask_b)[0]
            if cell_type_a == cell_type_b:
                nearest_indices = tree_a.query(coords_b, k=2, distance_upper_bound=d_max)[1][:, 1]
            else:
                nearest_indices = tree_a.query(coords_b, k=1, distance_upper_bound=d_max)[1]
            valid_pairs = nearest_indices != len(coords_a)
            cell_a_indices = indices_a[nearest_indices[valid_pairs]]
            cell_b_indices = indices_b[valid_pairs]
            if len(cell_a_indices) < min_pairs:
                print(f"Insufficient paired cells for cell types {cell_type_a} and {cell_type_b}.")
                continue
            for ligand, receptor in valid_lr_pairs:
                ligand_expr = expr_data.loc[spatial_data.iloc[cell_a_indices]['cell_ID'], ligand].reset_index(drop=True)
                receptor_expr = expr_data.loc[spatial_data.iloc[cell_b_indices]['cell_ID'], receptor].reset_index(drop=True)
                if any((ligand_expr.value_counts(normalize=True) > min_dup_pct).head(1)) or \
                        any((receptor_expr.value_counts(normalize=True) > min_dup_pct).head(1)):
                    print(f"Insufficient variability for ligand {ligand} and receptor {receptor} in {cell_type_a}-{cell_type_b} in {orig_ident}")
                    continue
                cor, pvalue = pearsonr(ligand_expr, receptor_expr)
                if plotting and pvalue < 0.05:
                    plt.figure(figsize=(6, 6))
                    sns.scatterplot(x=ligand_expr, y=receptor_expr, color='blue')
                    plt.title(f"{orig_ident}")
                    plt.xlabel(f"{ligand}   {cell_type_a}")
                    plt.ylabel(f"{receptor}   {cell_type_b}")
                    plt.axline(xy1=(0, 0), slope=cor, color='red', linestyle='--')
                    plt.text(0.05, 0.95, f'Pearson R: {cor:.2f}\nP-value: {pvalue:.2e}',
                             transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
                    if save_plots:
                        plot_filename = f"{save_plots}/{cell_type_a}_{cell_type_b}_{ligand}_{receptor}_{orig_ident}.png"
                        plt.savefig(plot_filename)
                    plt.show()
                results.append({'cell_type_a': cell_type_a,
                                'cell_type_b': cell_type_b,
                                'ligand': ligand,
                                'receptor': receptor,
                                'n': len(ligand_expr),
                                'cor': cor,
                                'pvalue': pvalue,
                                'd_max': d_max,
                                'condition': condition,
                                'sample_ID': orig_ident})
    gc.collect()
    return pd.DataFrame(results).reset_index(drop=True)
