
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

import pandas as pd
import numpy as np
from scipy.spatial import KDTree
from scipy.sparse import coo_array

def get_spatial_stats(spatial_data, cell_type_column, d_min_scale=0, d_max_scale=20):
    tree = KDTree(spatial_data[['x', 'y']])
    d_scale = np.median(tree.query(spatial_data[['x', 'y']], k=2)[0][:, 1])
    
    d_min = d_min_scale * d_scale
    d_max = d_max_scale * d_scale
    
    condition = spatial_data['Condition'].unique()[0]
    orig_ident = spatial_data['orig.ident'].unique()[0]
    print(f"d_scale {orig_ident}: {d_scale}")
    
    pairs = tree.query_pairs(d_max)
    if d_min > 0:
        pairs -= tree.query_pairs(d_min)
    mat = np.array(list(pairs))
    mat = np.concatenate((mat, mat[:, ::-1]))
    sparse_mat = coo_array((np.ones(len(mat), dtype=bool), mat.T),
        shape=(len(spatial_data), len(spatial_data))).tocsr()

    results = []
    for cell_type_b in spatial_data[cell_type_column].unique():
        cell_type_mask = spatial_data[cell_type_column].values == cell_type_b
        cell_b_count = sparse_mat[:, cell_type_mask].sum(axis=1)
        all_count = sparse_mat.sum(axis=1)
        with np.errstate(invalid='ignore'):
            cell_b_ratio = cell_b_count / all_count
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
    return results
