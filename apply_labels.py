import json
import os
import pandas as pd
from pandas import json_normalize


def read_json_file(file_path):     
    with open(file_path, 'r') as json_file:         
        data = json.load(json_file)     
        return data

def apply_taxonomy_tree_label_names(json_results, mapping_results_df, json_hierarchy):
    taxonomy_tree = json_normalize(json_results["taxonomy_tree"]["name_mapper"])
    taxonomy_tree = taxonomy_tree[[col for col in taxonomy_tree.columns if 'name' in col]]  # keep only name columns, ex CCN20250428_LEVEL_4.CS20250428_CLUST_0161.name
    taxonomy_tree.columns = taxonomy_tree.columns.str.split('.').str[1]    # rename the columns, ex df{CCN20250428_LEVEL_4.CS20250428_CLUST_0161.name: Human-143} => df{CS20250428_CLUST_0161: Human-143}. CCN20250428_LEVEL_4 (already used in the lookup key below CCN20250428_LEVEL_4.assignment) is the assignment variable name in json object and CS20250428_CLUST_0161 is the actual assigned label.
    taxonomy_tree = taxonomy_tree.to_dict('records')[0]
    #
    mapping_results_df[f"{json_hierarchy[0]}.assignment"] = mapping_results_df[f"{json_hierarchy[0]}.assignment"].map(taxonomy_tree)  # ex) for CCN20250428_LEVEL_4 swap all 'CS20250428_CLUST_0161' to 'Human-143'
    mapping_results_df[f"{json_hierarchy[1]}.assignment"] = mapping_results_df[f"{json_hierarchy[1]}.assignment"].map(taxonomy_tree)
    mapping_results_df[f"{json_hierarchy[2]}.assignment"] = mapping_results_df[f"{json_hierarchy[2]}.assignment"].map(taxonomy_tree)
    mapping_results_df[f"{json_hierarchy[3]}.assignment"] = mapping_results_df[f"{json_hierarchy[3]}.assignment"].map(taxonomy_tree)
    mapping_results_df[f"{json_hierarchy[4]}.assignment"] = mapping_results_df[f"{json_hierarchy[4]}.assignment"].map(taxonomy_tree)
    #
    return mapping_results_df

def combine_predicted_labels(working_dir, adata, hierarchy, json_hierarchy, n_folds, ref_dataset_name, mapping_method="Hierarchical", apply_tax_tree_label_names = False):
    os.makedirs(working_dir, exist_ok=True)
    output_dir = os.path.join(working_dir, f"results/{mapping_method}")
    os.makedirs(output_dir, exist_ok=True)
    folds_dir = os.path.join(working_dir, 'folds')
    # Initialize storage for each hierarchy level
    hierarchy_dfs = {h: [] for h in hierarchy}

    for fold in range(1, n_folds+1):
        # look for folds/fold_#/results_hann/flat.json file
        results_file_json = os.path.join(
            folds_dir,
            f"fold_{fold}",
            f"results_{'hann' if mapping_method=='Hierarchical' else 'flat'}.json"
        )
        print(f"Reading fold {fold}: {results_file_json}")
        # load the json file
        mapping_results = read_json_file(results_file_json)
        mapping_results_df = json_normalize(mapping_results["results"])
        mapping_results_df.set_index('cell_id', inplace=True)
        # check if labels were coded using the taxonomy tree lookup at the bottom of json file
        if apply_tax_tree_label_names:
            mapping_results_df = apply_taxonomy_tree_label_names(mapping_results, mapping_results_df, json_hierarchy) 
        # compile results df with original and predicted labels
        valid_cells = mapping_results_df.index.intersection(adata.obs.index)
        if all(mapping_results_df.index == valid_cells): # check if all mapping validation subset cells exist in adata
            print("All IDs match.")
            #true_pred_labels = adata.obs[hierarchy].copy()
            adata_subset = adata.obs.loc[valid_cells]
            for hierarchy_level, json_hierarchy_level in zip(hierarchy, json_hierarchy):
                hann_pred_labels = [f'{hierarchy_level}_{mapping_method}_{ref_dataset_name}', 
                                    f'bootstrapping_probability.{mapping_method}.{hierarchy_level}_{ref_dataset_name}', 
                                    f'avg_correlation.{mapping_method}.{hierarchy_level}_{ref_dataset_name}']
                df_fold = pd.DataFrame({
                    hierarchy_level: adata_subset[hierarchy_level].values,
                    hann_pred_labels[0]: mapping_results_df[f"{json_hierarchy_level}.assignment"],
                    hann_pred_labels[1]: mapping_results_df[f"{json_hierarchy_level}.bootstrapping_probability"],
                    hann_pred_labels[2]: mapping_results_df[f"{json_hierarchy_level}.avg_correlation"],
                    "fold": fold
                }, index=valid_cells)  # set index to cell_id
                hierarchy_dfs[hierarchy_level].append(df_fold)
    # Merge all hierarchy-level DataFrames into one combined DataFrame
    merged_df = None

    for hierarchy_level, dfs in hierarchy_dfs.items():
        for i, df in enumerate(dfs):
            if df.index.duplicated().any():
                print(f"Duplicate indexes in {hierarchy_level} fold {i+1}: {df.index[df.index.duplicated()]}")
        # Combine all folds for this hierarchy
        combined_df = pd.concat(dfs)
        
        if merged_df is None:
            merged_df = combined_df
        else:
            # Merge on 'cell_id' only, keep all columns from each hierarchy
            merged_df = merged_df.merge(combined_df.drop(columns=['fold']), left_index=True, right_index=True, how='outer') # dropping fold here since it's already at the initial df creation in the if statement above. Each cell_id (index) has only one fold it was ran in

    # Save the merged CSV
    out_path = os.path.join(output_dir, f"{mapping_method}_predictions.csv")
    merged_df.to_csv(out_path, index=True)
    print(f"Saved combined predictions for all hierarchies to {out_path}")
    print("Finished processing all folds.")        
