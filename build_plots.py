import pandas as pd
from pandas.api.types import CategoricalDtype
import matplotlib.pyplot as plt
import csv
import os
from matplotlib.backends.backend_pdf import PdfPages
import ast
import sys
sys.path.insert(1, '/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Inkar/scripts/bmark')

import bmark
from bmark.utils.config import load_config
from bmark.utils.analysis import order_cat,get_scores,plot_confidence, plot_confusion, plot_scores, plot_confusion

import seaborn as sns
sns.set()

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def generate_plots(figures_folder_path, 
                   df_mapping_predictions, mapping_prediction_scores, 
                   cluster_name, mapping_method, 
                   cluster_count_order, cluster_id_order, save_pdf=False):
    # plot scores and confidence
    df_mapping_predictions['true'] = df_mapping_predictions['true'].astype(cluster_count_order)
    fig_size=(30, 6) if cluster_name.lower() != 'cluster' else (70, 10) #50,20
    _ = plot_scores(df=mapping_prediction_scores, figsize=fig_size)
    _ = plot_confidence(df=df_mapping_predictions, figsize=fig_size)
    #
    # plot confusion
    df_mapping_predictions['true'] = df_mapping_predictions['true'].astype(cluster_id_order)
    fig_size=(20, 20) if cluster_name.lower() != 'cluster' else (80,80)
    _ = plot_confusion(df=df_mapping_predictions, figsize=(fig_size))
    #
    #save plots
    pdf_path = os.path.join(figures_folder_path, f"{mapping_method.lower()}_{cluster_name}_figures.pdf")
    plot_names = ["F1 scores", "Model Confidence", "Heatmap"]
    
    print("Saving plots to PNG and/or PDF")
    with PdfPages(pdf_path) as pdf:
        for i in plt.get_fignums():
            fig = plt.figure(i)
            ax = plt.gca()
            ax.set_xlabel(cluster_name.capitalize())
            fig_name = os.path.join(figures_folder_path, f'{mapping_method.lower()}_{cluster_name}_figure_{i}.png')
            plt.savefig(fig_name, bbox_inches='tight')
            if save_pdf:
                fig_title = f"{mapping_method.capitalize()} {cluster_name.capitalize()} {plot_names[i-1]}"
                fig.suptitle(fig_title, fontsize=16)
                # Save to PDF
                pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

def write_metrics_to_file(metrics_file_path, mapping_prediction_scores, mapping_method, cluster_name):
    with open(metrics_file_path, 'a', newline='') as file:
        writer = csv.writer(file)
        avg_f1 = mapping_prediction_scores["f1"].sum()/len(mapping_prediction_scores["f1"])
        avg_precision = mapping_prediction_scores["precision"].sum()/len(mapping_prediction_scores["precision"])
        avg_recall = mapping_prediction_scores["recall"].sum()/len(mapping_prediction_scores["recall"])
        
        print(f"writing to file {mapping_method.lower()}_{cluster_name}, f1 {avg_f1}, rec {avg_recall}, prec {avg_precision}. \n")
        writer.writerow([f"{mapping_method.lower()}_{cluster_name}", avg_f1, avg_recall, avg_precision])

def get_mapping_adata_pred_true_score_columns(mapping_results_csv, hierarchy_name, mapping_method, ref_dataset_name):
    """
    Selects only the three relevant columns for a hierarchy:
    assignment, bootstrapping_probability, avg_correlation
    """
    # Build the expected column names
    if ref_dataset_name is not None:
        assignment_col = f"{hierarchy_name}_{mapping_method}_{ref_dataset_name}"
        boot_col = f"bootstrapping_probability.{mapping_method}.{hierarchy_name}_{ref_dataset_name}"
        corr_col = f"avg_correlation.{mapping_method}.{hierarchy_name}_{ref_dataset_name}"
    else:
        # fallback if ref_dataset_name not provided
        assignment_col = f"{hierarchy_name}.assignment"
        boot_col = f"{hierarchy_name}.bootstrapping_probability"
        corr_col = f"{hierarchy_name}.avg_correlation"

    # Validate columns exist
    for col in [assignment_col, boot_col, corr_col]:
        if col not in mapping_results_csv.columns:
            raise KeyError(f"Column {col} not found in CSV")

    df_mapping_predictions = mapping_results_csv[[hierarchy_name, assignment_col, boot_col, corr_col]].copy()
    
    #
    # get bootstrapping_probability and avg_correlation column names for a given cluster
    numeric_columns = df_mapping_predictions.select_dtypes(include=['number']).columns.to_list()
    df_mapping_predictions.loc[:, 'score'] = (df_mapping_predictions[numeric_columns[0]] + df_mapping_predictions[numeric_columns[1]])/2
    df_mapping_predictions = df_mapping_predictions.drop(columns=numeric_columns)
    #
    df_mapping_predictions.columns=['true','pred','conf']
    return df_mapping_predictions

def get_cluster_and_algorithm_names(mapping_dataframe_results):
    #remove score columns, leave only label columns
    non_numeric_columns = mapping_dataframe_results.select_dtypes(exclude=['number']).columns.to_list()
    #
    # generalized lambda function for pick fist(cluster names) and last(mapping algorithms) from a list of strings
    # this could be a simple for loop
    extract_first_last = lambda input_string: (input_string.split('_')[0], input_string.split('_')[-1])
    #
    # Apply the lambda function to each string in the list
    first_words, last_words = zip(*map(extract_first_last, non_numeric_columns))
    #
    # Convert to lists with unique values, fist words are cluster names and last are mapping algorithms 
    first_words = list(first_words) # cluster names
    last_words = list(last_words) # mapping algorithms
    #
    return non_numeric_columns, first_words, last_words

# create folder with the taxonomy name and csv, plots folders inside it
def get_folder_paths(mapping_results_csv_file_path):
    # csv folder
    metrics_folder_path = os.path.join(os.path.dirname(mapping_results_csv_file_path), "csv")
    os.makedirs(metrics_folder_path, exist_ok=True)
    metrics_csv_file_path = os.path.join(metrics_folder_path, "performance_metrics.csv")  
    #
    # plots folder
    figures_folder_path = os.path.join(os.path.dirname(mapping_results_csv_file_path),"plots")
    os.makedirs(figures_folder_path, exist_ok=True)
    #
    return metrics_csv_file_path, figures_folder_path

def build_benchmarks(working_dir, hierarchy, ref_dataset_name, mapping_method, save_pdf):
    print("=========== Building Plots ===========")
    #
    mapping_results_csv_file_path = os.path.join(working_dir, f"results/{mapping_method}/{mapping_method}_predictions.csv")
    mapping_results_csv = pd.read_csv(mapping_results_csv_file_path)
    #
    # create folders for saving the plots and csv files
    metrics_csv_file_path, figures_folder_path = get_folder_paths(mapping_results_csv_file_path)
    #
    # Open the file for writing
    with open(metrics_csv_file_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Model Name", "F1 Score", "Recall", "Precision"])
    #
    for i in range(len(hierarchy)):
        print(f"Processing hierarchy level {hierarchy[i]} and method {mapping_method}.")
        # order of categorical labels for plotting convenience
        if hierarchy[i] == "Group":
            anno_dir = "./group_order.csv"
            metadata = pd.read_csv(anno_dir)
            col_order = metadata['Group']
            mapping_results_csv['Group'] = pd.Categorical(mapping_results_csv['Group'], categories=col_order, ordered=True)
        df = mapping_results_csv.rename(columns={f"{hierarchy[i]}": f"{hierarchy[i]}_label"}, inplace=False)
        df[f"{hierarchy[i]}_id"] = df[f"{hierarchy[i]}_label"]
        _, cluster_count_order = order_cat(df, cat=hierarchy[i], by='id') # count
        _, cluster_id_order = order_cat(df, cat=hierarchy[i], by='id')
        #
        # make a sub dataframe with true and pred columns
        df_mapping_predictions = get_mapping_adata_pred_true_score_columns(mapping_results_csv, hierarchy[i], mapping_method, ref_dataset_name)
        #
        # plot scores and confidence
        df_mapping_predictions['true'] = df_mapping_predictions['true'].astype(cluster_count_order)
        mapping_prediction_scores = get_scores(df=df_mapping_predictions)
        #
        write_metrics_to_file(metrics_csv_file_path, mapping_prediction_scores, mapping_method, hierarchy[i])
        generate_plots(figures_folder_path, 
                df_mapping_predictions, mapping_prediction_scores, 
                hierarchy[i], mapping_method, 
                cluster_count_order, cluster_id_order, save_pdf)
