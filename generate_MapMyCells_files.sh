#!/bin/bash

FOLD="$1"
HIERARCHY="$2"
CURRENT_DIR="$3"

echo "FOLD: $FOLD"
echo "HIERARCHY: $HIERARCHY"
echo "WORKING_DIR: $CURRENT_DIR"

exec > MapMyCells_Output_${FOLD}.txt 2>&1

python -m cell_type_mapper.cli.precompute_stats_scrattch \
--h5ad_path "${CURRENT_DIR}/folds/fold_${FOLD}/train.h5ad" \
--hierarchy "$HIERARCHY" \
--output_path "${CURRENT_DIR}/folds/fold_${FOLD}/precompute_stats.h5" \
--normalization log2CPM

python -m cell_type_mapper.cli.reference_markers \
--precomputed_path_list "['${CURRENT_DIR}/folds/fold_${FOLD}/precompute_stats.h5']" \
--output_dir "${CURRENT_DIR}/folds/fold_${FOLD}"

python -m cell_type_mapper.cli.query_markers \
--reference_marker_path_list "['${CURRENT_DIR}/folds/fold_${FOLD}/reference_markers.h5']" \
--output_path "${CURRENT_DIR}/folds/fold_${FOLD}/query_markers.json"

python -m cell_type_mapper.cli.from_specified_markers \
--query_path "${CURRENT_DIR}/folds/fold_${FOLD}/validation.h5ad" \
--type_assignment.normalization log2CPM \
--precomputed_stats.path "${CURRENT_DIR}/folds/fold_${FOLD}/precompute_stats.h5" \
--query_markers.serialized_lookup "${CURRENT_DIR}/folds/fold_${FOLD}/query_markers.json" \
--extended_result_path "${CURRENT_DIR}/folds/fold_${FOLD}/results_hann.json" \
--type_assignment.n_processors 16 \
--type_assignment.chunk_size 1000 \
--flatten False 

python -m cell_type_mapper.cli.from_specified_markers \
--query_path "${CURRENT_DIR}/folds/fold_${FOLD}/validation.h5ad" \
--type_assignment.normalization log2CPM \
--precomputed_stats.path "${CURRENT_DIR}/folds/fold_${FOLD}/precompute_stats.h5" \
--query_markers.serialized_lookup "${CURRENT_DIR}/folds/fold_${FOLD}/query_markers.json" \
--extended_result_path "${CURRENT_DIR}/folds/fold_${FOLD}/results_flat.json" \
--type_assignment.n_processors 16 \
--type_assignment.chunk_size 1000 \
--flatten True