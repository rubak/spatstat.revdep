# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

compute_split_frequencies <- function(forest_object, max_depth) {
    .Call('_drf_compute_split_frequencies', PACKAGE = 'drf', forest_object, max_depth)
}

compute_weights <- function(forest_object, train_matrix, sparse_train_matrix, test_matrix, sparse_test_matrix, num_threads) {
    .Call('_drf_compute_weights', PACKAGE = 'drf', forest_object, train_matrix, sparse_train_matrix, test_matrix, sparse_test_matrix, num_threads)
}

compute_weights_oob <- function(forest_object, test_matrix, sparse_test_matrix, num_threads) {
    .Call('_drf_compute_weights_oob', PACKAGE = 'drf', forest_object, test_matrix, sparse_test_matrix, num_threads)
}

merge <- function(forest_objects) {
    .Call('_drf_merge', PACKAGE = 'drf', forest_objects)
}

gini_train <- function(train_matrix, sparse_train_matrix, outcome_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, num_features, bandwidth, node_scaling) {
    .Call('_drf_gini_train', PACKAGE = 'drf', train_matrix, sparse_train_matrix, outcome_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, num_features, bandwidth, node_scaling)
}

fourier_train <- function(train_matrix, sparse_train_matrix, outcome_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, num_features, bandwidth, node_scaling) {
    .Call('_drf_fourier_train', PACKAGE = 'drf', train_matrix, sparse_train_matrix, outcome_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, num_features, bandwidth, node_scaling)
}

regression_predict <- function(forest_object, train_matrix, sparse_train_matrix, outcome_index, test_matrix, sparse_test_matrix, num_threads, estimate_variance) {
    .Call('_drf_regression_predict', PACKAGE = 'drf', forest_object, train_matrix, sparse_train_matrix, outcome_index, test_matrix, sparse_test_matrix, num_threads, estimate_variance)
}

regression_predict_oob <- function(forest_object, train_matrix, sparse_train_matrix, outcome_index, num_threads, estimate_variance) {
    .Call('_drf_regression_predict_oob', PACKAGE = 'drf', forest_object, train_matrix, sparse_train_matrix, outcome_index, num_threads, estimate_variance)
}

