.get_relevant_args_list <- function(
# .get_relevant_args_list <- function(out_dir=NULL,
#                                     HMM=FALSE,
#                                     HMM_type='i6',
#                                     num_ref_groups=NULL,
#                                     analysis_mode='samples',
#                                     tumor_subcluster_partition_method,
#                                     smooth_method='pyramidinal',
#                                     max_centered_threshold=3,
#                                     remove_genes_at_chr_ends=FALSE,
#                                     prune_outliers=FALSE,
#                                     BayesMaxPNormal=0.5,
#                                     mask_nonDE_genes=FALSE,
#                                     denoise=FALSE,
#                                     noise_filter=NA,
#                                     sd_amplifier=1.5,
#                                     noise_logistic=FALSE,
                                    ...) {
    run_arguments = list(...)

    # creation args  ## add check for matrix size, col/row names, and ref/obs indices
    relevant_args = vector("list", 21) # 21 = max steps
    expected_file_names = vector("character", 21)
    step_i = 1

    # 1 _incoming_data.infercnv_obj
    relevant_args[[step_i]] = c("chr_exclude", "max_cells_per_group", "min_max_counts_per_cell", "counts_md5")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_incoming_data.infercnv_obj", step_i))
    step_i = step_i + 1

    # 2 _reduced_by_cutoff.infercnv_obj
    relevant_args[[step_i]] = c("cutoff", "min_cells_per_gene")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_reduced_by_cutoff.infercnv_obj", step_i))
    step_i = step_i + 1

    # 3 _normalized_by_depth%s.infercnv_obj
    relevant_args[[step_i]] = c("HMM")
    if (run_arguments$HMM) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "HMM_type") 
        if (run_arguments$HMM_type == "i6") {
            relevant_args[[step_i]] = c(relevant_args[[step_i]], "sim_method", "hspike_aggregate_normals", "sim_foreground")
        }
    }
    resume_file_token = ifelse((run_arguments$HMM), paste0("HMM", run_arguments$HMM_type), "")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_normalized_by_depth%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 4 _logtransformed%s.infercnv_obj
    relevant_args[[step_i]] = c()
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_logtransformed%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 5 _scaled%s.infercnv_obj
    relevant_args[[step_i]] = c("scale_data")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_scaled%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 6 _split_%sf_refs%s.infercnv_obj
    relevant_args[[step_i]] = c("num_ref_groups")
    if (!is.null(run_arguments$num_ref_groups)) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "hclust_method")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_split_%sf_refs%s.infercnv_obj", step_i, resume_file_token, run_arguments$num_ref_groups))
    }
    step_i = step_i + 1

    # 7 _tumor_subclusters%s.%s.infercnv_obj
    relevant_args[[step_i]] = c("analysis_mode", "tumor_subcluster_partition_method")
    if (run_arguments$analysis_mode == 'subclusters' & run_arguments$tumor_subcluster_partition_method == 'random_trees') {
        resume_file_token = paste0(resume_file_token, ".rand_trees")
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "hclust_method", "tumor_subcluster_pval", "cluster_by_groups")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_tumor_subclusters%s.%s.infercnv_obj", step_i, resume_file_token, run_arguments$tumor_subcluster_partition_method))
    }
    step_i = step_i + 1

    # 8 _remove_ref_avg_from_obs_logFC%s.infercnv_obj
    relevant_args[[step_i]] = c("ref_subtract_use_mean_bounds")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_remove_ref_avg_from_obs_logFC%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1




    # 9 _apply_max_centered_expr_threshold%s.infercnv_obj
    relevant_args[[step_i]] = c("max_centered_threshold")
    if(!is.na(run_arguments$max_centered_threshold)) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "max_centered_threshold")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_apply_max_centered_expr_threshold%s.infercnv_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 10 _smoothed_by_chr%s.infercnv_obj
    relevant_args[[step_i]] = c("smooth_method", "window_length")
    if (run_arguments$smooth_method == 'pyramidinal') {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "smooth_ends")
    }
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_smoothed_by_chr%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 11 _recentered_cells_by_chr%s.infercnv_obj
    relevant_args[[step_i]] = c()
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_recentered_cells_by_chr%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 12 _remove_ref_avg_from_obs_adjust%s.infercnv_obj
    relevant_args[[step_i]] = c("ref_subtract_use_mean_bounds")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_remove_ref_avg_from_obs_adjust%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 13 _remove_gene_at_chr_ends%s.infercnv_obj
    relevant_args[[step_i]] = c("remove_genes_at_chr_ends")
    if (run_arguments$remove_genes_at_chr_ends == TRUE) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "window_length")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_remove_gene_at_chr_ends%s.infercnv_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 14 _invert_log_transform%s.infercnv_obj
    relevant_args[[step_i]] = c()
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_invert_log_transform%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 15 _tumor_subclusters%s.infercnv_obj
    relevant_args[[step_i]] = c("analysis_mode", "tumor_subcluster_partition_method", "z_score_filter")
    if (run_arguments$analysis_mode == 'subclusters' & run_arguments$tumor_subcluster_partition_method != 'random_trees') {
        resume_file_token = paste0(resume_file_token, '.', run_arguments$tumor_subcluster_partition_method)
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "k_nn", "leiden_resolution", "tumor_subcluster_pval", "leiden_method", "leiden_function", "leiden_method_per_chr", "leiden_function_per_chr", "leiden_resolution_per_chr", "per_chr_hmm_subclusters", "per_chr_hmm_subclusters_references", "hclust_method", "cluster_by_groups") 
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_tumor_subclusters%s.infercnv_obj", step_i, resume_file_token))
    }
    else if (run_arguments$analysis_mode != 'subclusters') {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "tumor_subcluster_pval", "hclust_method", "cluster_by_groups")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_no_subclustering%s.infercnv_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 16 _remove_outlier%s.infercnv_obj
    relevant_args[[step_i]] = c("prune_outliers")
    if (run_arguments$prune_outliers) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "outlier_method_bound", "outlier_lower_bound", "outlier_upper_bound")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_remove_outlier%s.infercnv_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 17 _HMM_pred%s.infercnv_obj
    relevant_args[[step_i]] = c("HMM")
    if (run_arguments$HMM) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "HMM_type", "analysis_mode", "HMM_transition_prob", "HMM_report_by")
        hmm_resume_file_token = paste0(resume_file_token, ".hmm_mode-", run_arguments$analysis_mode)
        # hmm.infercnv_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred%s.infercnv_obj", step_i, hmm_resume_file_token))  ########
        if (run_arguments$analysis_mode == 'subclusters' && run_arguments$tumor_subcluster_partition_method == 'random_trees') {
            relevant_args[[step_i]] = c(relevant_args[[step_i]], "hclust_method")
        }

        if (run_arguments$HMM_type == "i6") {
            # relevant_args[[step_i]] = c(relevant_args[[step_i]], )

        }
        else {
            relevant_args[[step_i]] = c(relevant_args[[step_i]], "HMM_i3_pval", "HMM_i3_use_KS")
        }
        expected_file_names[[step_i]] = hmm.infercnv_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred%s.infercnv_obj", step_i, hmm_resume_file_token))
    }


    step_i = step_i + 1

    # 18 _HMM_pred.Bayes_Net%s.mcmc_obj
    relevant_args[[step_i]] = c("HMM")
    if (run_arguments$HMM == TRUE & run_arguments$BayesMaxPNormal > 0) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "diagnostics", "reassignCNVs")
        # mcmc_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token))
        # mcmc.infercnv_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.Pnorm_%g.infercnv_obj", step_i, hmm_resume_file_token, BayesMaxPNormal))
        # expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token)), 
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token))
    }
    step_i = step_i + 1

    # 19 _HMM_pred.Bayes_Net%s.Pnorm_%g.infercnv_obj
    relevant_args[[step_i]] = c("HMM", "BayesMaxPNormal")
    if (run_arguments$HMM == TRUE & run_arguments$BayesMaxPNormal > 0) {
        # relevant_args[[step_i]] = c(relevant_args[[step_i]], "diagnostics", "reassignCNVs")
        # mcmc_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token))
        # mcmc.infercnv_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.Pnorm_%g.infercnv_obj", step_i, hmm_resume_file_token, BayesMaxPNormal))
        # expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token)), 
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.Pnorm_%g.infercnv_obj", step_i, hmm_resume_file_token, run_arguments$BayesMaxPNormal))
    }
    step_i = step_i + 1

    # 20 _HMM_pred.repr_intensities%s.Pnorm_%g.infercnv_obj
    relevant_args[[step_i]] = c("HMM")
    if (run_arguments$HMM) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "BayesMaxPNormal")
        # hmm.infercnv_obj_file = file.path(out_dir, sprintf("%02d_HMM_pred.repr_intensities%s.Pnorm_%g.infercnv_obj", step_i, hmm_resume_file_token, BayesMaxPNormal))
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.repr_intensities%s.Pnorm_%g.infercnv_obj", step_i, hmm_resume_file_token, run_arguments$BayesMaxPNormal))
    }
    step_i = step_i + 1

    # 21 _mask_nonDE%s.infercnv_obj
    relevant_args[[step_i]] = c("mask_nonDE_genes")
    if (run_arguments$mask_nonDE_genes) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "require_DE_all_normals", "test.use", "mask_nonDE_pval")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_mask_nonDE%s.infercnv_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 22 _denoise%s.NF_%s.SD_%g.NL_%s.infercnv_obj
    relevant_args[[step_i]] = c("denoise")
    if (run_arguments$denoise) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "noise_filter", "sd_amplifier", "noise_logistic")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_denoise%s.NF_%s.SD_%g.NL_%s.infercnv_obj", step_i, resume_file_token, run_arguments$noise_filter, run_arguments$sd_amplifier, run_arguments$noise_logistic))
    }
    step_i = step_i + 1

    return(list(relevant_args=relevant_args, expected_file_names=expected_file_names))
}


