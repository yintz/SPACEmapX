run <- function(SPACEmapX_obj,

                # gene filtering settings
                cutoff=1,
                min_cells_per_gene=3,

                out_dir=NULL,

                ## smoothing params
                window_length=101,
                smooth_method=c('pyramidinal', 'runmeans', 'coordinates'),

                num_ref_groups=NULL,
                ref_subtract_use_mean_bounds=TRUE,

                # observation cell clustering settings
                cluster_by_groups=TRUE,
                cluster_references=TRUE,
                k_obs_groups=1,

                hclust_method='ward.D2',

                max_centered_threshold=3, # or set to a specific value or "auto", or NA to turn off
                scale_data=FALSE,

                ## HMM opts
                HMM=FALSE, # turn on to auto-run the HMM prediction of CNV levels
                ## tumor subclustering opts
                HMM_transition_prob=1e-6,
                HMM_report_by=c("subcluster","consensus","cell"),
                HMM_type=c('i6', 'i3'),
                HMM_i3_pval=0.05,
                HMM_i3_use_KS=FALSE,
                BayesMaxPNormal=0.5,

                ## some experimental params
                #sim_method=c('meanvar', 'simple', 'splatter'), ## only meanvar supported, others experimental
                sim_method='meanvar',
                sim_foreground=FALSE, ## experimental
                reassignCNVs=TRUE,


                ## tumor subclustering options
                analysis_mode=c('subclusters', 'samples', 'cells'), # for filtering and HMM
                tumor_subcluster_partition_method=c('leiden', 'random_trees', 'qnorm', 'pheight', 'qgamma', 'shc'),
                tumor_subcluster_pval=0.1,
                k_nn=20,
                leiden_method=c("PCA", "simple"),
                leiden_function = c("CPM", "modularity"),
                leiden_resolution="auto",
                leiden_method_per_chr=c("simple", "PCA"),
                leiden_function_per_chr = c("modularity", "CPM"),
                leiden_resolution_per_chr = 1,
                per_chr_hmm_subclusters=FALSE,
                per_chr_hmm_subclusters_references=FALSE,
                z_score_filter = 0.8,


                ## noise settings
                denoise=FALSE,
                noise_filter=NA,
                sd_amplifier = 1.5,
                noise_logistic=FALSE, # if false, does complete 'noise' elimination.

                # outlier adjustment settings
                outlier_method_bound="average_bound",
                outlier_lower_bound=NA,
                outlier_upper_bound=NA,

                ## misc options
                final_scale_limits = NULL,
                final_center_val = NULL,
                debug=FALSE, #for debug level logging
                num_threads = 4,
                plot_steps=FALSE,
                inspect_subclusters = TRUE,
                resume_mode=TRUE,
                png_res=300,
                plot_probabilities = TRUE,
                save_rds = TRUE,
                save_final_rds = TRUE,
                diagnostics = FALSE,

                ## experimental options
                remove_genes_at_chr_ends=FALSE,
                prune_outliers=FALSE,

                mask_nonDE_genes=FALSE,
                mask_nonDE_pval=0.05, # use permissive threshold
                test.use='wilcoxon',
                require_DE_all_normals="any",


                hspike_aggregate_normals = FALSE,

                no_plot = FALSE,
                no_prelim_plot = FALSE,
                write_expr_matrix = FALSE,
                write_phylo = FALSE,
                output_format = "png",
                plot_chr_scale = FALSE,
                chr_lengths = NULL,
                useRaster = TRUE,

                up_to_step=100

) {


    smooth_method = match.arg(smooth_method)
    HMM_type = match.arg(HMM_type)
    if (HMM && HMM_type == 'i6' && smooth_method == 'coordinates') {
        flog.error("Cannot use 'coordinate' smoothing method with i6 HMM model at this time.")
        stop("Incompatible HMM mode and smoothing method.")
    }
    if (smooth_method == 'coordinates' && window_length < 10000) {
        window_length=10000000
        flog.warn(paste0("smooth_method set to 'coordinates', but window_length ",
                         "is less than 10.000, setting it to 10.000.000. Please ",
                         "set a different value > 10.000 if this default does not seem suitable."))
    }

    leiden_function = match.arg(leiden_function)
    leiden_function_per_chr = match.arg(leiden_function_per_chr)
    HMM_report_by = match.arg(HMM_report_by)
    analysis_mode = match.arg(analysis_mode)
    if(analysis_mode == "cells" && HMM_report_by != "cell") {
        flog.warn(paste0("analysis_mode is \"cells\" but HMM_report_by is not, "),
                         "changing HMM_report_by to \"cells\".")
        HMM_report_by = "cell"
    }
    tumor_subcluster_partition_method = match.arg(tumor_subcluster_partition_method)
    leiden_method = match.arg(leiden_method)
    leiden_method_per_chr = match.arg(leiden_method_per_chr)
    if (tumor_subcluster_partition_method != "leiden") {
        per_chr_hmm_subclusters = FALSE
    }
    
    if (debug) {
        flog.threshold(DEBUG)
    } else {
        flog.threshold(INFO)
    }
    
    flog.info(paste("::process_data:Start", sep=""))
    
    SPACEmapX$GLOBAL_NUM_THREADS <- 4
    if (is.null(out_dir)) {
        flog.error("Error, out_dir is NULL, please provide a path.")
        stop("out_dir is NULL")
    }

    call_match = match.call()

    arg_names = names(call_match)
    arg_names = arg_names[-which(arg_names == "" | arg_names == "SPACEmapX_obj")]

    # evaluate variables such as output_dir
    for (n in arg_names) {
        call_match[[n]] = get(n)
    }

    # add the arguments that are checked against a list of allowed values and autofilled with the first if none is given using match.arg()
    call_match$smooth_method = smooth_method
    call_match$HMM_type = HMM_type
    call_match$HMM_report_by = HMM_report_by
    call_match$analysis_mode = analysis_mode
    call_match$tumor_subcluster_partition_method = tumor_subcluster_partition_method
    call_match$hclust_method = hclust_method
    call_match$leiden_function = leiden_function
    # add argument needed to get relevant args list
    call_match$HMM = HMM
    call_match$z_score_filter = z_score_filter
    call_match$tumor_subcluster_pval = tumor_subcluster_pval
    call_match$k_nn = k_nn
    call_match$leiden_method = leiden_method
    call_match$per_chr_hmm_subclusters = per_chr_hmm_subclusters
    call_match$cluster_by_groups = cluster_by_groups
    call_match$max_centered_threshold = max_centered_threshold
    call_match$remove_genes_at_chr_ends = remove_genes_at_chr_ends
    call_match$prune_outliers = prune_outliers
    call_match$BayesMaxPNormal = BayesMaxPNormal
    call_match$mask_nonDE_genes = mask_nonDE_genes
    call_match$denoise = denoise
    call_match$noise_filter = noise_filter
    call_match$sd_amplifier = sd_amplifier
    call_match$noise_logistic = noise_logistic

    call_match$SPACEmapX_obj = NULL

    # pull the list of argument that are set by the user in case they differ from the default
    current_args = as.list(call_match[-which(names(call_match) == "" | names(call_match) == "SPACEmapX_obj")])
    # add the arguments that are checked against a list of allowed values and autofilled with the first if none is given using match.arg()
    # current_args$smooth_method = smooth_method
    # current_args$HMM_type = HMM_type
    # current_args$HMM_report_by = HMM_report_by
    # current_args$analysis_mode = analysis_mode
    # current_args$tumor_subcluster_partition_method = tumor_subcluster_partition_method

    if(out_dir != "." && !file.exists(out_dir)){
        flog.info(paste0("Creating output path ", out_dir))
        dir.create(out_dir)
    }

    SPACEmapX_obj@options = c(current_args, SPACEmapX_obj@options)

    #run_call <- match.call()
    call_match[[1]] <- as.symbol(".get_relevant_args_list")
    reload_info = eval(call_match)

    # reload_info$relevant_args
    # reload_info$expected_file_names

    skip_hmm = 0
    # skip_mcmc = 0
    skip_past = 0
    
    if (resume_mode) {
        flog.info("Checking for saved results.")
        for (i in rev(seq_along(reload_info$expected_file_names))) {
            if (file.exists(reload_info$expected_file_names[[i]])) {
                flog.info(paste0("Trying to reload from step ", i))
                # if ((i == 17 && skip_hmm == 0) || (i %in% (18:20))) {
                if ((i == 20) || (i %in% seq(17, 19) && skip_hmm == 0) || (i == 17 && skip_hmm == 2)) {   # step 17 appears in two conditions because if last found step is 18, step 19 requires the results of both step 17 and 18 to run
                    if (i == 18 && BayesMaxPNormal > 0) {  # mcmc_obj
                        mcmc_obj = readRDS(reload_info$expected_file_names[[i]])
                        if (!.compare_args(SPACEmapX_obj@options, unlist(reload_info$relevant_args[1:i]), mcmc_obj@options)) {
                            rm(mcmc_obj)
                            invisible(gc())
                        }
                        else {
                            mcmc_obj@options = SPACEmapX_obj@options
                            skip_hmm = max(skip_hmm, i - 16)
                            flog.info(paste0("Using backup MCMC from step ", i))
                        }
                    }
                    else {
                        hmm.SPACEmapX_obj = readRDS(reload_info$expected_file_names[[i]])
                        if (!.compare_args(SPACEmapX_obj@options, unlist(reload_info$relevant_args[1:i]), hmm.SPACEmapX_obj@options)) {
                            rm(hmm.SPACEmapX_obj)
                            invisible(gc())
                        }
                        else {
                            hmm.SPACEmapX_obj@options = SPACEmapX_obj@options
                            skip_hmm = max(skip_hmm, i - 16)  # max for case where step 18 is found, so can be skipped, but still needed to reload step 17 results to be able to run step 19
                            flog.info(paste0("Using backup HMM from step ", i))
                        }
                    }
                }
                else if (!(i %in% 17:20)) {
                    reloaded_SPACEmapX_obj = readRDS(reload_info$expected_file_names[[i]])
                    if (skip_past > i) { # in case denoise was found
                        if (21 > i) { # if 22/21 already found and checked HMM too, stop
                            break
                        }
                        else { # if 22 denoise found, don't check 21 maskDE
                            next
                        }
                    }
                    if (.compare_args(SPACEmapX_obj@options, unlist(reload_info$relevant_args[1:i]), reloaded_SPACEmapX_obj@options)) {
                        if (i ==  15 && per_chr_hmm_subclusters) {
                            tmp = nchar(reload_info$expected_file_names[[i]])
                            tmp = paste0(substr(reload_info$expected_file_names[[i]], 1, tmp-13), ".per_chr_subclusters", substr(reload_info$expected_file_names[[i]], tmp-12, tmp))
                            if (file.exists(tmp)) {
                                subclusters_per_chr = readRDS(tmp)
                            }
                            else {
                                next
                            }
                        }
                        options_backup = SPACEmapX_obj@options
                        SPACEmapX_obj = reloaded_SPACEmapX_obj # replace input SPACEmapX_obj
                        rm(reloaded_SPACEmapX_obj) # remove first (temporary) reference so there's no duplication when they would diverge
                        SPACEmapX_obj@options = options_backup
                        skip_past = i
                        flog.info(paste0("Using backup from step ", i))
                        if (i < 21) { # do not stop check right away if denoise was found to also allow to check for HMM
                            break
                        }
                    }
                    else {
                        rm(reloaded_SPACEmapX_obj)
                        invisible(gc())
                    }
                }
            }
        }
    }
 

    step_count = 0;
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 1
    flog.info(sprintf("\n\n\tSTEP %d: incoming data\n", step_count))

    if (skip_past < step_count && save_rds) {
        saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
    }
    
    ## #################################################
    ## Step: removing insufficiently expressed genes
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 2
    flog.info(sprintf("\n\n\tSTEP %02d: Removing lowly expressed genes\n", step_count))
    
    ## Remove genes that aren't sufficiently expressed, according to min mean count cutoff.
    ## Examines the original (non-log-transformed) data, gets mean for each gene, and removes genes
    ##  with mean values below cutoff.

    if (skip_past < step_count) {
    # }
        SPACEmapX_obj <- require_above_min_mean_expr_cutoff(SPACEmapX_obj, cutoff)
        
        ## require each gene to be present in a min number of cells for ref sets
        
        SPACEmapX_obj <- require_above_min_cells_ref(SPACEmapX_obj, min_cells_per_gene=min_cells_per_gene)
        
        if (save_rds) {
            saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
        }
    }
    
    
    ## #########################################
    ## # STEP: normalization by sequencing depth

    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 3
    flog.info(sprintf("\n\n\tSTEP %02d: normalization by sequencing depth\n", step_count))
    
    
    resume_file_token = ifelse( (HMM), paste0("HMM",HMM_type), "")

    if (skip_past < step_count) {
        SPACEmapX_obj <- normalize_counts_by_seq_depth(SPACEmapX_obj)
        
        if (HMM && HMM_type == 'i6') {
            ## add in the hidden spike needed by the HMM
            SPACEmapX_obj <- .build_and_add_hspike(SPACEmapX_obj, sim_method=sim_method, aggregate_normals=hspike_aggregate_normals)
            
            if (sim_foreground) {
                SPACEmapX_obj <- .sim_foreground(SPACEmapX_obj, sim_method=sim_method)
            }
        }
        
        if (save_rds) {
            saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
        }
    }
    

    ## #########################
    ## Step: log transformation

    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 4
    flog.info(sprintf("\n\n\tSTEP %02d: log transformation of data\n", step_count))

    if (skip_past < step_count) {
        SPACEmapX_obj <- log2xplus1(SPACEmapX_obj)

        if (save_rds) {
            saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
        }
        invisible(gc())

        ## Plot incremental steps.
        if (plot_steps){
            plot_cnv(SPACEmapX_obj=SPACEmapX_obj,
                     obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title=sprintf("%02d_log_transformed_data",step_count),
                     output_filename=sprintf("SPACEmapX.%02d_log_transformed",step_count),
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster
                     )
        }
    }
    
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 5
    if (scale_data) {
        
        flog.info(sprintf("\n\n\tSTEP %02d: scaling all expression data\n", step_count))

        if (skip_past < step_count) {
            SPACEmapX_obj <- scale_SPACEmapX_expr(SPACEmapX_obj)
            
            if (save_rds) {
                saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            ## Plot incremental steps.
            if (plot_steps){
                
                plot_cnv(SPACEmapX_obj,
                                              obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_scaled",step_count),
                         output_filename=sprintf("SPACEmapX.%02d_scaled",step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
                
            }
        }
    }
    
    
    ## #################################################
    ## Step: Split the reference data into groups if requested
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 6
    if (!is.null(num_ref_groups)) {
        
        if (! has_reference_cells(SPACEmapX_obj)) {
            stop("Error, no reference cells defined. Cannot split them into groups as requested")
        }
        
        flog.info(sprintf("\n\n\tSTEP %02d: splitting reference data into %d clusters\n", step_count, num_ref_groups))

        if (skip_past < step_count) {
            SPACEmapX_obj <- split_references(SPACEmapX_obj,
                                       num_groups=num_ref_groups,
                                       hclust_method=hclust_method)
            if (save_rds) {
                saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
            }
        }
        
    }
    
    
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 7
    if (analysis_mode == 'subclusters' & tumor_subcluster_partition_method == 'random_trees') {
        
        flog.info(sprintf("\n\n\tSTEP %02d: computing tumor subclusters via %s\n", step_count, tumor_subcluster_partition_method))
        
        resume_file_token = paste0(resume_file_token, ".rand_trees")

        if (skip_past < step_count) {
            SPACEmapX_obj <- define_signif_tumor_subclusters_via_random_smooothed_trees(SPACEmapX_obj,
                                                                                       p_val=tumor_subcluster_pval,
                                                                                       hclust_method=hclust_method,
                                                                                       cluster_by_groups=cluster_by_groups)
            if (save_rds) {
                saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            if (plot_steps) {
                plot_cnv(SPACEmapX_obj,
                                              obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_tumor_subclusters.%s", step_count, tumor_subcluster_partition_method),
                         output_filename=sprintf("SPACEmapX.%02d_tumor_subclusters.%s", step_count, tumor_subcluster_partition_method),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
                
            }

            if (inspect_subclusters) {
                plot_subclusters(SPACEmapX_obj,
                                 out_dir=out_dir,
                                 output_filename="SPACEmapX_subclusters")
            }
        }
    }
    
    
    ## ##################################
    ## Step: Subtract average reference
    ## Since we're in log space, this now becomes log(fold_change)
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 8
    flog.info(sprintf("\n\n\tSTEP %02d: removing average of reference data (before smoothing)\n", step_count))

    if (skip_past < step_count) {
        SPACEmapX_obj <- subtract_ref_expr_from_obs(SPACEmapX_obj, inv_log=FALSE, use_bounds=ref_subtract_use_mean_bounds)
        
        if (save_rds) {
            saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
        }
        invisible(gc())
        
        if (plot_steps) {
            plot_cnv(SPACEmapX_obj,
                                          obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title=sprintf("%02d_remove_average",step_count),
                     output_filename=sprintf("SPACEmapX.%02d_remove_average", step_count),
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster)
        }
    }
    
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 9
    if (! is.na(max_centered_threshold)) {
        
        ## #####################################################
        ## Apply maximum centered expression thresholds to data
        ## Cap values between threshold and -threshold, retaining earlier center_cell_expr_across_chromosome
        
        flog.info(sprintf("\n\n\tSTEP %02d: apply max centered expression threshold: %s\n", step_count, max_centered_threshold))

        if (skip_past < step_count) {
            threshold = max_centered_threshold
            if (is.character(max_centered_threshold) && max_centered_threshold == "auto") {
                threshold = mean(abs(get_average_bounds(SPACEmapX_obj)))
                flog.info(sprintf("Setting max centered threshoolds via auto to: +- %g", threshold))
            }
            
            SPACEmapX_obj <- apply_max_threshold_bounds(SPACEmapX_obj, threshold=threshold)
            
            if (save_rds) {
                saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            ## Plot incremental steps.
            if (plot_steps){
                
                plot_cnv(SPACEmapX_obj,
                                              obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_apply_max_centered_expr_threshold",step_count),
                         output_filename=sprintf("SPACEmapX.%02d_apply_max_centred_expr_threshold",step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
                
            }
        }
    }
    
    
    
    ## #########################################################################
    ## Step: For each cell, smooth the data along chromosome with gene windows
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 10
    flog.info(sprintf("\n\n\tSTEP %02d: Smoothing data per cell by chromosome\n", step_count))
    
    if (skip_past < step_count) {

        if (smooth_method == 'runmeans') {
            
            SPACEmapX_obj <- smooth_by_chromosome_runmeans(SPACEmapX_obj, window_length)
        } else if (smooth_method == 'pyramidinal') {
            
            SPACEmapX_obj <- smooth_by_chromosome(SPACEmapX_obj, window_length=window_length, smooth_ends=TRUE)
        } else if (smooth_method == 'coordinates') {
            SPACEmapX_obj <- smooth_by_chromosome_coordinates(SPACEmapX_obj, window_length=window_length)
        } else {
            stop(sprintf("Error, don't recognize smoothing method: %s", smooth_method))
        }
        
        if (save_rds) {
            saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
        }
        invisible(gc())
        
        ## Plot incremental steps.
        if (plot_steps){
            
            plot_cnv(SPACEmapX_obj,
                                          obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title=sprintf("%02d_smoothed_by_chr",step_count),
                     output_filename=sprintf("SPACEmapX.%02d_smoothed_by_chr", step_count),
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster)
        }
    }
    
    
    ##
    ## Step:
    ## Center cells/observations after smoothing. This helps reduce the
    ## effect of complexity.
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 11
    flog.info(sprintf("\n\n\tSTEP %02d: re-centering data across chromosome after smoothing\n", step_count))
    
    if (skip_past < step_count) {
        SPACEmapX_obj <- center_cell_expr_across_chromosome(SPACEmapX_obj, method="median")
        
        if (save_rds) {
            saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
        }
        invisible(gc())
        
        ## Plot incremental steps.
        if (plot_steps) {
            
            plot_cnv(SPACEmapX_obj,
                                          obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title=sprintf("%02d_centering_of_smoothed",step_count),
                     output_filename=sprintf("SPACEmapX.%02d_centering_of_smoothed", step_count),
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster)
            
        }
    }
    
    
    
    ## ##################################
    ## Step: Subtract average reference (adjustment)
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 12
    flog.info(sprintf("\n\n\tSTEP %02d: removing average of reference data (after smoothing)\n", step_count))
    
    if (skip_past < step_count) {
        SPACEmapX_obj <- subtract_ref_expr_from_obs(SPACEmapX_obj, inv_log=FALSE, use_bounds=ref_subtract_use_mean_bounds)
        
        if (save_rds) {
            saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
        }
        invisible(gc())
        
        if (plot_steps) {
            plot_cnv(SPACEmapX_obj,
                                          obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title=sprintf("%02d_remove_average",step_count),
                     output_filename=sprintf("SPACEmapX.%02d_remove_average", step_count),
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster)
        }
    }
    
    
    ## Step: Remove Ends
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 13
    if (remove_genes_at_chr_ends == TRUE && smooth_method != 'coordinates') {
        
        flog.info(sprintf("\n\n\tSTEP %02d: removing genes at chr ends\n", step_count))
        
        if (skip_past < step_count) {
            SPACEmapX_obj <- remove_genes_at_ends_of_chromosomes(SPACEmapX_obj, window_length)
            
            if (save_rds) {
                saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            ## Plot incremental steps.
            if (plot_steps){
                
                plot_cnv(SPACEmapX_obj,
                                              obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_remove_genes_at_chr_ends",step_count),
                         output_filename=sprintf("SPACEmapX.%02d_remove_genes_at_chr_ends",step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res)
                
            }
        }
    }
    
    
    ## ###########################
    ## Step: invert log transform  (convert from log(FC) to FC)
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 14
    flog.info(sprintf("\n\n\tSTEP %02d: invert log2(FC) to FC\n", step_count))

    if (skip_past < step_count) {
        
        SPACEmapX_obj <- invert_log2(SPACEmapX_obj)
        
        if (save_rds) {
            saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
        }
        invisible(gc())
        
        if (plot_steps) {
            plot_cnv(SPACEmapX_obj,
                                          obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title=sprintf("%02d_invert_log_transform log(FC)->FC",step_count),
                     output_filename=sprintf("SPACEmapX.%02d_invert_log_FC",step_count),
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster)
            
        }
    }
    
    
    ## ###################################################################
    ## Done restoring SPACEmapX_obj's from files now under resume_mode
    ## ###################################################################

    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 15
    if (analysis_mode == 'subclusters' & tumor_subcluster_partition_method != 'random_trees') {
        
        resume_file_token = paste0(resume_file_token, '.', tumor_subcluster_partition_method)

        flog.info(sprintf("\n\n\tSTEP %02d: computing tumor subclusters via %s\n", step_count, tumor_subcluster_partition_method))
        
        if (skip_past < step_count) {
            res <- define_signif_tumor_subclusters(SPACEmapX_obj = SPACEmapX_obj,
                                                   p_val = tumor_subcluster_pval,
                                                   k_nn = k_nn,
                                                   leiden_resolution = leiden_resolution,
                                                   leiden_method = leiden_method,
                                                   leiden_function = leiden_function,
                                                   leiden_method_per_chr = leiden_method_per_chr,
                                                   leiden_function_per_chr = leiden_function_per_chr,
                                                   leiden_resolution_per_chr = leiden_resolution_per_chr,
                                                   hclust_method = hclust_method,
                                                   cluster_by_groups = cluster_by_groups,
                                                   partition_method = tumor_subcluster_partition_method,
                                                   per_chr_hmm_subclusters = per_chr_hmm_subclusters,
                                                   z_score_filter=z_score_filter
                                                   )

            SPACEmapX_obj = res[[1]]
            subclusters_per_chr = res[[2]]
            rm(res)

            if (save_rds) {
                saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
                if(per_chr_hmm_subclusters) {
                    tmp = nchar(reload_info$expected_file_names[[step_count]])
                    tmp = paste0(substr(reload_info$expected_file_names[[step_count]], 1, tmp-13), ".per_chr_subclusters", substr(reload_info$expected_file_names[[step_count]], tmp-12, tmp))
                    saveRDS(subclusters_per_chr, tmp)
                }
            }
            invisible(gc())
            
            if (plot_steps) {
                
                plot_cnv(SPACEmapX_obj,
                                              obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_tumor_subclusters",step_count),
                         output_filename=sprintf("SPACEmapX.%02d_tumor_subclusters",step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
            }

            if (inspect_subclusters) {
                plot_subclusters(SPACEmapX_obj,
                                 out_dir=out_dir,
                                 output_filename="SPACEmapX_subclusters")
            }
        }
    }
    else if (analysis_mode != 'subclusters') {
        
        flog.info(sprintf("\n\n\tSTEP %02d: Clustering samples (not defining tumor subclusters)\n", step_count))
        
        if (skip_past < step_count) {
            
            
            SPACEmapX_obj <- define_signif_tumor_subclusters(SPACEmapX_obj,
                                                            p_val=tumor_subcluster_pval,
                                                            hclust_method=hclust_method,
                                                            cluster_by_groups=cluster_by_groups,
                                                            partition_method='none',
                                                            per_chr_hmm_subclusters=FALSE,
                                                            z_score_filter=z_score_filter
                                                            )[[1]]
            
            if (save_rds) {
                saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
            }
        }
        
    }
    
    
    ## This is a milestone step and results should always be examined here.
    if (skip_past < step_count) {
        if (save_rds) {
            saveRDS(SPACEmapX_obj, file=file.path(out_dir, "preliminary.SPACEmapX_obj"))
        }
        
        invisible(gc())
    
        if (! (no_prelim_plot | no_plot) ) {
            
            #prelim_heatmap_png = "SPACEmapX.preliminary.png"
            
            #if (! file.exists(file.path(out_dir, prelim_heatmap_png))) {
            plot_cnv(SPACEmapX_obj,
                                          obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title="Preliminary SPACEmapX (pre-noise filtering)",
                     output_filename="SPACEmapX.preliminary", # png ext auto added
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster)
            #}
        }
    }
    
    ## Below represent optional downstream analysis steps:
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 16
    if (prune_outliers) {
        
        ## ################################
        ## STEP: Remove outliers for viz
        
        flog.info(sprintf("\n\n\tSTEP %02d: Removing outliers\n", step_count))

        if (skip_past < step_count) {
            
            SPACEmapX_obj = remove_outliers_norm(SPACEmapX_obj,
                                                out_method=outlier_method_bound,
                                                lower_bound=outlier_lower_bound,
                                                upper_bound=outlier_upper_bound)
            
            if (save_rds) {
                saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            ## Plot incremental steps.
            if (plot_steps) {
                
                plot_cnv(SPACEmapX_obj,
                                              obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_removed_outliers",step_count),
                         output_filename=sprintf("SPACEmapX.%02d_removed_outliers", step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
            }
        }
    }
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 17
    hmm_resume_file_token = paste0(resume_file_token, ".hmm_mode-", analysis_mode)
    if (HMM) {

        if (HMM_type == 'i6') {
            hmm_center = 3
            hmm_state_range = c(0,6)
        } else {
            ## i3
            hmm_center = 2
            hmm_state_range = c(1,3)
        }
        
        if (skip_hmm < 1) {
            flog.info(sprintf("\n\n\tSTEP %02d: HMM-based CNV prediction\n", step_count))
            
            if (analysis_mode == 'subclusters') {
                
                if (HMM_type == 'i6') {
                    if (tumor_subcluster_partition_method == "leiden" && per_chr_hmm_subclusters) {
                        hmm.SPACEmapX_obj <- predict_CNV_via_HMM_on_tumor_subclusters_per_chr(SPACEmapX_obj = SPACEmapX_obj, 
                                                                                             subclusters_per_chr = subclusters_per_chr,
                                                                                             t=HMM_transition_prob)
                    }
                    else {
                        hmm.SPACEmapX_obj <- predict_CNV_via_HMM_on_tumor_subclusters(SPACEmapX_obj = SPACEmapX_obj,
                                                                                     t=HMM_transition_prob)
                    }
                } else if (HMM_type == 'i3') {
                    hmm.SPACEmapX_obj <- i3HMM_predict_CNV_via_HMM_on_tumor_subclusters(SPACEmapX_obj,
                                                                                       i3_p_val=HMM_i3_pval,
                                                                                       t=HMM_transition_prob,
                                                                                       use_KS=HMM_i3_use_KS)
                } else {
                    stop("Error, not recognizing HMM_type")
                }

                if (tumor_subcluster_partition_method == 'random_trees') {
                    ## need to redo the hierarchicial clustering, since the subcluster assignments dont always perfectly line up with the top-level dendrogram.
                    hmm.SPACEmapX_obj <- .redo_hierarchical_clustering(hmm.SPACEmapX_obj, hclust_method=hclust_method)
                }
                
            } else if (analysis_mode == 'cells') {
                
                if (HMM_type == 'i6') {
                    hmm.SPACEmapX_obj <- predict_CNV_via_HMM_on_indiv_cells(SPACEmapX_obj, t=HMM_transition_prob)
                } else if (HMM_type == 'i3') {
                    hmm.SPACEmapX_obj <- i3HMM_predict_CNV_via_HMM_on_indiv_cells(SPACEmapX_obj,
                                                                                 i3_p_val=HMM_i3_pval,
                                                                                 t=HMM_transition_prob,
                                                                                 use_KS=HMM_i3_use_KS)
                } else {
                    stop("Error, not recognizing HMM_type")
                }
                
                
            } else {
                ## samples mode
                
                if (HMM_type == 'i6') {
                    hmm.SPACEmapX_obj <- predict_CNV_via_HMM_on_whole_tumor_samples(SPACEmapX_obj,
                                                                                   cluster_by_groups=cluster_by_groups,
                                                                                   t=HMM_transition_prob)
                } else if (HMM_type == 'i3') {
                    hmm.SPACEmapX_obj <- i3HMM_predict_CNV_via_HMM_on_whole_tumor_samples(SPACEmapX_obj,
                                                                                         cluster_by_groups=cluster_by_groups,
                                                                                         i3_p_val=HMM_i3_pval,
                                                                                         t=HMM_transition_prob,
                                                                                         use_KS=HMM_i3_use_KS
                                                                                         )
                } else {
                    stop("Error, not recognizing HMM_type")
                }

            }
            
            
            ## report predicted cnv regions:
            generate_cnv_region_reports(hmm.SPACEmapX_obj,
                                        output_filename_prefix=sprintf("%02d_HMM_pred%s", step_count, hmm_resume_file_token),
                                        out_dir=out_dir,
                                        ignore_neutral_state=hmm_center,
                                        by=HMM_report_by)

            
            ## ##################################
            ## Note, HMM invercnv object is only leveraged here, but stored as file for future use:
            ## ##################################
            
            if (save_rds) {
                saveRDS(hmm.SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
            }

            invisible(gc())
                        
            if (! no_plot) {
                
                ## Plot HMM pred img
                plot_cnv(SPACEmapX_obj=hmm.SPACEmapX_obj,
                                              obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_HMM_preds",step_count),
                         output_filename=sprintf("SPACEmapX.%02d_HMM_pred%s", step_count, hmm_resume_file_token),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         x.center=hmm_center,
                         x.range=hmm_state_range,
                         png_res=png_res,
                         useRaster=useRaster
                         )
            }
        }
    }

        ## ############################################################
        ## Bayesian Network Mixture Model
        ## ############################################################
        
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 18
    if (skip_hmm < 2) {
        if (HMM == TRUE && BayesMaxPNormal > 0 && length(unique(apply(hmm.SPACEmapX_obj@expr.data,2,unique))) != 1 ) {
            flog.info(sprintf("\n\n\tSTEP %02d: Run Bayesian Network Model on HMM predicted CNVs\n", step_count))
            
            ## the MCMC  object
            
            mcmc_obj <- SPACEmapX::SPACEmapXBayesNet( SPACEmapX_obj     = SPACEmapX_obj,
                                                   HMM_states        = hmm.SPACEmapX_obj@expr.data,
                                                   file_dir          = out_dir,
                                                   no_plot           = no_plot,
                                                   postMcmcMethod    = "removeCNV",
                                                   out_dir           = file.path(out_dir, sprintf("BayesNetOutput.%s", hmm_resume_file_token)),
                                                   resume_file_token = hmm_resume_file_token,
                                                   quietly           = TRUE,
                                                   CORES             = num_threads,
                                                   plotingProbs      = plot_probabilities,
                                                   diagnostics       = diagnostics,
                                                   HMM_type          = HMM_type,
                                                   k_obs_groups      = k_obs_groups,
                                                   cluster_by_groups = cluster_by_groups,
                                                   reassignCNVs      = reassignCNVs,
                                                   useRaster         = useRaster)

            # mcmc_obj_file = file.path(out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj",
            #                                        step_count, hmm_resume_file_token))
            
            if (save_rds) {
                # saveRDS(mcmc_obj, file=mcmc_obj_file)
                saveRDS(mcmc_obj, reload_info$expected_file_names[[step_count]])
            }
        }
    }


        ## ############################################################
        ## Bayesian Filtering
        ## ############################################################

    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 19
    if (skip_hmm < 3) {
        if (HMM == TRUE && BayesMaxPNormal > 0 && length(unique(apply(hmm.SPACEmapX_obj@expr.data,2,unique))) != 1 ) {
            flog.info(sprintf("\n\n\tSTEP %02d: Filter HMM predicted CNVs based on the Bayesian Network Model results and BayesMaxPNormal\n", step_count))
            ## Filter CNV's by posterior Probabilities
            mcmc_obj_hmm_states_list <- SPACEmapX::filterHighPNormals( MCMC_SPACEmapX_obj = mcmc_obj,
                                                                     HMM_states         = hmm.SPACEmapX_obj@expr.data, 
                                                                     BayesMaxPNormal    = BayesMaxPNormal,
                                                                     useRaster          = useRaster)
            
            hmm_states_highPnormCNVsRemoved.matrix <- mcmc_obj_hmm_states_list[[2]]

            # replace states
            hmm.SPACEmapX_obj@expr.data <- hmm_states_highPnormCNVsRemoved.matrix
            
            ## Save the MCMC SPACEmapX object
            if (save_rds) {
                saveRDS(hmm.SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            if (! no_plot) {
                ## Plot HMM pred img after cnv removal
                plot_cnv(SPACEmapX_obj=hmm.SPACEmapX_obj,
                                              obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_HMM_preds_Bayes_Net",step_count),
                         output_filename=sprintf("SPACEmapX.%02d_HMM_pred.Bayes_Net.Pnorm_%g",step_count, BayesMaxPNormal),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         x.center=hmm_center,
                         x.range=hmm_state_range,
                         png_res=png_res,
                         useRaster=useRaster
                         )
            }
            ## write the adjusted CNV report files
            ## report predicted cnv regions:
            adjust_genes_regions_report(mcmc_obj_hmm_states_list[[1]],
                                        # input_filename_prefix=sprintf("%02d_HMM_preds", (step_count-1)),
                                        input_filename_prefix=sprintf("%02d_HMM_pred%s", (step_count-2), hmm_resume_file_token),
                                        output_filename_prefix=sprintf("HMM_CNV_predictions.%s.Pnorm_%g", hmm_resume_file_token, BayesMaxPNormal),
                                        out_dir=out_dir)
        }
    }
        
        
        ## convert from states to representative  intensity values
        
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 20
    if (skip_hmm < 4) {
        if (HMM) {
            flog.info(sprintf("\n\n\tSTEP %02d: Converting HMM-based CNV states to repr expr vals\n", step_count))
            
            if (HMM_type == 'i6') {
                hmm.SPACEmapX_obj <- assign_HMM_states_to_proxy_expr_vals(hmm.SPACEmapX_obj)
            } else if (HMM_type == 'i3') {
                hmm.SPACEmapX_obj <- i3HMM_assign_HMM_states_to_proxy_expr_vals(hmm.SPACEmapX_obj)
            }
            
            if (save_rds) {
                saveRDS(hmm.SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            ## Plot HMM pred img
            if (! no_plot) {
                plot_cnv(SPACEmapX_obj=hmm.SPACEmapX_obj,
                                              obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_HMM_preds.repr_intensities",step_count),
                         output_filename=sprintf("SPACEmapX.%02d_HMM_pred%s.Pnorm_%g.repr_intensities", step_count, hmm_resume_file_token, BayesMaxPNormal),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         x.center=1,
                         x.range=c(-1,3),
                         png_res=png_res,
                         useRaster=useRaster
                         )
            }
        }
    }
    
    ## all processes that are alternatives to the HMM prediction wrt DE analysis and/or denoising
    
    ## Step: Filtering significantly DE genes
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 21
    if (mask_nonDE_genes) {
        
        if (!has_reference_cells(SPACEmapX_obj)) {
            stop("Error, cannot mask non-DE genes when there are no normal references set")
        }
        
        flog.info(sprintf("\n\n\tSTEP %02d: Identify and mask non-DE genes\n", step_count))

        if (skip_past < step_count) {

            SPACEmapX_obj <- mask_non_DE_genes_basic(SPACEmapX_obj,
                                                    p_val_thresh=mask_nonDE_pval,
                                                    test.use = test.use,
                                                    center_val=mean(SPACEmapX_obj@expr.data),
                                                    require_DE_all_normals=require_DE_all_normals)
            
            if (save_rds) {
                saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
            }
            
            invisible(gc())
            
            ## Plot incremental steps.
            if (plot_steps) {
                
                plot_cnv(SPACEmapX_obj,
                                              obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_mask_nonDE",step_count),
                         output_filename=sprintf("SPACEmapX.%02d_mask_nonDE", step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
                
            }
        }
    }
    
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(SPACEmapX_obj)
    }
    step_count = step_count + 1 # 22
    if (denoise) {
        
        ## ##############################
        ## Step: de-noising
        
        flog.info(sprintf("\n\n\tSTEP %02d: Denoising\n", step_count))

        if (skip_past < step_count) {
            
            if (! is.na(noise_filter)) {
                
                if (noise_filter > 0) {
                    flog.info(paste("::process_data:Remove noise, noise threshold at: ", noise_filter))
                    SPACEmapX_obj <- clear_noise(SPACEmapX_obj,
                                                threshold=noise_filter,
                                                noise_logistic=noise_logistic)
                }
                else {
                    ## noise == 0 or negative...
                    ## don't remove noise.
                }
                
            }
            else {
                ## default, use quantiles, if NA
                flog.info(paste("::process_data:Remove noise, noise threshold defined via ref mean sd_amplifier: ", sd_amplifier))
                SPACEmapX_obj <- clear_noise_via_ref_mean_sd(SPACEmapX_obj,
                                                            sd_amplifier = sd_amplifier,
                                                            noise_logistic=noise_logistic)
            }
            
            if (save_rds) {
                saveRDS(SPACEmapX_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            if (plot_steps) {
                plot_cnv(SPACEmapX_obj,
                                              obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         color_safe_pal=FALSE,
                         title=sprintf("%02d_denoised", step_count),
                         output_filename=sprintf("SPACEmapX.%02d_denoised", step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
            }
            
        }
    }
    
    if (save_final_rds) {
        saveRDS(SPACEmapX_obj, file=file.path(out_dir, "run.final.SPACEmapX_obj"))
    }
    
    if (! no_plot) {
        if (is.null(final_scale_limits)) {
            final_scale_limits = "auto"
        }
        if (is.null(final_center_val)) {
            final_center_val = 1
        }
        
        
        flog.info("\n\n## Making the final SPACEmapX heatmap ##")
        invisible(gc())
        plot_cnv(SPACEmapX_obj,
                                      obs_title="Observations (Spots)",
                     ref_title="References (Spots)",
                     contig_cex = 2, # 2 is the best choice
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 cluster_references=cluster_references,
                 plot_chr_scale=plot_chr_scale,
                 chr_lengths=chr_lengths,
                 out_dir=out_dir,
                 x.center=final_center_val,
                 x.range=final_scale_limits,
                 title="SPACEmapX",
                 output_filename="SPACEmapX",
                 output_format=output_format,
                 write_expr_matrix=write_expr_matrix,
                 write_phylo=write_phylo,
                 png_res=png_res,
                 useRaster=useRaster)
    }
    
    return(SPACEmapX_obj)
    
}
