color_5 <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")
color_9 <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51", "#a8dadc", "#457b9d", "#1d3557", "#f1faee")

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

# specify user for paths
user <- Sys.info()[["user"]]

if (user == "choiyej") { # local version
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
  root_path <- ""
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else {
  stop("Unknown user: please specify root path for this user.")
}

# ATTEMPT analysis related functions

#----------------------------------------------------------#
# make_subsets
#----------------------------------------------------------#

make_subsets <- function(so, groups, celltype_col = "celltype", prefix = "croc_so_") {
  stopifnot(celltype_col %in% colnames(so@meta.data))
  
  objs <- lapply(names(groups), function(grp) {
    types <- as.character(groups[[grp]])           # <- ensure character vector
    cells <- rownames(so@meta.data)[so@meta.data[[celltype_col]] %in% types]
    if (length(cells) == 0) return(NULL)
    subset(so, cells = cells)
  })
  names(objs) <- names(groups)
  
  # Drop empty groups
  objs <- objs[!vapply(objs, is.null, logical(1))]
  
  # Assign into env
  for (nm in names(objs)) {
    obj_name <- paste0(prefix, tolower(gsub("[^A-Za-z0-9]+", "_", nm)))
    assign(obj_name, objs[[nm]], envir = .GlobalEnv)
  }
  
  # Robust cell counts (always integer)
  n_cells <- vapply(objs, function(o) as.integer(ncol(o)), integer(1))
  
  data.frame(
    group   = names(objs),
    n_cells = n_cells,
    object  = paste0(prefix, tolower(gsub("[^A-Za-z0-9]+","_", names(objs)))),
    row.names = NULL
  )
}

# ---- Analysis Functions ----

#----------------------------------------------------------#
# run_nebula_for_clinvar
#----------------------------------------------------------#

run_nebula_for_clinvar <- function(
    seurat_obj,
    celltype_values,                 # e.g., c("PT-S1/S2","PT-S3","aPT")
    celltype_var = "celltype",
    clinical_var,                    # e.g., "mgfr_jodal_bsa"
    suffix = NULL,                   # e.g., "kpmp"
    exclude_values = NULL,           # e.g., "PT_lowQuality"
    n_hvgs = 2000,
    workers = 50,                    # parallel workers for foreach
    s3,                              # your S3 client object (e.g., `s3`)
    s3_bucket = "attempt",                       # e.g., "attempt"
    s3_key_prefix = "associations/nebula",  # base path
    assay_layer = "counts",          # use "counts" layer (Seurat v5); change to slot="counts" if needed
    id_var = "subject_id",
    offset_var = "pooled_offset",
    treatment_var = "treatment",
    extra_covars = NULL              # character vec of extra covariate names in meta.data (optional)
) {
  # --- packages
  library(Seurat)
  library(Matrix)
  library(nebula)
  library(doParallel)
  library(foreach)
  library(rlang)
  
  # --- sanity checks
  md <- seurat_obj@meta.data
  needed_cols <- c(id_var, treatment_var, clinical_var, offset_var, celltype_var)
  missing_cols <- setdiff(needed_cols, colnames(md))
  if (length(missing_cols)) {
    stop("Missing required meta columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!all(celltype_values %in% unique(md[[celltype_var]]))) {
    warning("Some target celltype values not present in meta; continuing with those that exist.")
  }
  
  message("Subsetting cells where ", celltype_var, " %in% {", paste(celltype_values, collapse = ", "), "}")
  
  # subset to chosen cell types (and optionally drop any excluded labels)
  subset_expr <- md[[celltype_var]] %in% celltype_values
  if (!is.null(exclude_values)) {
    subset_expr <- subset_expr & !(md[[celltype_var]] %in% exclude_values)
  }
  so_sub <- seurat_obj[, subset_expr]
  
  # HVGs
  so_sub <- FindVariableFeatures(so_sub, selection.method = "vst", nfeatures = n_hvgs)
  hvgs <- VariableFeatures(so_sub)
  so_hvg <- subset(so_sub, features = hvgs)
  
  # keep only rows with non-missing clinical var
  md_hvg <- so_hvg@meta.data
  keep_cells <- rownames(md_hvg)[!is.na(md_hvg[[clinical_var]])]
  so_hvg <- so_hvg[, keep_cells, drop = FALSE]
  
  # counts matrix (Seurat v5 “layer” API)
  counts <- round(GetAssayData(so_hvg, layer = assay_layer))
  genes <- rownames(counts)
  
  # design matrix terms: clinical_var + treatment + optional extras
  rhs_terms <- c(clinical_var, treatment_var, extra_covars)
  rhs <- paste(rhs_terms, collapse = " + ")
  fml <- as.formula(paste("~", rhs))
  
  # parallel
  workers <- max(1, min(workers, parallel::detectCores(logical = TRUE)))
  cl <- parallel::makeCluster(workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  registerDoParallel(cl)
  
  message("Running nebula on ", length(genes), " genes with ", workers, " workers...")
  start_time <- Sys.time()
  
  res_list <- foreach(g = genes, .packages = c("nebula", "Matrix", "Seurat")) %dopar% {
    warn <- NULL; err <- NULL; res <- NULL
    # per-gene meta (keeps rows aligned to columns of count matrix)
    meta_g <- subset(so_hvg, features = g)@meta.data
    pred <- tryCatch(model.matrix(fml, data = meta_g), error = function(e) { err <<- conditionMessage(e); NULL })
    if (is.null(pred)) {
      return(list(gene = g, result = NULL, warning = warn, error = err))
    }
    
    # build grouped count structure
    count_g <- counts[g, , drop = FALSE]
    dat <- tryCatch({
      group_cell(count = count_g, id = meta_g[[id_var]], pred = pred)
    }, error = function(e) { err <<- conditionMessage(e); NULL })
    
    if (is.null(dat)) {
      return(list(gene = g, result = NULL, warning = warn, error = err))
    }
    
    # run nebula
    tryCatch({
      res <- withCallingHandlers({
        nebula(
          count = dat$count,
          id    = dat$id,
          pred  = dat$pred,
          ncore = 1,
          output_re = TRUE,
          covariance = TRUE,
          reml = 1,
          model = "NBLMM",
          offset = meta_g[[offset_var]]
        )
      }, warning = function(w) { warn <<- conditionMessage(w); invokeRestart("muffleWarning") })
    }, error = function(e) { err <<- conditionMessage(e); res <- NULL })
    
    list(gene = g, result = res, warning = warn, error = err)
  }
  
  # collect logs
  for (x in res_list) {
    if (!is.null(x$warning)) message(sprintf("[WARN] %s: %s", x$gene, x$warning))
    if (!is.null(x$error))   message(sprintf("[ERR ] %s: %s", x$gene, x$error))
  }
  
  names(res_list) <- vapply(res_list, function(x) x$gene, "", USE.NAMES = FALSE)
  res_clean <- lapply(res_list, `[[`, "result")
  res_clean <- Filter(Negate(is.null), res_clean)
  
  filtered_prop <- (length(genes) - length(res_clean)) / length(genes)
  message(sprintf("Filtered due to errors/low expr: %.2f%%", 100 * filtered_prop))
  
  # save to S3
  group_tag <- paste0(gsub("[^A-Za-z0-9]+", "_", paste(unique(celltype_values), collapse = "-")))
  file_stub <- paste0(clinical_var, "__", group_tag, if (!is.null(suffix)) paste0("__", suffix))
  out_key <- file.path(s3_key_prefix, clinical_var, paste0(file_stub, ".rds"))
  tmp <- tempfile(fileext = ".rds")
  saveRDS(res_clean, tmp)
  s3$upload_file(tmp, s3_bucket, out_key)
  unlink(tmp)
  
  total_min <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2)
  message("Done: ", clinical_var, " [", group_tag, "] in ", total_min, " min; saved to s3://", s3_bucket, "/", out_key)
  
  # return both results and where it went
  list(
    clinical_var = clinical_var,
    celltypes = celltype_values,
    n_genes = length(genes),
    kept_genes = length(res_clean),
    filtered_prop = filtered_prop,
    s3_bucket = s3_bucket,
    s3_key = out_key,
    results = res_clean
  )
}

#----------------------------------------------------------#
# run_nebula_by_groups
#----------------------------------------------------------#

run_nebula_by_groups <- function(
    seurat_obj,
    celltype_groups,        # named list: group -> vector of celltype values
    celltype_var = "celltype",
    clin_vars,              # character vec of clinical vars
    suffix = NULL,
    exclude_values = NULL,  # vector of labels to drop (e.g., "PT_lowQuality")
    n_hvgs = 2000,
    workers = 50,
    s3, s3_bucket, s3_key_prefix = "associations/nebula",
    assay_layer = "counts",
    id_var = "subject_id",
    offset_var = "pooled_offset",
    treatment_var = "treatment",
    extra_covars = NULL
) {
  out <- list()
  for (grp in names(celltype_groups)) {
    ct_vals <- celltype_groups[[grp]]
    for (cv in clin_vars) {
      message("=== Group: ", grp, " | ClinVar: ", cv, " ===")
      res <- run_nebula_for_clinvar(
        seurat_obj = seurat_obj,
        celltype_values = ct_vals,
        celltype_var = celltype_var,
        clinical_var = cv,
        suffix = suffix,
        exclude_values = exclude_values,
        n_hvgs = n_hvgs,
        workers = workers,
        s3 = s3,
        s3_bucket = s3_bucket,
        s3_key_prefix = s3_key_prefix,
        assay_layer = assay_layer,
        id_var = id_var,
        offset_var = offset_var,
        treatment_var = treatment_var,
        extra_covars = extra_covars
      )
      if (is.null(out[[grp]])) out[[grp]] <- list()
      out[[grp]][[cv]] <- res
    }
  }
  out
}



#----------------------------------------------------------#
# run_ancova
#----------------------------------------------------------#

run_ancova <- function(data, outcome_var, baseline_var, visit_week) {
  outcome_sym <- rlang::sym(outcome_var)
  baseline_sym <- rlang::sym(baseline_var)
  
  data_filtered <- data %>% 
    filter(visit == visit_week)
  
  ancova_results <- data_filtered %>%
    nest_by(visit) %>%
    dplyr::mutate(
      model = list(lm(!!outcome_sym ~ treatment + !!baseline_sym, data = data_filtered)),
      tidy_results = list(broom::tidy(model)),
      glance_results = list(broom::glance(model)),
      emm = list(emmeans::emmeans(model, "treatment")),
      emm_results = list(broom::tidy(emm, conf.int = TRUE)),
      contrasts = list(pairs(emm)),
      contrast_results = list(broom::tidy(contrasts))
    )
  
  ancova_summary <- ancova_results %>%
    dplyr::select(visit, tidy_results) %>%
    unnest(tidy_results) %>%
    filter(term == "treatmentDapagliflozin 5mg") %>%
    dplyr::select(visit, estimate, std.error, p.value, statistic)
  
  ancova_means <- ancova_results %>%
    dplyr::select(visit, emm_results) %>%
    unnest(emm_results)
  
  ancova_contrasts <- ancova_results %>%
    dplyr::select(visit, contrast_results) %>%
    unnest(contrast_results)
  
  return(list(summary = ancova_summary,
              means = ancova_means,
              contrasts = ancova_contrasts))
}


#----------------------------------------------------------#
# run_doubletfinder_custom
#----------------------------------------------------------#
# run_doubletfinder_custom runs Doublet_Finder() and returns a dataframe with the cell IDs and a column with either 'Singlet' or 'Doublet'
run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL){
  # for debug
  #seu_sample_subset <- samp_split[[1]]
  # Print sample number
  print(paste0("Sample ", unique(seu_sample_subset[['SampleID']]), '...........')) 
  
  if(is.null(multiplet_rate)){
    print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
    
    print(multiplet_rates_10x)
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
  }
  
  # Pre-process seurat object with standard seurat workflow --- 
  sample <- NormalizeData(seu_sample_subset)
  sample <- FindVariableFeatures(sample)
  sample <- ScaleData(sample)
  sample <- RunPCA(sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  
  # Finish pre-processing with min_pc
  sample <- RunUMAP(sample, dims = 1:min_pc)
  sample <- FindNeighbors(object = sample, dims = 1:min_pc)              
  sample <- FindClusters(object = sample, resolution = 0.1)
  
  # pK identification (no ground-truth) 
  #introduces artificial doublets in varying props, merges with real data set and 
  # preprocesses the data + calculates the prop of artficial neighrest neighbours, 
  # provides a list of the proportion of artificial nearest neighbours for varying
  # combinations of the pN and pK
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  ## Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
  homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
  
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  # run DoubletFinder
  sample <- doubletFinder(seu = sample, 
                          PCs = 1:min_pc, 
                          pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                          nExp = nExp.poi.adj) # number of expected real doublets
  # change name of metadata column with Singlet/Doublet information
  colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
  
  # Subset and save
  # head(sample@meta.data['doublet_finder'])
  # singlets <- subset(sample, doublet_finder == "Singlet") # extract only singlets
  # singlets$ident
  double_finder_res <- sample@meta.data['doublet_finder'] # get the metadata column with singlet, doublet info
  double_finder_res <- rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly
  return(double_finder_res)
}


# ===========================================================================
# Function: run_nebula_parallel
# ===========================================================================

run_nebula_parallel <- function(seurat_obj,
                                n_cores = 100,
                                layer = "counts",
                                subject_id_col = "subject_id",
                                offset_col = "pooled_offset",
                                formula = ~ treatment * visit,
                                model = "NBLMM",
                                reml = 1,
                                output_re = TRUE,
                                covariance = TRUE,
                                s3_bucket = "attempt",
                                s3_key = NULL,
                                verbose = TRUE, 
                                group = T) {
  
  # Extract counts and gene list
  counts_mat <- round(GetAssayData(seurat_obj, layer = layer))
  genes_list <- rownames(counts_mat)
  
  # Set up parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Ensure cleanup on exit
  on.exit({
    stopCluster(cl)
  }, add = TRUE)
  
  start_time <- Sys.time()
  
  # Run nebula in parallel
  nebula_results_list <- foreach(g = genes_list, 
                                 .packages = c("nebula", "Matrix"),
                                 .errorhandling = "pass") %dopar% {
                                   warn <- NULL
                                   err <- NULL
                                   res <- NULL
                                   
                                   tryCatch({
                                     # Subset data for single gene
                                     count_gene <- counts_mat[g, , drop = FALSE]
                                     meta_gene <- subset(seurat_obj, features = g)@meta.data
                                     
                                     # Create model matrix
                                     pred_gene <- model.matrix(formula, data = meta_gene)
                                     
                                     # Group cells
                                     if (group) {
                                       data_g_gene <- group_cell(count = count_gene, 
                                                                 id = meta_gene[[subject_id_col]], 
                                                                 pred = pred_gene)
                                     } else {
                                       data_g_gene <- list(count = count_gene, 
                                                           id = meta_gene[[subject_id_col]], 
                                                           pred = pred_gene)
                                     }
                                     
                                     # Run nebula with warning handling
                                     res <- withCallingHandlers({
                                       nebula(count = data_g_gene$count, 
                                              id = data_g_gene$id, 
                                              pred = data_g_gene$pred,
                                              ncore = 1, 
                                              output_re = output_re, 
                                              covariance = covariance, 
                                              reml = reml, 
                                              model = model, 
                                              offset = meta_gene[[offset_col]])
                                     }, warning = function(w) {
                                       warn <<- conditionMessage(w)
                                       invokeRestart("muffleWarning")
                                     })
                                   }, error = function(e) {
                                     err <<- conditionMessage(e)
                                   })
                                   
                                   list(gene = g, result = res, warning = warn, error = err)
                                 }
  
  # Report warnings and errors if verbose
  if (verbose) {
    for (res in nebula_results_list) {
      if (!is.null(res$warning)) {
        cat(sprintf("Warning for gene %s: %s\n", res$gene, res$warning))
      }
      if (!is.null(res$error)) {
        cat(sprintf("Error for gene %s: %s\n", res$gene, res$error))
      }
    }
  }
  
  # Clean up results
  names(nebula_results_list) <- sapply(nebula_results_list, function(x) x$gene)
  nebula_results_list <- lapply(nebula_results_list, function(x) x$result)
  nebula_results_list <- Filter(Negate(is.null), nebula_results_list)
  
  # Calculate runtime and non-convergence rate
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  
  nebula_nonconverged_percent <- (length(genes_list) - length(nebula_results_list)) / length(genes_list)
  
  if (verbose) {
    cat(sprintf("\nRuntime: %.2f minutes\n", runtime))
    cat(sprintf("%.2f%% of genes filtered due to low expression or convergence issues\n", 
                nebula_nonconverged_percent * 100))
  }
  
  # Save to S3 if requested
  if (!is.null(s3_key) && !is.null(s3_bucket)) {
    temp_file <- tempfile(fileext = ".rds")
    saveRDS(nebula_results_list, temp_file)
    s3$upload_file(temp_file, s3_bucket, s3_key)
    unlink(temp_file)  # Clean up temp file
    
    if (verbose) {
      cat(sprintf("Results uploaded to s3://%s/%s\n", s3_bucket, s3_key))
    }
  }
  
  # Return results and metadata
  return(list(
    results = nebula_results_list,
    n_genes_tested = length(genes_list),
    n_genes_converged = length(nebula_results_list),
    nonconverged_percent = nebula_nonconverged_percent,
    runtime_minutes = as.numeric(runtime)
  ))
}

# ===========================================================================
# Function: process_nebula_results
# ===========================================================================

process_nebula_results <- function(nebula_list, 
                                   pval_col = "p_treatmentDapagliflozin:visitPOST", 
                                   convergence_cut = -10) {
  # Extract convergence codes
  convergence_df <- purrr::map_dfr(names(nebula_list), function(gene_name) {
    convergence_code <- nebula_list[[gene_name]]$convergence
    data.frame(Gene = gene_name, Convergence_Code = convergence_code)
  })
  
  # Filter to converged models
  converged_genes <- convergence_df %>%
    filter(Convergence_Code >= convergence_cut) %>%
    pull(Gene)
  
  # Combine model summary results
  summary_df <- purrr::map_dfr(converged_genes, function(gene_name) {
    nebula_list[[gene_name]]$summary %>%
      dplyr::mutate(Gene = gene_name)
  })
  
  # Add FDR adjustment
  if (pval_col %in% names(summary_df)) {
    summary_df <- summary_df %>%
      dplyr::mutate(fdr = p.adjust(.data[[pval_col]], method = "fdr"))
  } else {
    warning(paste("Column", pval_col, "not found in summary data. FDR not computed."))
    summary_df$fdr <- NA
  }
  
  # Extract overdispersion estimates
  overdisp_df <- purrr::map_dfr(names(nebula_list), function(gene_name) {
    od <- nebula_list[[gene_name]]$overdispersion
    od$Gene <- gene_name
    od
  })
  
  return(list(
    convergence     = convergence_df,
    results         = summary_df,
    overdispersion  = overdisp_df
  ))
}

# ===========================================================================
# Function: process_nebula_results_clin
# ===========================================================================

process_nebula_results_clin <- function(cell_type, clinical_var, cell_subtype,
                                        bucket = "attempt", region = "",
                                        output_dir = "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/ATTEMPT/Results/nebula/") {
  
  # Construct file path
  clean_subtype <- paste(
    gsub("[/\\-]", "_",                      # replace / and - with _
         gsub("\\+", "",                   # remove +
              gsub("\\s+", "_", cell_subtype) # replace spaces with _
         )
    ),
    collapse = "_"
  )
  
  cat(clean_subtype)
  file_path <- paste0('associations/nebula/', clinical_var, '/', clinical_var, '__', clean_subtype, '.rds')
  
  # Read in nebula results
  nebula_res <- s3readRDS(file_path, bucket = bucket, region = region)
  
  # Check convergence
  convergence_df <- map_dfr(
    names(nebula_res),
    function(gene_name) {
      converged <- nebula_res[[gene_name]]$convergence
      df <- data.frame(Gene = gene_name,
                       Convergence_Code = converged)
      return(df)
    }
  )
  
  # Filter for converged genes
  converged_genes <- convergence_df %>%
    filter(Convergence_Code >= -10)
  
  # Combine results for converged genes
  res_combined <- map_dfr(
    names(nebula_res),
    function(gene_name) {
      df <- nebula_res[[gene_name]]$summary
      df <- df %>% 
        filter(gene_name %in% converged_genes$Gene) %>%
        dplyr::mutate(Gene = gene_name)
      return(df)
    }
  )
  
  # Adjust for multiple testing
  # Create dynamic column names based on clinical variable
  p_col <- paste0("p_", clinical_var)
  
  res_combined <- res_combined %>%
    ungroup() %>%
    dplyr::mutate(fdr_interaction = p.adjust(!!sym(p_col), method = "fdr"))
  
  # Extract overdispersion information
  res_combined_disp <- map_dfr(
    names(nebula_res),
    function(gene_name) {
      df <- nebula_res[[gene_name]]$overdispersion
      df <- df %>%dplyr::mutate(Gene = gene_name)
      return(df)
    }
  )
  
  # Add treatment direction columns if the data contains treatment interaction terms
  if ("p_treatmentDapagliflozin:visitPOST" %in% names(res_combined)) {
    res_combined <- res_combined %>%
      dplyr::mutate(
        treatment_direction = case_when(
          `p_treatmentDapagliflozin:visitPOST` < 0.05 &
            `logFC_treatmentDapagliflozin:visitPOST` > 0 ~ "Positive",
          `p_treatmentDapagliflozin:visitPOST` < 0.05 &
            `logFC_treatmentDapagliflozin:visitPOST` < 0 ~ "Negative",  # Fixed: should be < 0
          TRUE ~ "NS"
        ),
        treatment_direction_fdr = case_when(
          `fdr_interaction` < 0.05 & 
            `logFC_treatmentDapagliflozin:visitPOST` > 0 ~ "Positive",
          `fdr_interaction` < 0.05 & 
            `logFC_treatmentDapagliflozin:visitPOST` < 0 ~ "Negative",  # Fixed: should be < 0
          TRUE ~ "NS"
        )
      )
  }
  
  # Save results
  output_file <- file.path(output_dir, paste0(cell_type, "_", clinical_var, ".csv"))
  write.csv(res_combined, output_file, row.names = FALSE)
  
  # Return a list with all results
  return(list(
    results_combined = res_combined,
    convergence = convergence_df,
    converged_genes = converged_genes,
    overdispersion = res_combined_disp
  ))
}

# Example usage:
# pt_mgfr_jodal_results <- process_nebula_results(cell_type = "pt", clinical_var = "mgfr_jodal")
# cd4_glucose_results <- process_nebula_results(cell_type = "cd4", clinical_var = "glucose_level")


# ===========================================================================
# Function: run_nebula_attempt
# ===========================================================================

run_nebula_attempt <- function(so,
                               trait            = "",
                               extra_covars     = "treatment",
                               subject_var      = "subject_id",
                               offset_var       = "pooled_offset",
                               assay_layer      = "counts",
                               n_cores          = max(parallel::detectCores() - 1, 1),
                               aws_s3           = NULL,
                               s3_bucket        = NULL,
                               s3_key           = NULL) {
  
  stopifnot(trait %in% colnames(so@meta.data))
  
  # ── 1. Keep only cells with non-missing trait ──────────────────────────────
  md            <- so@meta.data
  keep_cells    <- rownames(md)[!is.na(md[[trait]])]
  so_subset     <- so[, keep_cells]
  
  # ── 2. Pull counts and gene list ──────────────────────────────────────────
  counts_mat    <- round(GetAssayData(so_subset, layer = assay_layer))
  genes_list    <- rownames(counts_mat)
  
  # ── 3. Spin up parallel backend ───────────────────────────────────────────
  cl            <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  on.exit({          # make *sure* we clean up
    parallel::stopCluster(cl)
  }, add = TRUE)
  
  start_time <- Sys.time()
  
  # ── 4. Per-gene nebula fits ───────────────────────────────────────────────
  nebula_res <- foreach::foreach(
    g = genes_list,
    .packages      = c("nebula", "Matrix"),
    .errorhandling = "pass"
  ) %dopar% {
    
    warn <- err <- NULL
    res  <- NULL
    
    tryCatch({
      count_gene <- counts_mat[g, , drop = FALSE]
      meta_gene  <- subset(so_subset, features = g)@meta.data
      
      pred_formula <- reformulate(c(trait, extra_covars), response = NULL) # ~ trait  extras
      pred_gene    <- model.matrix(pred_formula, data = meta_gene)
      
      data_g       <- list(count = count_gene,
                           id    = meta_gene[[subject_var]],
                           pred  = pred_gene)
      
      res <- withCallingHandlers(
        nebula::nebula(count      = data_g$count,
                       id         = data_g$id,
                       pred       = data_g$pred,
                       model      = "NBLMM",
                       output_re  = TRUE,
                       covariance = TRUE,
                       reml       = TRUE,
                       offset     = if (!is.null(offset_var)) meta_gene[[offset_var]] else NULL,
                       ncore      = 1),
        warning = function(w) { warn <<- conditionMessage(w); invokeRestart("muffleWarning") }
      )
      
    }, error = function(e) {
      err <<- conditionMessage(e)
    })
    
    list(gene = g, result = res, warning = warn, error = err)
  }
  
  # ── 5. Collate warnings / errors ──────────────────────────────────────────
  for (x in nebula_res) {
    if (!is.null(x$warning)) message(sprintf("⚠️  Warning for %s: %s", x$gene, x$warning))
    if (!is.null(x$error))   message(sprintf("⛔ Error   for %s: %s", x$gene, x$error))
  }
  
  # ── 6. Keep successful fits only & name the list ──────────────────────────
  fits <- lapply(nebula_res, `[[`, "result")
  names(fits) <- vapply(nebula_res, `[[`, "", "gene")
  fits        <- Filter(Negate(is.null), fits)
  
  # ── 7. Report & (optionally) persist to S3 ────────────────────────────────
  drop_pct <- 100 * (1 - length(fits) / length(genes_list))
  message(sprintf("%0.2f%% of genes were dropped (low expression / errors).", drop_pct))
  
  if (!is.null(aws_s3) && !is.null(s3_bucket) && !is.null(s3_key)) {
    tmp <- tempfile(fileext = ".rds")
    saveRDS(fits, tmp)
    aws_s3$upload_file(tmp, Bucket = s3_bucket, Key = s3_key)
    unlink(tmp)
    message(sprintf("⬆️  Results uploaded to s3://%s/%s", s3_bucket, s3_key))
  }
  
  end_time <- Sys.time()
  message(sprintf("Finished in %.1f minutes.", as.numeric(difftime(end_time, start_time, units = "mins"))))
  
  invisible(fits)
}

# ---- Plot Functions ----

# ===========================================================================
# Function: plot_delta_by_category
# ===========================================================================

plot_delta_by_category <- function(data,
                                   id_var = "subject_id",
                                   treatment_var = "treatment",
                                   value_var = "value",
                                   x_var = NULL,  # NEW: variable to bin on (optional)
                                   visit_var = "visit",
                                   visit_baseline,
                                   visit_pre,
                                   visit_post,
                                   bin_by_tertiles = FALSE,
                                   bin_breaks = c(-Inf, 7.5, 8, Inf),
                                   bin_labels = c("T1 (lowest)", "T2", "T3 (highest)"),
                                   bin_var_name = "value_bin",
                                   proper_name = "Value",
                                   y_label = NULL,
                                   x_label = NULL,
                                   treatment_levels = c("Dapagliflozin 5mg", "Placebo"),
                                   fill_colors = c("Placebo" = "#f8ae9d", 
                                                   "Dapagliflozin 5mg" = "#a7b298",
                                                   "DiD" = "#669bbc"),
                                   y_limits = NULL,
                                   output_path = NULL,
                                   top_padding = 0) {
  
  library(dplyr)
  library(ggplot2)
  library(ggbeeswarm)
  library(rlang)
  
  if (is.null(y_label)) y_label <- paste0("\u0394 ", proper_name)
  if (is.null(x_label)) x_label <- paste0("\nBaseline ", proper_name)
  
  value_sym <- sym(value_var)
  x_sym <- if (!is.null(x_var)) sym(x_var) else value_sym
  
  # Step 1: Prepare raw data
  summary_df <- data %>%
    group_by(.data[[id_var]]) %>%
    dplyr::mutate(
      baseline_value = .data[[value_var]][.data[[visit_var]] == visit_baseline],
      baseline_x_value = .data[[as_string(x_sym)]][.data[[visit_var]] == visit_baseline],
      pre_value = .data[[value_var]][.data[[visit_var]] == visit_pre],
      post_value = .data[[value_var]][.data[[visit_var]] == visit_post],
      delta_value = post_value - pre_value
    ) %>%
    ungroup() %>%
    dplyr::select(-.data[[visit_var]], -.data[[value_var]]) %>%
    distinct(.data[[id_var]], .keep_all = TRUE) %>%
    filter(!is.na(delta_value), !is.na(baseline_x_value))
  
  # Step 2: Binning
  if (bin_by_tertiles) {
    tertiles <- quantile(summary_df$baseline_x_value, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
    if (length(bin_labels) != 3) {
      bin_labels <- c("T1 (lowest)", "T2", "T3 (highest)")
    }
    summary_df[[bin_var_name]] <- cut(
      summary_df$baseline_x_value,
      breaks = tertiles,
      include.lowest = TRUE
    )
  } else {
    summary_df[[bin_var_name]] <- cut(
      summary_df$baseline_x_value,
      breaks = bin_breaks,
      labels = bin_labels,
      right = FALSE
    )
  }
  
  summary_df <- summary_df %>%
    filter(!is.na(.data[[bin_var_name]]))
  
  # Step 3: Summary for bars
  plot_df <- summary_df %>%
    group_by(.data[[bin_var_name]], .data[[treatment_var]]) %>%
    dplyr::summarise(
      mean_delta = mean(delta_value, na.rm = TRUE),
      se = sd(delta_value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      ci_lower = mean_delta - 1.96 * se,
      ci_upper = mean_delta + 1.96 * se
    )
  
  
  # Step 4: Paired t-test per bin × treatment
  sig_df <- summary_df %>%
    group_by(.data[[bin_var_name]], .data[[treatment_var]]) %>%
    dplyr::summarise(
      p_val = tryCatch(t.test(post_value, pre_value, paired = TRUE)$p.value, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      stars = case_when(
        is.na(p_val) ~ "",
        p_val < 0.001 ~ "***",
        p_val < 0.01 ~ "**",
        p_val < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  plot_df <- plot_df %>%
    left_join(sig_df, by = c(bin_var_name, treatment_var)) %>%
    dplyr::mutate(y_position = ifelse(ci_upper > 0, ci_upper + 0.5, 0.5)) 
  
  # Adding DiD bars
  # Ensure treatment factor levels (matches function's colors)
  summary_df[[treatment_var]] <- factor(summary_df[[treatment_var]], levels = treatment_levels)
  
  # Model and pairwise contrasts by bin
  # Dynamically build formula
  model_formula <- as.formula(
    paste("delta_value ~", treatment_var, "*", bin_var_name)
  )
  
  # Fit model and compute pairwise comparisons
  bin_contrasts <- pairs(
    emmeans(
      lm(model_formula, data = summary_df),
      specs = as.formula(paste0("~", treatment_var, "|", bin_var_name))
    )
  )
  
  # Convert to dataframe and add CI + stars
  bin_contrasts_df <- as.data.frame(bin_contrasts) %>%
    dplyr::mutate(
      ci_lower = estimate - 1.96 * SE,
      ci_upper = estimate + 1.96 * SE,
      stars = case_when(
        is.na(p.value) ~ "",
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""),
      contrast = "DiD",
      y_position = ifelse(ci_upper > 0, ci_upper + 0.5, 0.5)
    ) %>%
    dplyr::select(all_of(bin_var_name), contrast, estimate, SE, ci_lower, ci_upper, p.value, stars, y_position)
  
  # Standardize names to match plotting code expectations
  names(bin_contrasts_df) <- c(
    bin_var_name,  # dynamic name for binning variable
    "treatment",   # contrast column
    "mean_delta",  # estimate
    "se", 
    "ci_lower", 
    "ci_upper", 
    "p_val", 
    "stars", 
    "y_position"
  )  
  
  # join contrast df with rest
  plot_df <- plot_df %>%
    rbind(bin_contrasts_df)
  
  # Step 5: Fill colors
  if (is.null(fill_colors)) {
    treatments <- unique(summary_df[[treatment_var]])
    fill_colors <- setNames(RColorBrewer::brewer.pal(length(treatments), "Dark2"), treatments)
  }
  
  # Step 6: Plot
  # Combine: DiD first, then treatment_levels
  
  final_levels <- c("DiD", rev(treatment_levels))
  
  # Apply the factor
  plot_df[[treatment_var]] <- factor(plot_df[[treatment_var]], 
                                     levels = final_levels)
  
  p <- ggplot(plot_df, aes(x = .data[[bin_var_name]], y = mean_delta, fill = .data[[treatment_var]])) +
    geom_col(color = "black") +
    geom_quasirandom(
      data = summary_df,
      aes(x = .data[[bin_var_name]], y = delta_value),
      width = 0.15, alpha = 0.15, size = 2
    ) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, size = 0.5) +
    geom_hline(yintercept = 0, size = 0.2) +
    geom_text(aes(label = stars, y = y_position), size = 7, vjust = 0) +
    facet_wrap(as.formula(paste("~", treatment_var))) +
    theme_minimal(base_size = 15) +
    theme(
      axis.text.x = element_text(angle = ifelse(bin_by_tertiles, 60, 0),
                                 hjust = ifelse(bin_by_tertiles, 1, 0)),
      axis.ticks = element_line(color = "grey"),
      panel.grid = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = x_label, y = y_label) +
    scale_fill_manual(values = fill_colors)
  
  # Adjust y-limits for padding
  if (is.null(y_limits)) {
    y_min <- min(plot_df$ci_lower, na.rm = TRUE)
    y_max <- max(plot_df$y_position, na.rm = TRUE)
    extra_space <- ifelse(top_padding <= 1, y_max * top_padding, top_padding) # percent or fixed
    p <- p + ylim(y_min, y_max + extra_space)
  } else {
    # if user also sets y_limits, just add padding to top
    p <- p + ylim(y_limits[1], y_limits[2] + top_padding)
  }
  
  # Output
  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 7, height = 5)
    return(p)
  } else {
    return(p)
  }
}

# ===========================================================================
# Function: make_comp_plot
# ===========================================================================

make_comp_plot <- function(attempt_df,
                           croc_df,
                           attempt_p_cut = 0.05,
                           croc_p_cut   = 0.05,
                           save_path    = NULL,
                           csv_path     = NULL,  
                           width        = 10,
                           height       = 5,
                           caption      = NULL,
                           FC = "logFC_treatmentDapagliflozin:visitPOST") {
  
  library(ggplot2)
  library(ggrounded)
  library(ggtext)
  
  ## 1. Harmonize & flag significance ------------------------------------------
  df_attempt <- attempt_df %>%
    dplyr::rename(p_fdr = fdr,
                  FC    = FC) %>%
    dplyr::mutate(direction = case_when(p_fdr < attempt_p_cut & FC < 0 ~ "-",
                                        p_fdr < attempt_p_cut & FC > 0 ~ "+"),
                  source = "ATTEMPT") %>%
    dplyr::select(Gene, direction, p_fdr, source, FC)
  
  df_croc <- croc_df %>%
    dplyr::rename(p_fdr = p_groupType_1_Diabetes,
                  FC    = logFC_groupType_1_Diabetes) %>%
    dplyr::mutate(direction = case_when(p_fdr < croc_p_cut & FC < 0 ~ "-",
                                        p_fdr < croc_p_cut & FC > 0 ~ "+"),
                  source = "CROCODILE") %>%
    dplyr::select(Gene, direction, p_fdr, source, FC)
  
  pt_comp_df <- bind_rows(df_attempt, df_croc) %>%
    filter(!is.na(direction))
  
  ## 2. Wide table → effect category ------------------------------------------
  pt_comp_df_wide <- pt_comp_df %>%
    dplyr::select(Gene, source, direction) %>%
    pivot_wider(names_from  = source,
                values_from = direction) %>%
    filter(!is.na(ATTEMPT) & !is.na(CROCODILE)) %>%
    dplyr::mutate(effect = case_when(
      (ATTEMPT == "+" & CROCODILE == "+") | 
        (ATTEMPT == "-" & CROCODILE == "-") ~ "Inconsistent",
      ATTEMPT == "+" & CROCODILE == "-" ~ "Reversed towards +\n(T1D depleted)",
      ATTEMPT == "-" & CROCODILE == "+" ~ "Reversed towards -\n(T1D elevated)"
    )) %>%
    dplyr::mutate(effect = factor(effect,
                                  levels = c("Inconsistent",
                                             "Reversed towards +\n(T1D depleted)",
                                             "Reversed towards -\n(T1D elevated)")))
  
  ## Count inconsistent genes for caption
  n_inconsistent <- sum(pt_comp_df_wide$effect == "Inconsistent", na.rm = TRUE)
  n_consistent <- sum(pt_comp_df_wide$effect != "Inconsistent", na.rm = TRUE)
  
  ## Save inconsistent genes to CSV if requested
  if (!is.null(csv_path)) {
    inconsistent_genes <- pt_comp_df_wide %>%
      filter(effect == "Inconsistent") %>%
      dplyr::mutate(consistency_type = case_when(
        ATTEMPT == "+" & CROCODILE == "+" ~ "Both_upregulated",
        ATTEMPT == "-" & CROCODILE == "-" ~ "Both_downregulated"
      )) %>%
      dplyr::arrange(Gene)
    
    write.csv(inconsistent_genes, csv_path, row.names = FALSE)
  }
  
  ## 3. Bar plot ---------------------------------------------------------------
  effect_cols <- c("Reversed towards +\n(T1D depleted)" = "#3a5a40",
                   "Reversed towards -\n(T1D elevated)" = "#a3b18a")
  
  # Update caption to include inconsistent count
  if (is.null(caption)) {
    caption <- paste0("Number of non-reversed genes: ", n_inconsistent,
                      "\nNumber of reversed genes: ", n_consistent)
  } else {
    caption <- paste0(caption, "\nNumber of non-reversed genes: ", n_inconsistent,
                      "\nNumber of reversed genes: ", n_consistent)
  }
  
  pt_comp_bar <- pt_comp_df_wide %>%
    filter(!is.na(effect) & effect != "Inconsistent") %>%  # Filter out inconsistent genes
    ggplot(aes(x = effect, fill = effect)) +
    geom_bar_rounded() +
    geom_text(stat  = "count",
              aes(label = after_stat(count)),
              vjust = 1.5, hjust = 0.5, color = "white") +
    theme_bw() +
    theme(panel.grid  = element_blank(),
          panel.border = element_blank(),
          text         = element_text(size = 15),
          legend.position = "none",
          axis.text.x  = element_text(angle = 50, hjust = 1),
          axis.ticks   = element_blank(),
          plot.caption = element_text(hjust = 0.5, size = 12)) +
    labs(y = "Count", x = NULL,
         caption = caption) +
    scale_fill_manual(values = effect_cols)
  
  ## 4. Label colouring for reversed genes ------------------------------------
  reversed_cols <- effect_cols[c(
    "Reversed towards +\n(T1D depleted)",
    "Reversed towards -\n(T1D elevated)"
  )]
  
  pt_gene_color_df <- pt_comp_df_wide %>%
    filter(str_detect(effect, "Reversed")) %>%
    dplyr::mutate(Gene_label = ifelse(
      str_detect(effect, "depleted"),
      paste0("<span style='color:", reversed_cols[1], "'>", Gene, "</span>"),
      paste0("<span style='color:", reversed_cols[2], "'>", Gene, "</span>")
    )) %>%
    dplyr::select(Gene, Gene_label)
  
  ## 5. Dot plot ---------------------------------------------------------------
  if (nrow(pt_comp_df_wide)>0) {
    pt_comp_dot <- pt_comp_df %>%
      left_join(pt_gene_color_df, by = "Gene") %>%
      filter(Gene %in% pt_gene_color_df$Gene) %>%
      dplyr::mutate(group = ifelse(source == "ATTEMPT",
                                   "Dapa vs. Placebo", "T1D vs. HC"),
                    Gene_label = dplyr::coalesce(Gene_label, Gene)) %>%
      ggplot(aes(x = group, y = Gene_label,
                 color = direction, size = abs(FC))) +
      geom_point() +
      theme_bw() +
      scale_color_manual(values = c("+" = "#f28482", "-" = "#457b9d")) +
      theme(panel.border = element_blank(),
            panel.grid   = element_blank(),
            legend.position = "top",
            legend.box      = "vertical",
            axis.text.x  = element_markdown(angle = 50, hjust = 1),
            axis.text.y  = element_markdown(angle = 0, hjust = 1),
            axis.ticks   = element_blank(),
            text = element_text(size = 15)) +
      scale_y_discrete(position = "right") +
      labs(x = NULL, y = NULL, color = "Direction",
           caption = caption)
    
    ## 6. Combine & optionally save ---------------------------------------------
    combined <-pt_comp_dot # + pt_comp_bar
    
    if (!is.null(save_path)) {
      ggsave(save_path, combined, width = width, height = height)
    }
    
    return(combined)
  }
  
}

# ===========================================================================
# Function: make_plot_df
# ===========================================================================

make_plot_df <- function(res_df, celltype_label, pos_genes, neg_genes) {
  res_df %>%
    dplyr::select(Gene, 
                  logFC = `logFC_treatmentDapagliflozin:visitPOST`, 
                  pval = `p_treatmentDapagliflozin:visitPOST`, 
                  fdr = fdr) %>%
    dplyr::mutate(celltype = celltype_label,
                  direction = case_when(
                    logFC > 0 ~ "+",
                    logFC < 0 ~ "-"
                  )) %>%
    filter(Gene %in% c(pos_genes, neg_genes))
}


# ===========================================================================
# Function: plot_volcano_proteomics (for proteomics)
# ===========================================================================
plot_volcano_proteomics <- function(data, fc, p_col, title = NULL, x_axis, y_axis, file_suffix, p_thresh = 0.05,
                                    positive_text = "Positive with Dapagliflozin", 
                                    negative_text = "Negative with Dapagliflozin",
                                    formula = "group", legend_position = c(0.8, 0.9),
                                    text_size = 15, caption_size = 8.5,
                                    top_n = 20,
                                    output_base_path = file.path(root_path, "ATTEMPT/Results/Figures/Proteomics/Volcano Plots/Volcano"),
                                    geom_text_size = 3,
                                    arrow_padding = 0.09,
                                    arrow_text_padding = 0.14,
                                    legend_text_size = 9,
                                    caption_padding = 15,
                                    x_title_padding_t = 32) {
  set.seed(1)
  top_pos <- data %>%
    dplyr::filter(!!sym(fc) > 0 & !!sym(p_col) < p_thresh) %>%
    dplyr::arrange(!!sym(p_col))
  
  n_pos <- nrow(top_pos)
  
  top_pos_n <- top_pos %>%
    slice_head(n=top_n)
  
  top_neg <- data %>%
    dplyr::filter(!!sym(fc) < 0 & !!sym(p_col) < p_thresh) %>%
    dplyr::arrange(!!sym(p_col))
  
  n_neg <- nrow(top_neg)
  
  top_neg_n <- top_neg %>%
    slice_head(n=top_n)
  
  data <- data %>%
    dplyr::mutate(top_color = case_when(data$AptName %in% top_pos$AptName ~ "#f28482",
                                        data$AptName %in% top_neg$AptName ~ "#457b9d",
                                        TRUE ~ "#ced4da"),
                  top_size = if_else(data$AptName %in% c(top_pos$AptName, top_neg$AptName), 1.3, 1),
                  EntrezGeneSymbol = if_else(data$EntrezGeneSymbol == "", data$Target, data$EntrezGeneSymbol),
                  top_lab  = if_else(data$AptName %in% c(top_pos_n$AptName, top_neg_n$AptName), EntrezGeneSymbol, ""))
  
  # Max and min for annotation arrows
  max_fc <- max(data[[fc]], na.rm = TRUE)
  min_fc <- min(data[[fc]], na.rm = TRUE)
  
  # Get y-axis max for dynamic scaling
  y_max <- max(-log10(data[[p_col]]), na.rm = TRUE) * 1.1
  
  
  p <- ggplot(data, aes(x = !!sym(fc), y = -log10(!!sym(p_col)))) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "darkgrey") +
    geom_point(alpha = 0.5, aes(color = top_color, size = top_size)) +
    geom_label_repel(aes(label = top_lab, color = top_color),
                    size = geom_text_size, max.overlaps = Inf,
                    force = 6, segment.alpha = 0.3, segment.size = 0.3,
                    fontface = "bold",
                    fill = fill_alpha("white", 0.7),
                    label.size = 0) +
    labs(title = paste(title),
         x = paste(x_axis),
         y = paste(y_axis),
         caption = paste0("Formula: log2(protein) ~ ", formula, 
                          "\nTotal analyzed proteins n: ", nrow(data),
                          " | Positive n = ", n_pos, " | Negative n = ", n_neg)) +
    scale_size_continuous(range = c(1, 1.3)) + 
    scale_color_manual(values = c("#457b9d"="#457b9d", "#ced4da"="#ced4da", "#f28482"="#f28482")) +
    theme_minimal() +
    guides(color = "none", size = "none")  +
    annotate("segment", 
             x=max_fc/8, 
             xend=(max_fc*7)/8, 
             y=-y_max * arrow_padding,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(max_fc/8, (max_fc*7)/8)), 
             y=-y_max * arrow_text_padding, 
             label=positive_text,
             size=geom_text_size, color="#343a40") +
    annotate("segment", 
             x=min_fc/8, 
             xend=(min_fc*7)/8, 
             y=-y_max * arrow_padding,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(min_fc/8, (min_fc*7)/8)), 
             y=-y_max * arrow_text_padding, 
             label=negative_text,
             size=geom_text_size, color="#343a40") +
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(ylim = c(0, y_max), clip="off") +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = text_size, family = "Arial"),
          title = element_text(size = legend_text_size, family = "Arial"),
          legend.position = legend_position,
          legend.justification = if(is.character(legend_position) && legend_position == "bottom") c(0.5, 0) else c(1, 1),
          legend.direction = "horizontal",
          legend.spacing.x = unit(0.05, 'cm'),
          plot.margin = margin(t = 10, r = 20, b = 25, l = 20),
          axis.title.x = element_text(margin = margin(t = x_title_padding_t)),
          plot.caption = element_text(size = caption_size, hjust = 0.5, margin = margin(t = caption_padding)),
          legend.margin = margin(t = 5, b = 5))
  ggsave(paste0(output_base_path, file_suffix, ".png"), plot = p, width = 7, height = 5)
  
  return(p)
}


# ===========================================================================
# Function: create_directional_volcano (for transcripts)
# ===========================================================================

create_directional_volcano <- function(rna_data, protein_data, cell_type_name, 
                                       rna_logfc_col = "logFC_treatmentDapagliflozin:visitPOST",
                                       rna_p_col = "p_treatmentDapagliflozin:visitPOST",
                                       protein_logfc_col = "logFC_log2(protein)", 
                                       protein_p_col = "p_log2(protein)",
                                       p_threshold = 0.05,
                                       top_n = 20) {
  # Add epsilon for log transformation
  epsilon <- 1e-300
  
  # Calculate -log10(p) with epsilon
  protein_data <- protein_data %>%
    mutate(neg_log_p = -log10(!!sym(protein_p_col) + epsilon))
  
  theme_transparent <- theme(
    plot.background   = element_rect(fill = "transparent", color = NA),
    panel.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
  # Create gene direction from RNA data
  gene_dir <- rna_data %>%
    mutate(gene_dir = case_when(!!sym(rna_logfc_col) < 0 & !!sym(rna_p_col) < p_threshold ~ "-",
                                !!sym(rna_logfc_col) > 0 & !!sym(rna_p_col) < p_threshold ~ "+")) %>%
    dplyr::select(Gene, gene_dir) %>%
    filter(!is.na(gene_dir))
  
  # Create protein direction from protein data
  protein_dir <- protein_data %>%
    mutate(protein_dir = case_when(!!sym(protein_logfc_col) < 0 & !!sym(protein_p_col) < p_threshold ~ "-",
                                   !!sym(protein_logfc_col) > 0 & !!sym(protein_p_col) < p_threshold ~ "+")) %>%
    dplyr::select(Gene, protein_dir) %>%
    filter(!is.na(protein_dir))
  
  # Merge and find matches
  dir_merged <- full_join(gene_dir, protein_dir, by = "Gene")
  
  dir_match <- dir_merged %>%
    mutate(match = case_when(protein_dir == gene_dir ~ "match", 
                             TRUE ~ "non-match")) %>%
    filter(match == "match")
  
  dir_match_genes <- dir_match$Gene
  
  # Calculate positive and negative counts
  pos_n <- dir_match %>%
    filter(protein_dir == "+") %>%
    distinct(Gene) %>%
    nrow()
  
  neg_n <- dir_match %>%
    filter(protein_dir == "-") %>%
    distinct(Gene) %>%
    nrow()
  
  # Prepare data for plotting
  plot_data <- protein_data %>%
    filter(Gene %in% dir_merged$Gene) %>%
    group_by(Gene) %>%
    slice_min(!!sym(protein_p_col), n = 1) %>%
    ungroup()
  
  # Create label data - top N from each direction
  label_data <- plot_data %>%
    filter(Gene %in% dir_match_genes) %>%
    mutate(direction = case_when(
      !!sym(protein_logfc_col) < 0 ~ "down",
      !!sym(protein_logfc_col) > 0 ~ "up"
    )) %>%
    group_by(direction) %>%
    slice_min(!!sym(protein_p_col), n = top_n) %>%
    ungroup()
  
  # Create the plot
  plot <- plot_data %>%
    ggplot(aes(x = !!sym(protein_logfc_col), y = neg_log_p)) +
    geom_hline(yintercept = -log10(p_threshold), color = "darkgrey", linetype = "dashed") +
    geom_point(aes(color = case_when(
      Gene %in% dir_match_genes & !!sym(protein_logfc_col) < 0 ~ "match_down",
      Gene %in% dir_match_genes & !!sym(protein_logfc_col) > 0 ~ "match_up",
      TRUE ~ "non-match"
    )), alpha = 0.3, size = 2) +
    geom_text_repel(data = label_data,
                    aes(label = Gene, color = case_when(
                      !!sym(protein_logfc_col) < 0 ~ "match_down",
                      !!sym(protein_logfc_col) > 0 ~ "match_up",
                      TRUE ~ "non-match"
                    )), alpha = 0.8,
                    max.overlaps = Inf) +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = 15, family = "Arial"),
          legend.position = "none",
          plot.caption = element_text(size = 8.5, hjust = 0.5, margin = margin(t = 15)),
          panel.background = element_blank()) + 
    scale_color_manual(values = c("match_down" = "#457b9d",
                                  "match_up" = "#f28482",
                                  "non-match" = "#e5e5e5")) +
    labs(x = "logFC log2(protein)", 
         y = "-log10(p.value)",
         caption = paste0("Formula: log2(protein) + treatment + (1|subject)\nCell type: ", cell_type_name, 
                          " | Positive n = ", pos_n, " | Negative n = ", neg_n)) +
    theme_transparent
  
  return(plot)
}

# ===========================================================================
# Function: plot_volcano (for transcripts)
# ===========================================================================
plot_volcano <- function(data, fc, p_col, title = NULL, x_axis, y_axis, file_suffix, p_thresh = 0.05,
                         positive_text = "Positive with Dapagliflozin", 
                         negative_text = "Negative with Dapagliflozin",
                         formula = "group", legend_position = c(0.8, 0.9),
                         full_formula = F,
                         text_size = 15, caption_size = 8.5,
                         cell_type = "",
                         output_base_path = file.path(root_path, "ATTEMPT/Results/Figures/Volcano Plots/"),
                         geom_text_size = 3,
                         arrow_padding = 0.09,
                         arrow_text_padding = 0.14,
                         legend_text_size = 9,
                         caption_padding = 15,
                         x_title_padding_t = 32,
                         genes_to_label= NULL,
                         volcano_force = 6,
                         volcano_box_padding = 0,
                         off_chart_threshold = 0.95,  # New parameter: proportion of y_max to consider "off chart"
                         off_chart_y_position = 0.85,  # New parameter: where to place off-chart labels
                         off_chart_arrow_length = 0.02) {  # New parameter: length of off-chart arrows
  
  set.seed(1)
  
  # Add epsilon for log transformation
  epsilon <- 1e-300
  
  # Calculate -log10(p) with epsilon
  data <- data %>%
    mutate(neg_log_p = -log10(!!sym(p_col) + epsilon))
  
  # Get y-axis max for dynamic scaling
  y_max <- max(data$neg_log_p, na.rm = TRUE) * 1.1
  y_cutoff <- y_max * off_chart_threshold  # Points above this are "off chart"
  
  top_pos <- data %>%
    dplyr::filter(!!sym(fc) > 0 & !!sym(p_col) < p_thresh) %>%
    dplyr::arrange(!!sym(p_col))
  
  n_pos <- nrow(top_pos)
  
  top_pos_n <- top_pos %>%
    filter(if (!is.null(genes_to_label)) Gene %in% genes_to_label else TRUE) %>%
    slice_head(n=20)
  
  top_neg <- data %>%
    dplyr::filter(!!sym(fc) < 0 & !!sym(p_col) < p_thresh) %>%
    dplyr::arrange(!!sym(p_col))
  
  n_neg <- nrow(top_neg)
  
  top_neg_n <- top_neg %>%
    filter(if (!is.null(genes_to_label)) Gene %in% genes_to_label else TRUE) %>%
    slice_head(n=20)
  
  # Identify off-chart genes
  off_chart_genes <- data %>%
    filter(Gene %in% c(top_pos_n$Gene, top_neg_n$Gene) & neg_log_p > y_cutoff) %>%
    mutate(
      is_positive = !!sym(fc) > 0,
      # Spread out x positions for off-chart labels
      x_position = if_else(is_positive,
                           !!sym(fc) + seq(from = 0.1, by = 0.2, length.out = n()),
                           !!sym(fc) - seq(from = 0.1, by = 0.2, length.out = n())),
      y_position = y_max * off_chart_y_position
    )
  
  # Separate on-chart and off-chart genes for labeling
  on_chart_genes <- c(top_pos_n$Gene, top_neg_n$Gene)[!c(top_pos_n$Gene, top_neg_n$Gene) %in% off_chart_genes$Gene]
  
  data <- data %>%
    dplyr::mutate(
      top_color = case_when(
        Gene %in% top_pos$Gene ~ "#f28482",
        Gene %in% top_neg$Gene ~ "#457b9d",
        TRUE ~ "#ced4da"
      ),
      top_size = if_else(Gene %in% c(top_pos$Gene, top_neg$Gene), 1.3, 1),
      # Only label on-chart genes normally
      top_lab  = if_else(Gene %in% on_chart_genes, Gene, ""),
      # Cap display values at y_cutoff for plotting
      display_neg_log_p = pmin(neg_log_p, y_cutoff)
    ) %>%
    filter(abs(!!sym(fc)) < 10)
  
  # Max and min for annotation arrows
  max_fc <- max(data[[fc]], na.rm = TRUE)
  min_fc <- min(data[[fc]], na.rm = TRUE)
  
  p <- ggplot(data, aes(x = !!sym(fc), y = display_neg_log_p)) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "darkgrey") +
    geom_point(alpha = 0.5, aes(color = top_color), size = 3) +
    # Regular labels for on-chart genes
    geom_label_repel(seed = 1, 
      data = filter(data, top_lab != ""),
      aes(label = top_lab, color = top_color),
      fontface = "bold",
      size = geom_text_size, max.overlaps = Inf,
      force = volcano_force, segment.alpha = 0.3, segment.size = 0.3,
      box.padding = volcano_box_padding,
      fill = fill_alpha("white", 0.7),     # background color of the label box
      label.size = 0
      # bg.color = "black",
      # bg.r = .15
    ) +
    # Add arrows for off-chart genes
    {if(nrow(off_chart_genes) > 0) {
      list(
        geom_segment(
          data = off_chart_genes,
          aes(x = !!sym(fc), y = y_cutoff * 0.98,
              xend = !!sym(fc), yend = y_cutoff - (y_max * off_chart_arrow_length)),
          arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
          color = "black",
          size = 0.6
        ),
        geom_text(
          data = off_chart_genes,
          aes(x = x_position, y = y_position, label = Gene),
          size = geom_text_size * 0.9,
          hjust = if_else(off_chart_genes$is_positive, 0, 1),
          color = "black",
          fontface = "italic"
        ),
        # Add p-value annotation for off-chart genes
        geom_text(
          data = off_chart_genes,
          aes(x = x_position, y = y_position - (y_max * 0.03), 
              label = paste0("p=", format(!!sym(p_col), scientific = TRUE, digits = 2))),
          size = geom_text_size * 0.7,
          hjust = if_else(off_chart_genes$is_positive, 0, 1),
          color = "darkgrey"
        )
      )
    }} +
    labs(title = paste(title),
         x = paste(x_axis),
         y = paste(y_axis),
         caption = if (!is.null(cell_type) & !full_formula & !is.null(formula)) {
           paste0("Formula: ~ ", formula, " + (1|subject)", 
                  "\n\nCell type: ", 
                  cell_type, if(cell_type != "") " | " else "", 
                  "Positive n = ", n_pos, " | Negative n = ", n_neg,
                  if(nrow(off_chart_genes) > 0) paste0("\n", nrow(off_chart_genes), " gene(s) with p-values near zero (arrows indicate off-scale values)") else "")
         } else if(!is.null(cell_type) & full_formula) {
           bquote(
             atop(
               Expression[ij] == (beta[0] + b[0*i]) +
                 beta[1]*visit[ij] +
                 beta[2]*treatment[i] +
                 beta[3]*(visit[ij]*treatment[i]) +
                 epsilon[ij],
               "Cell type:" ~ .(cell_type) ~ "|" ~ 
                 "Positive n =" ~ .(n_pos) ~ "|" ~ 
                 "Negative n =" ~ .(n_neg) ~
                 .(if(nrow(off_chart_genes) > 0) paste0(" | ", nrow(off_chart_genes), " off-scale") else "")
             )
           )
         } else {
           paste0("\nCell type: ", 
                  cell_type, if(cell_type != "") " | " else "", 
                  "Positive n = ", n_pos, " | Negative n = ", n_neg,
                  if(nrow(off_chart_genes) > 0) paste0("\n", nrow(off_chart_genes), " gene(s) with p-values near zero (arrows indicate off-scale values)") else "")
         }) +
    scale_size_continuous(range = c(1, 1.3)) + 
    scale_color_manual(values = c("#457b9d"="#457b9d", "#ced4da"="#ced4da", "#f28482"="#f28482")) +
    # theme_minimal() +
    guides(color = "none", size = "none")  +
    annotate("segment", 
             x=max_fc/8, 
             xend=(max_fc*7)/8, 
             y=-y_max * arrow_padding,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(max_fc/8, (max_fc*7)/8)), 
             y=-y_max * arrow_text_padding, 
             label=positive_text,
             size=geom_text_size, color="#343a40") +
    annotate("segment", 
             x=min_fc/8, 
             xend=(min_fc*7)/8, 
             y=-y_max * arrow_padding,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(min_fc/8, (min_fc*7)/8)), 
             y=-y_max * arrow_text_padding, 
             label=negative_text,
             size=geom_text_size, color="#343a40") +
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(ylim = c(0, y_max), clip="off") +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = text_size, family = "Arial"),
          title = element_text(size = legend_text_size, family = "Arial"),
          legend.position = legend_position,
          legend.justification = if(is.character(legend_position) && legend_position == "bottom") c(0.5, 0) else c(1, 1),
          legend.direction = "horizontal",
          legend.spacing.x = unit(0.05, 'cm'),
          plot.margin = margin(t = 10, r = 20, b = 25, l = 20),
          axis.title.x = element_text(margin = margin(t = x_title_padding_t)),
          plot.caption = element_text(size = caption_size, hjust = 0.5, margin = margin(t = caption_padding)),
          legend.margin = margin(t = 5, b = 5),
          panel.background = element_blank(),
          legend.background = element_blank())
  
  if (!is.null(output_base_path)) {
    ggsave(paste0(output_base_path, file_suffix, ".png"), bg = "transparent", plot = p, width = 7, height = 5)
  }
  return(p)
}

# ===========================================================================
# Function: plot_volcano_associations
# ===========================================================================
plot_volcano_associations <- function(clin_results, fc, p_col, title_suffix, 
                                      x_axis, y_axis, file_suffix, p_thresh = 0.05,  treatment_results, treatment_p_col, 
                                      positive_text = "Positive with Dapagliflozin", 
                                      negative_text = "Negative with Dapagliflozin",
                                      formula = "\u0394 treatment", color_by = "treatment_direction",
                                      legend_position = c(0.3,0.8)) {
  
  set.seed(1)
  
  # Add treatment effect directions
  treatment_results <- treatment_results %>%
    dplyr::mutate(treatment_direction = case_when(`p_treatmentDapagliflozin:visitPOST` < 0.05 &
                                                    `logFC_treatmentDapagliflozin:visitPOST` > 0 ~ "Up w/ Dapagliflozin",
                                                  `p_treatmentDapagliflozin:visitPOST` < 0.05 &
                                                    `logFC_treatmentDapagliflozin:visitPOST` < 0 ~ "Down w/ Dapagliflozin",
                                                  TRUE ~ "NS"))
  
  # Genes significant in both models
  sig_clin_genes <- clin_results %>%
    filter(!!sym(p_col) < p_thresh) %>%
    pull(Gene)
  
  sig_treat_genes <- treatment_results %>%
    filter(!!sym(treatment_p_col) < p_thresh) %>%
    pull(Gene)
  
  sig_both_genes <- intersect(sig_clin_genes, sig_treat_genes)
  
  # Top 5 left pos / neg and right pos / neg
  top_5_left_pos <- clin_results %>%
    filter(Gene %in% sig_both_genes, !!sym(fc) < 0) %>%
    left_join(treatment_results %>% 
                dplyr::select(Gene, `logFC_treatmentDapagliflozin:visitPOST`), by = "Gene") %>%
    filter(`logFC_treatmentDapagliflozin:visitPOST` > 0) %>%
    dplyr::arrange(!!sym(p_col)) %>%
    slice_head(n = 5) %>%
    pull(Gene)
  
  top_5_left_neg <- clin_results %>%
    filter(Gene %in% sig_both_genes, !!sym(fc) < 0) %>%
    left_join(treatment_results %>% 
                dplyr::select(Gene, `logFC_treatmentDapagliflozin:visitPOST`), by = "Gene") %>%
    filter(`logFC_treatmentDapagliflozin:visitPOST` < 0) %>%
    dplyr::arrange(!!sym(p_col)) %>%
    slice_head(n = 5) %>%
    pull(Gene)
  
  top_5_right_pos <- clin_results %>%
    filter(Gene %in% sig_both_genes, !!sym(fc) > 0) %>%
    left_join(treatment_results %>% 
                dplyr::select(Gene, `logFC_treatmentDapagliflozin:visitPOST`), by = "Gene") %>%
    filter(`logFC_treatmentDapagliflozin:visitPOST` > 0) %>%
    dplyr::arrange(!!sym(p_col)) %>%
    slice_head(n = 5) %>%
    pull(Gene)
  
  top_5_right_neg <- clin_results %>%
    filter(Gene %in% sig_both_genes, !!sym(fc) > 0) %>%
    left_join(treatment_results %>% 
                dplyr::select(Gene, `logFC_treatmentDapagliflozin:visitPOST`), by = "Gene") %>%
    filter(`logFC_treatmentDapagliflozin:visitPOST` < 0) %>%
    dplyr::arrange(!!sym(p_col)) %>%
    slice_head(n = 5) %>%
    pull(Gene)
  
  # Combine for plotting
  top_label_genes <- c(top_5_left_pos, top_5_left_neg, top_5_right_pos, top_5_right_neg)
  
  message("Number of genes significant in both: ", length(sig_both_genes))
  message("Number of labeled genes: ", length(top_label_genes))
  
  clin_results <- clin_results %>%
    dplyr::mutate(shape_var_plot = case_when(
      Gene %in% sig_both_genes & !!sym(fc) < 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Up w/ Dapagliflozin" ~ "Left_Pos",
      Gene %in% sig_both_genes & !!sym(fc) < 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Down w/ Dapagliflozin" ~ "Left_Neg",
      Gene %in% sig_both_genes & !!sym(fc) > 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Up w/ Dapagliflozin" ~ "Right_Pos",
      Gene %in% sig_both_genes & !!sym(fc) > 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Down w/ Dapagliflozin" ~ "Right_Neg",
      TRUE ~ "NS"
    ),
    color_var_plot = case_when(
      Gene %in% sig_both_genes ~ treatment_results$treatment_direction[match(Gene, treatment_results$Gene)],
      TRUE ~ "NS"
    ),
    top_lab = if_else(Gene %in% top_label_genes, Gene, "")
    )
  
  clin_results$shape_var_plot <- factor(clin_results$shape_var_plot, 
                                        levels = c("Left_Pos", "Left_Neg", "Right_Pos", "Right_Neg", "NS"))
  
  clin_results$color_var_plot <- factor(clin_results$color_var_plot, 
                                        levels = c("Up w/ Dapagliflozin", "Down w/ Dapagliflozin", "NS"))
  
  # Max and min for annotation arrows
  max_fc <- max(clin_results[[fc]], na.rm = TRUE)
  min_fc <- min(clin_results[[fc]], na.rm = TRUE)
  
  # Get y-axis max for dynamic scaling
  y_max <- max(-log10(clin_results[[p_col]]), na.rm = TRUE) * 1.1
  
  # Plot
  p <- ggplot(clin_results, aes(x = !!sym(fc), y = -log10(!!sym(p_col)))) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "darkgrey") +
    geom_point(alpha = 0.3, aes(fill = color_var_plot, 
                                color = color_var_plot,
                                shape = shape_var_plot), size = 2) +
    scale_shape_manual(values = c("Left_Pos" = 23, "Left_Neg" = 23, "Right_Pos" = 22, "Right_Neg" = 22, "NS" = 21),
                       labels = c("Left_Pos" = "Left Pos", "Left_Neg" = "Left Neg", 
                                  "Right_Pos" = "Right Pos", "Right_Neg" = "Right Neg", "NS" = "NS")) +
    scale_color_manual(values = c("Up w/ Dapagliflozin" = "#f28482", "Down w/ Dapagliflozin" = "#457b9d", "NS" = "#ced4da")) +
    scale_fill_manual(values = c("Up w/ Dapagliflozin" = "#f28482", "Down w/ Dapagliflozin" = "#457b9d", "NS" = "#ced4da")) +
    guides(color = guide_legend(title = NULL), fill = "none") +
    geom_text_repel(aes(label = top_lab, color = color_var_plot),
                    size = 3, max.overlaps = Inf, force = 10,
                    segment.alpha = 0.5, segment.size = 0.4,
                    min.segment.length = 0, box.padding = 0.6, point.padding = 0.4,
                    segment.color = "#ced4da") +
    labs(title = paste(title_suffix),
         x = paste(x_axis),
         y = paste(y_axis)) +
    theme_minimal() +
    annotate("segment", 
             x=max_fc/8, 
             xend=(max_fc*7)/8, 
             y=-y_max * 0.09,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(max_fc/8, (max_fc*7)/8)), 
             y=-y_max * 0.14, 
             label=positive_text,
             size=3, color="#343a40") +
    annotate("segment", 
             x=min_fc/8, 
             xend=(min_fc*7)/8, 
             y=-y_max * 0.09,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(min_fc/8, (min_fc*7)/8)), 
             y=-y_max * 0.14, 
             label=negative_text,
             size=3, color="#343a40") +
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(ylim = c(0, y_max), clip="off") +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = 9, family = "Arial"),
          legend.position = legend_position,
          legend.justification = if(is.character(legend_position) && legend_position == "bottom") c(0.5, 0) else c(1, 1),
          legend.direction = "horizontal",
          legend.spacing.x = unit(0.3, 'cm'),
          plot.margin = margin(t = 10, r = 20, b = 25, l = 20),
          axis.title.x = element_text(margin = margin(t = 28)),
          plot.caption = element_text(size = 8.5, hjust = 0.5, margin = margin(t = 8), family = "Arial"),
          legend.margin = margin(t = 5, b = 5)) + 
    guides(shape = "none")
  
  # Save
  ggsave(paste0(file.path(root_path, "ATTEMPT/Results/Figures/Volcano Plots/"), file_suffix, ".png"), plot = p, width = 7, height = 5)
  
  return(p)
}

# ===========================================================================
# Function: plot_volcano_concordance
# ===========================================================================

plot_volcano_concordance <- function(clin_results, fc, p_col, 
                                     x_axis, y_axis, file_suffix, p_thresh = 0.05,  treatment_results, treatment_p_col, 
                                     positive_text = "Positive with Dapagliflozin", 
                                     negative_text = "Negative with Dapagliflozin",
                                     arrow_label = "",
                                     formula = "\u0394 treatment", 
                                     color_by = "treatment_direction",
                                     legend_position = "top",
                                     clinical_direction = NULL,
                                     caption_text = paste0("Point position reflects association with clinical variable; point color indicates treatment effect direction.\nPoints are colored if concordant with clinical variable direction after treatment. \nUp to top ", top_n, " from each direction are labeled."),
                                     cell_type = "", top_n = 50,
                                     seed = 1, 
                                     force = 50,
                                     force_pull = 1,
                                     box.padding = 0.6,
                                     point.padding = 0.4, 
                                     text_size = 3,
                                     title_size = 9,
                                     arrow_padding = 0.09,
                                     arrow_text_padding = 0.18,
                                     x_axis_padding = 38,
                                     # New size parameters
                                     label_size = 3,           # Size for gene labels (default same as text_size)
                                     arrow_symbol_size = 6,    # Size for up/down arrow symbols
                                     main_arrow_size = 10,     # Size for the central arrow
                                     axis_title_size = 9,      # Size for axis titles (default same as title_size)
                                     axis_text_size = 9,       # Size for axis tick labels
                                     legend_text_size = 9,     # Size for legend text
                                     caption_size = 9,         # Size for caption text
                                     annotation_text_size = 3, # Size for annotation text (positive/negative text)
                                     nudge_y = 0.3,
                                     nudge_x = 0.2,
                                     y_expand = 0,
                                     ylim = 0.2,
                                     save = T) {
  
  set.seed(seed)
  # Add treatment effect directions
  treatment_results <- treatment_results %>%
    dplyr::mutate(treatment_direction = case_when(`p_treatmentDapagliflozin:visitPOST` < 0.05 &
                                                    `logFC_treatmentDapagliflozin:visitPOST` > 0 ~ "Up w/ Dapagliflozin",
                                                  `p_treatmentDapagliflozin:visitPOST` < 0.05 &
                                                    `logFC_treatmentDapagliflozin:visitPOST` < 0 ~ "Down w/ Dapagliflozin",
                                                  TRUE ~ "NS/NC"))
  
  # Calculate concordant genes based on clinical direction
  concordant_genes <- c()
  if (!is.null(clinical_direction)) {
    if (clinical_direction == "-") {
      # For negative clinical direction: clinical - and treatment +
      concordant_genes <- clin_results %>%
        left_join(treatment_results %>% 
                    dplyr::select(Gene, `logFC_treatmentDapagliflozin:visitPOST`, `p_treatmentDapagliflozin:visitPOST`), by = "Gene") %>%
        filter((!!sym(fc) > 0 & !!sym(p_col) < 0.05) & (`logFC_treatmentDapagliflozin:visitPOST` < 0 & `p_treatmentDapagliflozin:visitPOST` < 0.05)|
                 (!!sym(fc) < 0 & !!sym(p_col) < 0.05) & (`logFC_treatmentDapagliflozin:visitPOST` > 0 & `p_treatmentDapagliflozin:visitPOST` < 0.05)) %>%
        pull(Gene)
    } else if (clinical_direction == "+") {
      # For positive clinical direction: clinical + and treatment -
      concordant_genes <- clin_results %>%
        left_join(treatment_results %>% 
                    dplyr::select(Gene, `logFC_treatmentDapagliflozin:visitPOST`, `p_treatmentDapagliflozin:visitPOST`), by = "Gene") %>%
        filter((!!sym(fc) > 0 & !!sym(p_col) < 0.05) & (`logFC_treatmentDapagliflozin:visitPOST` > 0 & `p_treatmentDapagliflozin:visitPOST` < 0.05)|
                 (!!sym(fc) < 0 & !!sym(p_col) < 0.05) & (`logFC_treatmentDapagliflozin:visitPOST` < 0 & `p_treatmentDapagliflozin:visitPOST` < 0.05)) %>%
        pull(Gene)
    }
  }
  
  n_concordant <- length(concordant_genes)
  sig_clin_genes <- clin_results %>%
    filter(!!sym(p_col) < p_thresh) %>%
    pull(Gene)
  
  sig_treat_genes <- treatment_results %>%
    filter(!!sym(treatment_p_col) < p_thresh) %>%
    pull(Gene)
  
  sig_both_genes <- intersect(sig_clin_genes, sig_treat_genes)
  
  # Label top 50 concordant genes on each side
  top_50_left <- clin_results %>%
    filter(Gene %in% concordant_genes, !!sym(fc) < 0) %>%
    dplyr::arrange(!!sym(p_col)) %>%
    slice_head(n = top_n) %>%
    pull(Gene)
  
  top_50_right <- clin_results %>%
    filter(Gene %in% concordant_genes, !!sym(fc) > 0) %>%
    dplyr::arrange(!!sym(p_col)) %>%
    slice_head(n = top_n) %>%
    pull(Gene)
  
  # Combine for plotting
  top_label_genes <- c(top_50_left, top_50_right)
  
  message("Number of genes significant in both: ", length(sig_both_genes))
  message("Number of concordant genes: ", n_concordant)
  message("Number of labeled genes: ", length(top_label_genes))
  
  clin_results <- clin_results %>%
    dplyr::mutate(
      color_var_plot = case_when(
        Gene %in% concordant_genes ~ treatment_results$treatment_direction[match(Gene, treatment_results$Gene)],
        TRUE ~ "NS/NC"
      ),
      shape_var_plot = case_when(
        Gene %in% sig_both_genes & !!sym(fc) < 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Up w/ Dapagliflozin" ~ "Left_Pos",
        Gene %in% sig_both_genes & !!sym(fc) < 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Down w/ Dapagliflozin" ~ "Left_Neg",
        Gene %in% sig_both_genes & !!sym(fc) > 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Up w/ Dapagliflozin" ~ "Right_Pos",
        Gene %in% sig_both_genes & !!sym(fc) > 0 & treatment_results$treatment_direction[match(Gene, treatment_results$Gene)] == "Down w/ Dapagliflozin" ~ "Right_Neg",
        TRUE ~ "NS/NC"
      ),
      top_lab = if_else(Gene %in% top_label_genes, Gene, ""),
      # Updated arrow logic based on clinical direction - only for concordant genes
      arrow_symbol = case_when(
        Gene %in% concordant_genes & Gene %in% top_label_genes & !is.null(clinical_direction) & clinical_direction == "+" & !!sym(fc) < 0 ~ "\u2193",
        Gene %in% concordant_genes & Gene %in% top_label_genes & !is.null(clinical_direction) & clinical_direction == "+" & !!sym(fc) > 0 ~ "\u2191",
        Gene %in% concordant_genes & Gene %in% top_label_genes & !is.null(clinical_direction) & clinical_direction == "-" & !!sym(fc) < 0 ~ "\u2191",
        Gene %in% concordant_genes & Gene %in% top_label_genes & !is.null(clinical_direction) & clinical_direction == "-" & !!sym(fc) > 0 ~ "\u2193",
        TRUE ~ ""
      )
    )
  
  clin_results$shape_var_plot <- factor(clin_results$shape_var_plot, 
                                        levels = c("Left_Pos", "Left_Neg", "Right_Pos", "Right_Neg", "NS/NC"))
  
  clin_results$color_var_plot <- factor(clin_results$color_var_plot, 
                                        levels = c("Up w/ Dapagliflozin", "Down w/ Dapagliflozin", "NS/NC"))
  
  # Max and min for annotation arrows
  max_fc <- max(clin_results[[fc]], na.rm = TRUE)
  min_fc <- min(clin_results[[fc]], na.rm = TRUE)
  
  # Get y-axis max for dynamic scaling
  y_max <- max(-log10(clin_results[[p_col]]), na.rm = TRUE) * 1.1
  
  # Plot (no background shading)
  p <- ggplot(clin_results, aes(x = !!sym(fc), y = -log10(!!sym(p_col)))) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "darkgrey") +
    geom_point(alpha = 0.4, aes(fill = color_var_plot, 
                                color = color_var_plot), size = 3) +
    geom_text(aes(label = arrow_symbol, color = color_var_plot,
                  vjust = ifelse(arrow_symbol == "\u2193", 1, -0.5)),
              size = arrow_symbol_size, family = "Arial", alpha = 0.3) +
    scale_shape_manual(values = c("Left_Pos" = 22, "Left_Neg" = 22, "Right_Pos" = 22, "Right_Neg" = 22, "NS/NC" = 22),
                       labels = c("Left_Pos" = "Left Pos", "Left_Neg" = "Left Neg", 
                                  "Right_Pos" = "Right Pos", "Right_Neg" = "Right Neg", "NS/NC" = "NS/NC")) +
    scale_color_manual(values = c("Up w/ Dapagliflozin" = "#f28482", "Down w/ Dapagliflozin" = "#457b9d", "NS/NC" = "#ced4da")) +
    scale_fill_manual(values = c("Up w/ Dapagliflozin" = "#f28482", "Down w/ Dapagliflozin" = "#457b9d", "NS/NC" = "#ced4da")) +
    guides(color = guide_legend(title = NULL, 
                                override.aes = list(label = ""), 
                                nrow = 1), 
           fill = "none",
           shape = "none") +
    annotate("segment",
             x = 0, xend = 0,
             y = if(!is.null(clinical_direction) && clinical_direction == "+") y_max * 0.4 else y_max * 0.9, 
             yend = if(!is.null(clinical_direction) && clinical_direction == "+") y_max * 0.9 else y_max * 0.4,
             size = main_arrow_size, linejoin = "mitre",
             color = if(!is.null(clinical_direction) && clinical_direction == "+") "#f7c1bf" else "#a9c9dd", 
             arrow = arrow(type = "closed")) +
    geom_vline(xintercept = 0, color = "darkgrey", linetype = "dashed") +
    annotate("text", x = 0, y = ((y_max*0.9 + y_max*0.4)/2),
             label = arrow_label,
             fontface = "bold",
             color = if(!is.null(clinical_direction) && clinical_direction == "+") "#f28482" else "#457b9d",
             size = text_size) +
    
    # Split geom_text_repel for positive logFC (right side)
    geom_label_repel(data = clin_results[clin_results[[fc]] > 0 & clin_results$top_lab != "", ],
                    aes(label = top_lab, color = color_var_plot),
                    nudge_x = max_fc * nudge_x,
                    nudge_y = nudge_y,
                    xlim = c(max_fc * 0.15, NA),
                    ylim = c(y_max*ylim, NA),
                    segment.size = 0.4,
                    segment.color = "darkgrey",
                    segment.alpha = 0.3,
                    hjust = 0,
                    size = label_size,
                    max.overlaps = Inf,
                    force = force,
                    force_pull = force_pull,
                    box.padding = box.padding,
                    point.padding = point.padding,
                    min.segment.length = 0,
                    seed = seed,
                    fontface = "bold",
                    fill = fill_alpha("white", 0.7),
                    label.size = 0) +
    
    # Split geom_text_repel for negative logFC (left side)
    geom_label_repel(data = clin_results[clin_results[[fc]] < 0 & clin_results$top_lab != "", ],
                    aes(label = top_lab, color = color_var_plot),
                    nudge_x = min_fc * nudge_x,
                    nudge_y = nudge_y,
                    xlim = c(NA, (min_fc * 0.15)),
                    ylim = c(y_max*ylim, NA),
                    segment.size = 0.4,
                    segment.color = "darkgrey",
                    segment.alpha = 0.3,
                    hjust = 1,
                    size = label_size,
                    max.overlaps = Inf,
                    force = force,
                    force_pull = force_pull,
                    box.padding = box.padding,
                    point.padding = point.padding,
                    min.segment.length = 0,
                    seed = seed,
                    fontface = "bold",
                    fill = fill_alpha("white", 0.7),
                    label.size = 0) +
    
    labs(title = NULL,
         x = paste(x_axis),
         y = paste(y_axis),
         caption = if (!is.null(clinical_direction) & !is.null(formula)) {
           paste0("Formula: ~ ", formula, " + (1|subject)", 
                  "\n\nCell type: ", 
                  cell_type, if(cell_type != "") " | " else "", 
                  "Concordant genes: n = ", n_concordant,
                  "\n\n", caption_text)
         } else if (!is.null(clinical_direction) & is.null(formula) & is.null(caption_text)) {
           paste0("Cell type: ", 
                  cell_type, if(cell_type != "") " | " else "", 
                  "Concordant genes: n = ", n_concordant)
         } else NULL) +
    theme_minimal() +
    annotate("segment", 
             x=max_fc/10, 
             xend=(max_fc*9)/10, 
             y=-y_max * arrow_padding,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(max_fc/10, (max_fc*9)/10)), 
             y=-y_max * arrow_text_padding, 
             label=positive_text,
             size=annotation_text_size, color="#343a40") +
    annotate("segment", 
             x=min_fc/10, 
             xend=(min_fc*9)/10, 
             y=-y_max * arrow_padding,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(min_fc/10, (min_fc*9)/10)), 
             y=-y_max * arrow_text_padding, 
             label=negative_text,
             size=annotation_text_size, color="#343a40") +
    scale_y_continuous(expand=c(0,y_expand)) +
    scale_x_continuous(expand = expansion(mult = c(0.25, 0.25))) + 
    coord_cartesian(ylim = c(0, y_max), clip="off") +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = title_size, family = "Arial"),
          legend.position = legend_position,
          legend.justification = if(is.character(legend_position) && legend_position == "bottom") c(0.5, 0) else c(1, 1),
          legend.direction = "horizontal",
          legend.spacing.x = unit(0.3, 'cm'),
          plot.margin = margin(t = 10, r = 20, b = 10, l = 20),
          axis.title.x = element_text(margin = margin(t = x_axis_padding), size = axis_title_size),
          axis.title.y = element_text(size = axis_title_size),
          axis.text.x = element_text(size = axis_text_size),
          axis.text.y = element_text(size = axis_text_size),
          legend.text = element_text(size = legend_text_size),
          plot.caption = element_text(size = caption_size, hjust = 0.5, margin = margin(t = 8), family = "Arial"),
          legend.margin = margin(t = 5, b = 5)) + 
    guides(shape = "none")
  
  # Save
  if (save) {
    ggsave(paste0(file.path(root_path, "ATTEMPT/Results/Figures/Volcano Plots/"), file_suffix, "_concordance.png"), plot = p, width = 7, height = 5)
    
    return(p)
    
  } else return(p)
}
# ===========================================================================
# Function: plot_mean_ci_stars
# ===========================================================================
plot_mean_ci_stars <- function(data, y_var, y_axis_title, 
                               legend_position = c(0.6, 0.8), 
                               baseline_visit = 0, 
                               visits_to_plot = c(-4, 0, 4, 16, 18),
                               covariates = NULL,
                               test_method = "lmer_covars") { # test method can be "lmer", "ancova", "lmer_covars"
  
  dodge_val <- 0.08
  y_sym <- rlang::ensym(y_var)
  y_name <- rlang::as_name(y_sym)
  baseline_var <- paste0(y_name, "_bl")
  
  # Add baseline values
  baseline_df <- data %>%
    filter(visit == baseline_visit) %>%
    dplyr::select(subject_id, !!y_sym) %>%
    rename_with(~ baseline_var, all_of(y_name))
  
  data <- data %>%
    filter(visit %in% visits_to_plot) %>%
    left_join(baseline_df, by = "subject_id")
  
  # --- Mean ± 95% CI ---
  mean_dat <- data %>%
    group_by(treatment, visit) %>%
    dplyr::summarise(
      n = sum(!is.na(!!y_sym)),
      mean_y = mean(!!y_sym, na.rm = TRUE),
      sd_y = sd(!!y_sym, na.rm = TRUE),
      t_crit = qt(0.975, df = n - 1),
      lower = mean_y - (t_crit * sd_y / sqrt(n)),
      upper = mean_y + (t_crit * sd_y / sqrt(n)),
      .groups = "drop"
    )
  
  # --- Choose significance test ---
  contrast_results <- purrr::map_dfr(visits_to_plot, function(v) {
    tryCatch({
      if (test_method == "ancova") {
        ancova <- run_ancova(data, outcome_var = y_name, baseline_var = baseline_var, visit_week = v)
        ancova$contrasts %>%
          filter(contrast == "Placebo - Dapagliflozin 5mg") %>%
          dplyr::mutate(visit = v)
        
      } else if (test_method == "lmer_covars") {
        model <- lmer(
          formula(paste0(y_name, " ~ treatment * factor(visit) + age + sex + diabetes_dx_duration + bmi + (1|subject_id)")),
          data = data
        )
        emm <- emmeans(model, ~ treatment | visit)
        contrast(emm, method = "revpairwise", by = "visit", adjust = "bonferroni") %>%
          as_tibble() %>% 
          filter(visit == v) %>% 
          dplyr::mutate(contrast = contrast, estimate = estimate, p.value = p.value)
        
      } else if (test_method == "lmer_covars_site") {
        model <- lmer(
          formula(paste0(y_name, " ~ treatment * factor(visit) + age + sex + diabetes_dx_duration + bmi + site + (1|subject_id)")),
          data = data
        )
        emm <- emmeans(model, ~ treatment | visit)
        contrast(emm, method = "revpairwise", by = "visit", adjust = "bonferroni") %>%
          as_tibble() %>% 
          filter(visit == v) %>% 
          dplyr::mutate(contrast = contrast, estimate = estimate, p.value = p.value)
        
      } else if (test_method == "lmer") {
        model <- lmer(
          formula(paste0(y_name, " ~ treatment * factor(visit) + (1|subject_id)")),
          data = data
        )
        emm <- emmeans(model, ~ treatment | visit)
        contrast(emm, method = "revpairwise", by = "visit", adjust = "bonferroni") %>%
          as_tibble() %>% 
          filter(visit == v) %>% 
          dplyr::mutate(contrast = contrast, estimate = estimate, p.value = p.value)
      }
    }, error = function(e) {
      message(glue::glue("Significance test failed for visit {v}: {e$message}"))
      tibble(visit = v, p.value = NA_real_, estimate = NA_real_,
             conf.low = NA, conf.high = NA, contrast = NA)
    })
  }) %>%
    dplyr::mutate(
      stars = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    )
  
  # --- Y position for stars ---
  ttest_results <- contrast_results %>%
    left_join(
      mean_dat %>% 
        group_by(visit) %>%
        dplyr::summarise(max_upper = max(upper, na.rm = TRUE), .groups = "drop"),
      by = "visit"
    ) %>%
    dplyr::mutate(y_position = max_upper + 0.05 * max_upper)
  
  print(contrast_results)
  # --- Plot ---
  p <- data %>%
    ggplot(aes(x = factor(visit, levels = visits_to_plot), y = !!y_sym, color = treatment)) +
    geom_errorbar(data = mean_dat, 
                  aes(x = factor(visit, levels = visits_to_plot), ymin = lower, ymax = upper, group = treatment), 
                  inherit.aes = FALSE,
                  width = 0.1, size = 0.3, position = position_dodge(width = dodge_val),
                  color = "black") +
    geom_line(data = mean_dat, aes(y = mean_y, group = treatment), 
              position = position_dodge(width = dodge_val),
              size = 1.2) +
    geom_point(data = mean_dat, aes(y = mean_y, shape = treatment), 
               size = 11, position = position_dodge(width = dodge_val), color = "white") +
    geom_point(data = mean_dat, aes(y = mean_y, shape = treatment), 
               size = 7, position = position_dodge(width = dodge_val)) +
    geom_text(data = ttest_results, 
              aes(x = factor(visit, levels = visits_to_plot), y = y_position, label = stars), 
              inherit.aes = FALSE,
              size = 6, vjust = 0) +
    scale_shape_manual(values = c("Placebo" = 15, "Dapagliflozin 5mg" = 18)) +
    scale_color_manual(values = c("Placebo" = "#f8ae9d", "Dapagliflozin 5mg" = "#a7b298")) +
    scale_x_discrete(expand = expansion(mult = c(0.3, 0.3))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_minimal(base_size = 15) +
    theme(panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = legend_position,
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = "transparent", color = NA),
          panel.grid = element_blank()) +
    labs(x = "Visit (weeks)", 
         y = y_axis_title,
         color = NULL, shape = NULL) +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  return(p)
}

# ===========================================================================
# Function: create_gene_expression_plots
# ===========================================================================

create_gene_expression_plots <- function(main_results,
                                         subtype_results_list,
                                         cell_type_labels = NULL,
                                         cell_type_order = cell_type_labels,
                                         cell_type_prefix = "PT",
                                         volcano_text_size = 15,
                                         volcano_caption_size = 15,
                                         dot_text_size = 15,
                                         n_top_genes = 20,
                                         output_dir = ".",
                                         save_plots = TRUE,
                                         output_prefix = cell_type_prefix,
                                         logfc_col = "logFC_treatmentDapagliflozin:visitPOST",
                                         fdr_col = "fdr",
                                         volcano_pval_col = "fdr",
                                         formula = "",
                                         full_formula = F,
                                         fdr_threshold = 0.05,
                                         heatmap_caption = T, 
                                         arrow_padding = 0.05,
                                         arrow_text_padding = 0.08,
                                         caption_padding = 10, 
                                         x_title_padding_t = 0,
                                         geom_text_size = 5,
                                         volcano_force = 6,
                                         volcano_box_padding = 0) { 
  
  # Load required libraries
  require(tidyverse)
  require(ggtext)
  require(patchwork)
  theme_transparent <- theme(
    plot.background   = element_rect(fill = "transparent", color = NA),
    panel.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
  
  main_results <- main_results %>%
    filter(abs(!!sym(logfc_col)) < 10)
  
  subtype_results_list <- lapply(subtype_results_list, function(df) {
    df %>% filter(abs(!!sym(logfc_col)) < 10)
  })
  
  # Determine cell type labels and mapping (keep existing code)
  if (is.null(cell_type_labels)) {
    cell_type_labels <- paste0(cell_type_prefix, "-", names(subtype_results_list))
  } else {
    if (length(cell_type_labels) != length(subtype_results_list)) {
      stop("Length of cell_type_labels must match length of subtype_results_list")
    }
  }
  
  label_mapping <- setNames(cell_type_labels, names(subtype_results_list))
  
  if (is.null(cell_type_order)) {
    cell_type_order <- rev(cell_type_labels)
  } else {
    if (!all(cell_type_labels %in% cell_type_order) || !all(cell_type_order %in% cell_type_labels)) {
      stop("cell_type_order must contain exactly the same labels as cell_type_labels")
    }
  }
  
  if (is.null(output_prefix)) {
    output_prefix <- ifelse(is.null(cell_type_labels), cell_type_prefix, "celltypes")
  }
  
  # MODIFIED: Select top positive genes (up to 20, but only if significant)
  top_pos_genes <- main_results %>%
    filter(!!sym(fdr_col) < fdr_threshold) %>%  # Only significant genes
    filter(!!sym(logfc_col) > 0) %>%
    dplyr::arrange(!!sym(fdr_col)) %>%
    head(n_top_genes) %>%
    pull(Gene)
  
  # MODIFIED: Select top negative genes (up to 20, but only if significant)
  top_neg_genes <- main_results %>%
    filter(!!sym(fdr_col) < fdr_threshold) %>%  # Only significant genes
    filter(!!sym(logfc_col) < 0) %>%
    dplyr::arrange(!!sym(fdr_col)) %>%
    head(n_top_genes) %>%
    pull(Gene)
  
  # All selected genes will be labeled in volcano
  volcano_genes <- c(top_neg_genes, top_pos_genes)
  
  # Create combined plot data
  combined_plot_dat <- purrr::imap_dfr(subtype_results_list, function(res_df, subtype_key) {
    label <- label_mapping[subtype_key]
    make_plot_df(res_df, label, top_pos_genes, top_neg_genes)
  })
  
  # Set factor levels for ordering
  combined_plot_dat$celltype <- factor(combined_plot_dat$celltype, levels = cell_type_order)
  combined_plot_dat$Gene <- factor(combined_plot_dat$Gene, levels = c(top_neg_genes, top_pos_genes))
  
  # Handle missing genes
  if (length(unique(combined_plot_dat$Gene)) < length(volcano_genes)) {
    missing_genes <- volcano_genes[!(volcano_genes %in% combined_plot_dat$Gene)]
    if (length(missing_genes) > 0) {
      missing_df <- data.frame(
        Gene = rep(missing_genes, each = length(cell_type_order)),
        celltype = rep(cell_type_order, length(missing_genes)),
        direction = "+"
      )
      combined_plot_dat <- dplyr::bind_rows(combined_plot_dat, missing_df)
    }
  }
  
  # Create color vector for x-axis labels
  gene_levels <- c(top_neg_genes, top_pos_genes)
  x_colors <- setNames(
    c(rep("#457b9d", length(top_neg_genes)), 
      rep("#f28482", length(top_pos_genes))),
    gene_levels
  )
  
  # Create dot plot
  dot_plot <- create_dot_plot(combined_plot_dat, x_colors, text_size = dot_text_size)  + theme_transparent
  
  # Get main logFC values for the selected genes
  main_logfc_values <- main_results %>%
    filter(Gene %in% volcano_genes) %>%
    dplyr::select(Gene, main_logFC = !!sym(logfc_col))
  
  # Create consistency score data with weighted scores
  consistency_scores <- calculate_weighted_consistency_scores(
    combined_plot_dat, 
    main_logfc_values,
    top_neg_genes, 
    top_pos_genes, 
    cell_type_labels
  )
  
  # Create bar plot
  bar_plot <- create_bar_plot(consistency_scores) + theme_transparent
  
  # Create volcano plot
  volcano_plot <- NULL
  if (!is.null(volcano_pval_col)) {
    volcano_plot <- plot_volcano(
      main_results, 
      logfc_col, 
      volcano_pval_col,
      title = NULL,
      paste0("logFC ", gsub("logFC_", "", logfc_col)), 
      "-log10(FDR adjusted p-value)",
      formula = formula,
      text_size = volcano_text_size,
      caption_size = volcano_caption_size,
      full_formula = full_formula,
      output_base_path = NULL,
      cell_type = cell_type_prefix,
      geom_text_size = geom_text_size,
      arrow_padding = arrow_padding,
      arrow_text_padding = arrow_text_padding,
      legend_text_size = volcano_text_size,
      caption_padding = caption_padding,
      x_title_padding_t = x_title_padding_t,
      genes_to_label = volcano_genes,  # All selected significant genes
      volcano_force = volcano_force,
      volcano_box_padding = volcano_box_padding
    )  + theme_transparent
  }
  
  # Create heatmap with only significant genes
  heatmap_plot <- plot_treatment_heatmap(
    data = main_results,
    heatmap_caption = heatmap_caption,
    genes_to_show = volcano_genes  # Only show significant genes
  )
  
  # Combine plots (rest of the code remains the same)
  if (!is.null(volcano_plot)) {
    design <- "
ACD
BCD"
    combined_plot <- bar_plot + dot_plot + volcano_plot + heatmap_plot + 
      plot_layout(design = design,
                  widths  = c(1, 1,  0.1),
                  heights = c(0.25, 1, 0.5)) &
      theme(plot.background = element_rect(fill = "transparent", color = NA))
    
    if (save_plots) {
      output_path <- file.path(output_dir, paste0(output_prefix, "_subtypes_nebula_scores_volcano_heat.png"))
      ggsave(output_path, width = 20, height = 10, plot = combined_plot, bg = "transparent")
    }
  }
  
  # Save individual plots if requested
  if (save_plots) {
    ggsave(file.path(output_dir, paste0(output_prefix, "_subtypes_nebula.png")), 
           width = 7, height = 5, plot = dot_plot)
  }
  
  # Return all components
  return(list(
    dot_plot = dot_plot,
    bar_plot = bar_plot,
    volcano_plot = volcano_plot,
    heatmap_plot = heatmap_plot,
    combined_plot = combined_plot,
    top_pos_genes = top_pos_genes,
    top_neg_genes = top_neg_genes,
    consistency_scores = consistency_scores,
    combined_plot_data = combined_plot_dat,
    volcano_genes = volcano_genes  # Also return the genes that were labeled
  ))
}
# ===========================================================================
# Function: create_dot_plot
# ===========================================================================

#' Create dot plot for gene expression
#' Create dot plot for gene expression
create_dot_plot <- function(plot_data, x_colors, text_size = 10) {
  plot_data$Gene <- factor(plot_data$Gene, levels = (names(x_colors)))
  
  # Add this line to count negative genes
  n_neg_genes <- sum(x_colors == "#457b9d")
  
  plot_data %>%
    ggplot(aes(x = Gene, y = celltype, size = abs(logFC), color = direction)) +
    annotate("rect", 
             xmin = -Inf, 
             xmax = 0.5 + n_neg_genes,  # Changed this line
             ymin = -Inf, 
             ymax = Inf, 
             fill = "#dceef5", 
             alpha = 0.5) +  
    annotate("rect", 
             xmin = 0.5 + n_neg_genes,  # Changed this line
             xmax = Inf, 
             ymin = -Inf, 
             ymax = Inf, 
             fill = "#f9dada", 
             alpha = 0.5) +
    geom_point() +
    scale_color_manual(values = c("+" = "#f28482", "-" = "#457b9d")) +
    scale_x_discrete(
      guide = guide_axis(angle = 70), 
      labels = function(x) {
        mapply(function(label, col) {
          glue::glue("<span style='color:{col}'>{label}</span>")
        }, x, x_colors[x], SIMPLIFY = TRUE)
      }
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_markdown(angle = 70, hjust = 1),
      panel.border = element_blank(),
      text = element_text(size = text_size),
      legend.key = element_rect(fill = NA),
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.position = c(0.5, .95),
      legend.background = element_rect(fill = NA, color = NA),
      legend.title = element_text(size = 8),
      plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
    ) +
    labs(x = NULL, y = NULL, color = "Direction")
}
# ===========================================================================
# Function: calculate_weighted_consistency_scores
# ===========================================================================

#' Calculate weighted consistency scores for genes based on logFC differences
calculate_weighted_consistency_scores <- function(plot_data, main_logfc_values, neg_genes, pos_genes, cell_type_labels) {
  # Join with main logFC values
  plot_data_with_main <- plot_data %>%
    left_join(main_logfc_values, by = "Gene")
  
  # Calculate weighted scores and mismatch flag
  weighted_scores <- plot_data_with_main %>%
    group_by(Gene) %>%
    dplyr::summarize(
      # Get main logFC
      main_logFC = dplyr::first(main_logFC),
      # Track mismatch
      has_mismatch = any(sign(logFC) != sign(main_logFC), na.rm = TRUE),
      score = {
        if (abs(main_logFC) < 0.001) {
          0
        } else {
          raw_weights <- logFC / main_logFC
          penalty_factor <- 2
          weights <- ifelse(
            sign(logFC) == sign(main_logFC),
            raw_weights,
            raw_weights * penalty_factor
          )
          weights <- pmax(pmin(weights, 2), -2)
          mean_weight <- mean(weights, na.rm = TRUE)
          prop_opposite <- sum(sign(logFC) != sign(main_logFC), na.rm = TRUE) / length(logFC)
          mean_weight <- mean_weight * (1 - 0.5 * prop_opposite)
          mean_weight
        }
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      direction = case_when(main_logFC < 0 ~ "-", main_logFC > 0 ~ "+"),
      bar_color = ifelse(has_mismatch, "gray", direction),
      Gene = factor(Gene, levels = c(neg_genes, pos_genes))
    ) %>%
    dplyr::arrange(Gene)
  
  return(weighted_scores)
}

# ===========================================================================
# Function: create_bar_plot
# ===========================================================================

# Create bar plot for consistency scores
create_bar_plot <- function(scores_data) {
  ggplot(scores_data, aes(x = Gene, y = score, fill = bar_color)) +
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.4, color = "#a7c957") +
    geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.4, color = "#f4a261") +
    geom_col_rounded(width = 0.8) +
    geom_text(
      aes(label = Gene),
      vjust = 0.5, 
      hjust = ifelse(scores_data$score >= 0, 1.3, -0.3),
      size = 2,
      angle = 90,
      color = "white"
    ) +
    scale_fill_manual(values = c("+" = "#f28482", "-" = "#457b9d", "gray" = "gray60")) +
    scale_y_continuous(
      breaks = c(-2, -1, 0, 1, 2),
      labels = c("-2", "-1", "0", "1", "2"),
      limits = c(min(-0.1, min(scores_data$score, na.rm = TRUE) * 1.1), 
                 max(0.1, max(scores_data$score, na.rm = TRUE) * 1.1))
    ) +
    theme_bw() + 
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      panel.border = element_blank(),
      text = element_text(size = 10),
      axis.text.y = element_text(size = 6),
      axis.title.y = element_text(size = 8),
      legend.position = "none",
      legend.background = element_rect(fill = NA, color = NA),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
    ) +
    labs(x = NULL, y = "Weighted \nconsistency score")
}

# ---- Pathway Plot Functions ----

# ===========================================================================
# Function: matrix_to_list
# ===========================================================================
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}
# ===========================================================================
# Function: prepare_gmt
# ===========================================================================
library(fgsea)
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}

# List of acronyms to uppercase
acronyms <- c("rna", "dna", "mtor", "foxo", "ppar", "nmd", "fgfr", "robo", 
              "bhl", "cov", "jak", "stat", "wnt", "hiv", "bcl", "mapk",
              "pt", "tal", "pc", "ic", "ec", "fibvsmcp")
special_mixed <- c("rrna", "mrna", "trna", "gtpase", "atpase", "robos", "slits", "fibvsmcp")
special_replacements <- c("rRNA", "mRNA", "tRNA", "GTPase", "ATPase", "ROBOs", "SLITs", "FIB/VSMC/P")


# ===========================================================================
# Function: replace_mixed_case
# ===========================================================================

replace_mixed_case <- function(text, from, to) {
  for (i in seq_along(from)) {
    pattern <- paste0("\\b", from[i], "\\b")
    text <- str_replace_all(text, regex(pattern, ignore_case = TRUE), to[i])
  }
  return(text)
}


# ===========================================================================
# Function: capitalize_acronyms
# ===========================================================================

capitalize_acronyms <- function(text, terms) {
  for (term in terms) {
    pattern <- paste0("\\b", term, "\\b")
    replacement <- toupper(term)
    text <- str_replace_all(text, regex(pattern, ignore_case = TRUE), replacement)
  }
  return(text)
}

# ===========================================================================
# Function: clean_pathway_names
# ===========================================================================

clean_pathway_names <- function(pathways) {
  # Remove REACTOME_ prefix
  cleaned <- gsub("^REACTOME_", "", pathways)
  cleaned <- gsub("^GOBP_", "", cleaned)
  cleaned <- gsub("^KEGG_", "", cleaned)
  cleaned <- gsub("^HALLMARK_", "", cleaned)
  
  # Replace underscores with spaces
  cleaned <- gsub("_", " ", cleaned)
  
  # Convert to title case (capitalize first letter of each word)
  cleaned <- tools::toTitleCase(tolower(cleaned))
  
  # Define patterns that should be in uppercase
  uppercase_words <- c(
    # Roman numerals
    "\\bI\\b", "\\bIi\\b", "\\bIii\\b", "\\bIv\\b", "\\bV\\b", "\\bVi\\b", 
    "\\bVii\\b", "\\bViii\\b", "\\bIx\\b", "\\bX\\b",
    # Common biological abbreviations
    "\\bTca\\b", "\\bAtp\\b", "\\bAdp\\b", "\\bAmp\\b", "\\bGtp\\b", "\\bGdp\\b",
    "\\bNad\\b", "\\bNadh\\b", "\\bFad\\b", "\\bFadh2\\b", "\\bCoa\\b",
    "\\bDna\\b", "\\bRna\\b", "\\bMrna\\b", "\\bTrna\\b", "\\bRrna\\b",
    "\\bEr\\b", "\\bUpr\\b", "\\bNf\\b", "\\bHif\\b", "\\bMhc\\b",
    "\\bTgf\\b", "\\bEgf\\b", "\\bVegf\\b", "\\bPdgf\\b", "\\bFgf\\b",
    "\\bRos\\b", "\\bRns\\b", "\\bNo\\b", "\\bNos\\b", "\\bInos\\b", "\\bEnos\\b", "\\bNnos\\b",
    "\\Mapk\\b"
  )
  
  # Apply uppercase replacements
  for (pattern in uppercase_words) {
    # Extract the word without the word boundaries
    word <- gsub("\\\\b", "", pattern)
    # Convert to uppercase
    replacement <- toupper(word)
    # Replace in the cleaned string
    cleaned <- gsub(pattern, replacement, cleaned, ignore.case = TRUE)
  }
  
  # Fix specific cases that might need custom handling
  cleaned <- gsub("\\bOf\\b", "of", cleaned)  # lowercase 'of'
  cleaned <- gsub("\\bBy\\b", "by", cleaned)  # lowercase 'by'
  cleaned <- gsub("\\bThe\\b", "the", cleaned)  # lowercase 'the' (except at start)
  cleaned <- gsub("^the", "The", cleaned)  # capitalize 'the' at start
  cleaned <- gsub("\\bO Linked\\b", "O-Linked", cleaned, ignore.case = TRUE)
  
  return(cleaned)
}

# ===========================================================================
# Function: filter_redundant_pathways
# ===========================================================================

filter_redundant_pathways <- function(gsea_result, overlap_pct = 0.3) {
  # Check required columns
  if (!all(c("pathway", "leadingEdge", "padj") %in% colnames(gsea_result))) {
    stop("Input data frame must have 'pathway', 'leadingEdge', and 'padj' columns.")
  }
  
  # Extract leading edge and name them
  leading_edges <- gsea_result$leadingEdge
  names(leading_edges) <- gsea_result$pathway
  
  # Compute Jaccard overlap matrix
  overlap_matrix <- sapply(leading_edges, function(x) 
    sapply(leading_edges, function(y) 
      length(intersect(x, y)) / length(union(x, y))
    )
  )
  
  # Identify pairs with high overlap
  overlap_pairs <- which(overlap_matrix > overlap_pct & lower.tri(overlap_matrix), arr.ind = TRUE)
  
  # If no overlaps found, return original
  if (nrow(overlap_pairs) == 0) {
    message("No redundant pathways found with overlap above threshold.")
    return(gsea_result)
  }
  
  # Create overlap graph
  library(igraph)
  edges <- data.frame(
    from = rownames(overlap_matrix)[overlap_pairs[,1]],
    to   = colnames(overlap_matrix)[overlap_pairs[,2]]
  )
  g <- graph_from_data_frame(edges, directed = FALSE)
  components <- components(g)
  
  # Cluster terms and select representative from each
  redundant_clusters <- split(names(components$membership), components$membership)
  representative_terms <- sapply(redundant_clusters, function(cluster_terms) {
    sub <- gsea_result[gsea_result$pathway %in% cluster_terms, ]
    sub[which.min(sub$padj), "pathway"]
  })
  
  # Keep only representative + non-overlapping terms
  all_overlapping <- unique(unlist(redundant_clusters))
  non_overlapping <- setdiff(gsea_result$pathway, all_overlapping)
  final_terms <- c(non_overlapping, representative_terms)
  
  # Return filtered GSEA result
  gsea_result[gsea_result$pathway %in% final_terms, ]
}

# ===========================================================================
# Function: plot_fgsea_transpose
# ===========================================================================

plot_fgsea_transpose <- function(fgsea_res,
                                 top_n = 30,
                                 title = "Top Enriched Pathways",
                                 xmin = 0,
                                 xmax = 3,
                                 xnudge = (xmax - xmin)/100,
                                 text1 = 6.5,
                                 text2 = 18,
                                 text3 = 20) {
  
  fgsea_res <- fgsea_res %>%
    dplyr::arrange(pval) %>%
    head(top_n) %>%
    dplyr::mutate(
      direction = case_when((NES < 0 & pval <= 0.05 ~ "Negative"), 
                            (NES > 0 & pval <= 0.05 ~ "Positive"),
                            (NES < 0 & pval > 0.05 ~ "Negative p > 0.05"), 
                            (NES > 0 & pval > 0.05 ~ "Positive p > 0.05")),
      face = case_when((NES < 0 & pval <= 0.05 ~ "bold"), 
                       (NES > 0 & pval <= 0.05 ~ "bold"),
                       (NES < 0 & pval > 0.05 ~ "plain"), 
                       (NES > 0 & pval > 0.05 ~ "plain")),
      pathway_clean = str_remove(pathway, "^KEGG_"), 
      pathway_clean = str_remove(pathway_clean, "^REACTOME_"), 
      pathway_clean = str_remove(pathway_clean, "^GOBP_"), 
      pathway_clean = str_remove(pathway_clean, "^GOMF_"), 
      pathway_clean = str_replace_all(pathway_clean, "_", " "),
      pathway_clean = str_to_sentence(pathway_clean),
      pathway_clean = str_replace_all(pathway_clean, "\\bi\\b", "I"),
      pathway_clean = str_replace_all(pathway_clean, "\\bii\\b", "II"),
      pathway_clean = str_replace_all(pathway_clean, "\\biii\\b", "III"),
      pathway_clean = str_replace_all(pathway_clean, "\\biv\\b", "IV"),
      pathway_clean = str_replace_all(pathway_clean, "\\bv\\b", "V"),
      pathway_clean = str_replace_all(pathway_clean, regex("\\(immune\\)", ignore_case = TRUE), "(IMMUNE)"),
      pathway_clean = capitalize_acronyms(pathway_clean, acronyms),
      pathway_clean = replace_mixed_case(pathway_clean, special_mixed, special_replacements),
      pathway_clean = paste0(pathway_clean, " (", size, ")")
    ) %>%
    dplyr::arrange(pval)
  
  fgsea_res$pathway_clean <- reorder(fgsea_res$pathway_clean, -abs(fgsea_res$NES))
  
  fgsea_res %>%
    ggplot(aes(x = abs(NES), y = fct_rev(pathway_clean), label = pathway_clean)) +
    geom_point(aes(size = -log10(pval), color = direction, alpha = 0.8)) +
    # geom_vline(xintercept = 2, linetype = "dashed") +
    geom_text(aes(group = pathway_clean, color = direction, fontface = face), 
              hjust = 0, size = text1, nudge_x = xnudge) +
    scale_size_binned() +
    scale_color_manual(values = c("Positive" = "#c75146", "Negative" = "#2c7da0", 
                                  "Positive p > 0.05" = "#e18c80", "Negative p > 0.05" = "#7ab6d1")) +
    scale_x_continuous(limits = c(xmin, xmax), expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(expand = expansion(add = 1)) +
    labs(
      x = "NES",
      y = "Pathways",
      color = "Direction",
      size = "-log(p-value)",
      title = title
    ) +
    guides(alpha = "none") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = text3),
      axis.title = element_text(size = text3),
      axis.ticks.y = element_blank(), 
      legend.position = c(0.9, 0.2),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(color = "black"),
      legend.title = element_text(size = text2),
      legend.text = element_text(size = text2),
      title = element_text(size = text3)
    )
}

# ---- Summary Functions ----

# ===========================================================================
# Function: trt_run_cell_type_analysis
# ===========================================================================

# Function to run complete analysis for a given cell type
trt_run_cell_type_analysis <- function(cell_type, 
                                       input_path, 
                                       input_suffix,
                                       output_base_path,
                                       output_prefix,
                                       bg_path = file.path(root_path, "GSEA/"),
                                       bucket = "attempt",
                                       region = "") {
  
  # Print status
  cat("Starting analysis for cell type:", cell_type, "\n")
  
  # Create cell type specific paths
  cell_type_lower <- tolower(cell_type)
  
  # Read in nebula results
  input_file <- file.path(input_path, cell_type, "nebula", 
                          paste0(cell_type_lower, input_suffix))
  
  res <- s3readRDS(input_file, bucket = bucket, region = region)
  
  # Check convergence
  res_convergence <- map_dfr(
    names(res),
    function(gene_name) {
      converged <- res[[gene_name]]$convergence
      df <- data.frame(Gene = gene_name,
                       Convergence_Code = converged)
      return(df)
    }
  )
  
  res_converged <- res_convergence %>%
    filter(Convergence_Code >= -10)
  
  # Combine results
  res_combined <- map_dfr(
    names(res),
    function(gene_name) {
      df <- res[[gene_name]]$summary
      df <- df %>% 
        filter(gene_name %in% res_converged$Gene) %>%
        dplyr::mutate(Gene = gene_name)
      return(df)
    }
  )
  
  # Calculate FDR
  res_combined <- res_combined %>%
    ungroup() %>%
    dplyr::mutate(fdr = p.adjust(`p_treatmentDapagliflozin:visitPOST`, method = "fdr"))
  
  res_combined <- subset(res_combined, abs(`logFC_treatmentDapagliflozin:visitPOST`) < 10)
  
  # Save results
  output_csv <- file.path(output_base_path, "Results", "NEBULA", 
                          paste0(output_prefix, cell_type_lower, "_nebula_res.csv"))
  write.csv(res_combined, output_csv, row.names = FALSE)
  
}
# ===========================================================================
# Function: t1dhc_run_cell_type_analysis
# ===========================================================================

# Function to run complete analysis for a given cell type
t1dhc_run_cell_type_analysis <- function(cell_type, 
                                         input_path, 
                                         input_suffix,
                                         output_base_path,
                                         output_prefix,
                                         bg_path = file.path(root_path, "GSEA/"),
                                         plot_title = "Volcano Plot",
                                         bucket = "attempt",
                                         region = "") {
  
  # Print status
  cat("Starting analysis for cell type:", cell_type, "\n")
  
  # Create cell type specific paths
  cell_type_lower <- tolower(cell_type)
  
  # Read in nebula results
  input_file <- file.path(input_path, cell_type, "nebula", 
                          paste0(cell_type_lower, input_suffix))
  
  res <- s3readRDS(input_file, bucket = bucket, region = region)
  
  # Check convergence
  res_convergence <- map_dfr(
    names(res),
    function(gene_name) {
      converged <- res[[gene_name]]$convergence
      df <- data.frame(Gene = gene_name,
                       Convergence_Code = converged)
      return(df)
    }
  )
  
  res_converged <- res_convergence %>%
    filter(Convergence_Code >= -10)
  
  # Combine results
  res_combined <- map_dfr(
    names(res),
    function(gene_name) {
      df <- res[[gene_name]]$summary
      df <- df %>% 
        filter(gene_name %in% res_converged$Gene) %>%
        dplyr::mutate(Gene = gene_name)
      return(df)
    }
  )
  
  # Calculate FDR
  res_combined <- res_combined %>%
    ungroup() %>%
    dplyr::mutate(fdr = p.adjust(p_groupType_1_Diabetes, method = "fdr"))
  
  res_combined <- subset(res_combined, abs(logFC_groupType_1_Diabetes) < 10)
  # Save results
  
  output_csv <- file.path(output_base_path, "Results", "NEBULA", 
                          paste0(output_prefix, cell_type_lower, "_nebula_res.csv"))
  write.csv(res_combined, output_csv, row.names = FALSE)
  
  # Add direction columns for visualization
  res_combined <- res_combined %>%
    dplyr::mutate(group_direction = case_when(
      p_groupType_1_Diabetes < 0.05 & logFC_groupType_1_Diabetes > 0 ~ "Positive",
      p_groupType_1_Diabetes < 0.05 & logFC_groupType_1_Diabetes < 0 ~ "Negative",
      TRUE ~ "NS"
    ),
    group_direction_fdr = case_when(
      fdr < 0.05 & logFC_groupType_1_Diabetes > 0 ~ "Positive",
      fdr < 0.05 & logFC_groupType_1_Diabetes < 0 ~ "Negative",
      TRUE ~ "NS"
    ))
  
  # Create volcano plots
  plot_volcano(res_combined, 
               "logFC_groupType_1_Diabetes", 
               "p_groupType_1_Diabetes",
               NULL,
               "logFC group", 
               "-log10(p-value)",
               "_deg",
               formula = "group",
               positive_text = "Upregulated in T1D compared to LC",
               negative_text = "Downregulated in T1D compared to LC",
               cell_type = cell_type,
               output_base_path = file.path(output_base_path, "Results", "Figures", "Volcano", "nebula",
                                            paste0(output_prefix, cell_type_lower)))
  
  plot_volcano(res_combined, 
               "logFC_groupType_1_Diabetes", 
               "fdr",
               NULL,
               "logFC group", 
               "-log10(FDR adjusted p-value)",
               "_deg_fdr",
               formula = "group",
               positive_text = "Upregulated in T1D compared to LC",
               negative_text = "Downregulated in T1D compared to LC",
               cell_type = cell_type,
               output_base_path = file.path(output_base_path, "Results", "Figures", "Volcano", "nebula",
                                            paste0(output_prefix, cell_type_lower)))
  
  # GSEA Analysis
  # Load pathway files
  gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
  kegg_legacy <- prepare_gmt(gmt_files[1], unique(res_combined$Gene), savefile = FALSE)
  reactome <- prepare_gmt(gmt_files[3], unique(res_combined$Gene), savefile = FALSE)
  go <- prepare_gmt(gmt_files[4], unique(res_combined$Gene), savefile = FALSE)
  
  # Rank genes by logFC
  rankings <- res_combined$logFC_groupType_1_Diabetes
  names(rankings) <- res_combined$Gene
  rankings <- sort(rankings, decreasing = TRUE)
  
  # Run GSEA
  set.seed(1234)
  
  kegg_legacy_res <- fgsea(pathways = kegg_legacy,
                           stats = rankings,
                           scoreType = 'std', 
                           minSize = 3,
                           maxSize = 500,
                           nproc = 1)
  
  reactome_res <- fgsea(pathways = reactome,
                        stats = rankings,
                        scoreType = 'std', 
                        minSize = 3,
                        maxSize = 500,
                        nproc = 1)
  
  go_res <- fgsea(pathways = go,
                  stats = rankings,
                  scoreType = 'std', 
                  minSize = 5,
                  maxSize = 500,
                  nproc = 1)
  
  # Summary table
  fgsea_summary <- data.frame(
    "KEGG_Legacy" = c(sum(kegg_legacy_res[, padj < 0.05], na.rm = TRUE), 
                      sum(kegg_legacy_res[, pval < 0.05], na.rm = TRUE)),
    "REACTOME" = c(sum(reactome_res[, padj < 0.05], na.rm = TRUE), 
                   sum(reactome_res[, pval < 0.05], na.rm = TRUE)),
    "GO" = c(sum(go_res[, padj < 0.05], na.rm = TRUE), 
             sum(go_res[, pval < 0.05], na.rm = TRUE))
  )
  rownames(fgsea_summary) <- c("adj.pval", "p.val")
  
  # Plot pathways
  # KEGG
  plot_fgsea_transpose(kegg_legacy_res, 
                       title = paste(cell_type, "Top 30 KEGG Pathways"), 
                       xmin = 1, xmax = 3)
  
  ggsave(file.path(output_base_path, "Results", "Figures", "Pathways", "nebula",
                   paste0(output_prefix, cell_type_lower, "_res_top30_kegg_pathways.png")),
         width = 27.5, height = 14, scale = 1)
  
  # REACTOME
  plot_fgsea_transpose(reactome_res, 
                       title = paste(cell_type, "Top 30 REACTOME Pathways"), 
                       xmin = 1.45, xmax = 2.6)
  
  ggsave(file.path(output_base_path, "Results", "Figures", "Pathways", "nebula",
                   paste0(output_prefix, cell_type_lower, "_res_top30_reactome_pathways.png")),
         width = 27.5, height = 14, scale = 1)
  
  # GO
  plot_fgsea_transpose(go_res, 
                       title = paste(cell_type, "Top 30 GO Pathways"), 
                       xmin = 1.4, xmax = 5)
  
  ggsave(file.path(output_base_path, "Results", "Figures", "Pathways", "nebula",
                   paste0(output_prefix, cell_type_lower, "_res_top30_go_pathways.png")),
         width = 27.5, height = 14, scale = 1)
  
  # Return results
  return(list(
    nebula_results = res_combined,
    fgsea_summary = fgsea_summary,
    kegg_results = kegg_legacy_res,
    reactome_results = reactome_res,
    go_results = go_res,
    rankings = rankings
  ))
}

# slingshot
# ===========================================================================
# Function: slingshot_setup
# ===========================================================================

slingshot_setup <- function(seurat_obj, celltype_prefix) {
  library(SingleCellExperiment)
  library(slingshot)
  library(uwot)
  library(mclust)
  library(ggplot2)
  library(RColorBrewer)
  
  # Extract counts and convert to SCE
  counts <- GetAssayData(seurat_obj, layer = "counts")
  sce <- SingleCellExperiment(assays = List(counts = counts), colData = seurat_obj@meta.data)
  
  # Gene filtering
  gene_filter <- apply(assays(sce)$counts, 1, function(x) sum(x >= 3) >= 10)
  sce <- sce[gene_filter, ]
  
  # Normalization using quantile normalization (FQ)
  FQnorm <- function(counts){
    rk <- apply(counts, 2, rank, ties.method = 'min')
    counts.sort <- apply(counts, 2, sort)
    refdist <- apply(counts.sort, 1, median)
    norm <- apply(rk, 2, function(r) { refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }
  assays(sce)$norm <- FQnorm(assays(sce)$counts)
  
  # PCA
  pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
  
  # Variance explained
  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  
  # Elbow plot as ggplot
  elbow_df <- data.frame(PC = 1:50, VarianceExplained = var_explained[1:50])
  elbow_plot <- ggplot(elbow_df, aes(x = PC, y = VarianceExplained)) +
    geom_point(size = 2) +
    geom_line() +
    labs(title = paste("Elbow Plot:", celltype_prefix),
         x = "Principal Component",
         y = "Proportion of Variance Explained") +
    theme_bw()
  
  # Return necessary outputs
  return(list(
    sce = sce,
    pca = pca,
    var_explained = var_explained,
    elbow_plot = elbow_plot
  ))
}

# ===========================================================================
# Function: run_slingshot
# ===========================================================================

run_slingshot <- function(sce, pca_obj, n_pcs = 6, start_cluster = NULL, end_cluster = NULL, cluster_label = "celltype") {
  library(slingshot)
  library(uwot)
  
  # Get top PCs
  rd1 <- pca_obj$x[, 1:n_pcs]
  
  # Visualize PCA (optional)
  plot(rd1, col = rgb(0, 0, 0, 0.5), pch = 16, asp = 1, main = "PCA")
  
  # Run UMAP
  umap_mat <- uwot::umap(t(log1p(assays(sce)$norm)))
  colnames(umap_mat) <- c("UMAP1", "UMAP2")
  
  # Visualize UMAP (optional)
  plot(umap_mat, col = rgb(0, 0, 0, 0.5), pch = 16, asp = 1, main = "UMAP")
  
  # Add to reducedDims
  reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = umap_mat)
  
  # Run Slingshot
  sce_sl <- slingshot(sce, clusterLabels = cluster_label, reducedDim = "PCA", 
                      start.clus = start_cluster,
                      end.clus = end_cluster)
  
  return(sce_sl)
}

# ===========================================================================
# Function: plot_slingshot_trajectory
# ===========================================================================

plot_slingshot_trajectory <- function(sce_sl, 
                                      celltype_levels, 
                                      custom_colors = color_5, 
                                      cluster_label = "celltype",
                                      bucket = "attempt",
                                      celltype_suffix = NULL,
                                      title = "Slingshot trajectory",
                                      lineage = 1) {
  library(ggplot2)
  library(RColorBrewer)
  library(slingshot)
  
  # Extract PCA and pseudotime
  pca_df <- as.data.frame(reducedDims(sce_sl)$PCA)
  
  # Dynamically extract cluster labels
  pca_df$cluster_label <- factor(colData(sce_sl)[[cluster_label]], levels = celltype_levels)
  pca_df$pseudotime <- slingPseudotime(sce_sl)[, lineage]
  
  # Get curve
  curves <- slingCurves(sce_sl)
  curve_df <- as.data.frame(curves[[lineage]]$s)
  lineage_obj <- slingshot::slingLineages(sce_sl)
  
  # Generate ggplot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster_label)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_path(data = curve_df, aes(x = PC1, y = PC2), 
              color = "black", size = 1.2, inherit.aes = FALSE) +
    # theme_bw() + 
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank(),
          plot.background = element_blank()) + 
    labs(color = NULL, title = title) +
    scale_color_manual(values = custom_colors)
  
  # Save if bucket is provided
  if (!is.null(bucket)) {
    temp_file <- tempfile(fileext = ".png")
    ggsave(filename = temp_file, plot = p, width = 7, height = 5, bg = "transparent")
    s3$upload_file(temp_file, bucket, paste0("slingshot/attempt_pca_", tolower(celltype_suffix), "_slingshot.png"))
  }
  
  return(list(
    pca_plot = p,
    curves = curves,
    pca_df = pca_df,
    curve_df = curve_df,
    lineage = lineage_obj
  ))
}
# ===========================================================================
# Function: create_pseudotime_df
# ===========================================================================

create_pseudotime_df <- function(sce_obj, lineage = 1) {
  # Check if slingPseudotime exists
  if (!"slingPseudotime_1" %in% names(colData(sce_obj))) {
    stop("The input object does not contain slingPseudotime. Run Slingshot first.")
  }
  
  # Extract relevant metadata and pseudotime
  pseudotime_df <- data.frame(
    pseudotime = slingPseudotime(sce_obj)[, lineage],
    treatment = sce_obj$treatment,
    visit = sce_obj$visit,
    celltype = sce_obj$celltype,
    visit_treatment = paste(sce_obj$visit, sce_obj$treatment)
  )
  
  # Reorder factor levels
  pseudotime_df$visit_treatment <- factor(pseudotime_df$visit_treatment, 
                                          levels = c("PRE Placebo", "POST Placebo", 
                                                     "PRE Dapagliflozin", "POST Dapagliflozin"))
  return(pseudotime_df)
}


# ===========================================================================
# Function: plot_pseudotime_violin
# ===========================================================================

plot_pseudotime_violin <- function(df, s3_folder = "slingshot", 
                                   celltype_suffix = "celltype") {
  pseudotime_df <- df
  # Calculate medians
  median_df <- pseudotime_df %>%
    group_by(visit_treatment) %>%
    dplyr::summarise(median_pt = median(pseudotime, na.rm = TRUE), .groups = "drop")
  
  # Define colors
  fill_colors <- c("PRE Placebo" = "#fbc4ab", "POST Placebo" = "#f4978e",
                   "PRE Dapagliflozin" = "#ccd5ae", "POST Dapagliflozin" = "#828e82")
  
  # Create plot
  p <- ggplot(pseudotime_df, aes(x = visit_treatment, y = pseudotime)) +
    geom_violin(aes(fill = visit_treatment, color = visit_treatment), trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, aes(color = visit_treatment)) +
    geom_line(data = subset(median_df, grepl("Placebo", visit_treatment)), 
              aes(x = visit_treatment, y = median_pt, group = 1), 
              color = "#343a40", linewidth = 0.5, linetype = "dashed") +
    geom_line(data = subset(median_df, grepl("Dapagliflozin", visit_treatment)), 
              aes(x = visit_treatment, y = median_pt, group = 1), 
              color = "#343a40", linewidth = 0.5, linetype = "dashed") +
    theme_classic() +
    labs(x = NULL, y = "Pseudotime", fill = NULL, color = NULL) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors) +
    theme(legend.position = "none",
          text = element_text(size = 15),
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  # Save and upload to S3
  temp_file <- tempfile(fileext = ".png")
  ggsave(filename = temp_file, plot = p, width = 7, height = 5)
  s3$upload_file(temp_file, "attempt", file.path(s3_folder, paste0(tolower(celltype_suffix), "_attempt_slingshot_violin.png")))
  return(p)
}

# ===========================================================================
# Function: plot_slingshot_3d
# ===========================================================================

library(plotly)
library(htmlwidgets)

plot_slingshot_3d <- function(pca_df, curve_df, custom_colors = color_5,
                              s3_folder = "slingshot", 
                              cluster_label = "celltype",
                              celltype_suffix = "celltype") {
  
  # Create plot
  pt_plotly <- plot_ly() %>%
    add_trace(
      data = pca_df,
      x = ~PC1, y = ~PC2, z = ~PC3,
      type = "scatter3d",
      mode = "markers",
      color = ~cluster_label,
      colors = custom_colors,
      marker = list(size = 2),
      name = "Cells") %>%
    add_trace(
      data = curve_df,
      x = ~PC1, y = ~PC2, z = ~PC3,
      type = "scatter3d",
      mode = "lines",
      line = list(color = "black", width = 4),
      name = "Trajectory")
  
  # Save and upload to S3
  temp_file <- tempfile(fileext = ".html")
  saveWidget(pt_plotly, file = temp_file, selfcontained = TRUE)
  s3$upload_file(temp_file, "attempt", file.path(s3_folder, paste0("attempt_", tolower(celltype_suffix), "_slingshot_plotly.html")))
  return(pt_plotly)
}

# ===========================================================================
# Function: plot_and_test_pseudotime_distribution
# ===========================================================================

plot_and_test_pseudotime_distribution <- function(df,
                                                  sce_object,
                                                  pseudotime_var = "pseudotime",
                                                  visit_treatment_var = "visit_treatment",
                                                  s3_folder = "slingshot",
                                                  filename_suffix = "celltype") {
  library(condiments)
  BiocParallel::register(BiocParallel::SerialParam())
  # Set colors
  group_colors <- c("PRE Placebo" = "#fbc4ab",
                    "POST Placebo" = "#f4978e",
                    "PRE Dapagliflozin" = "#ccd5ae",
                    "POST Dapagliflozin" = "#828e82")
  
  # Create plot
  p <- ggplot(df, aes(x = .data[[pseudotime_var]],
                      fill = .data[[visit_treatment_var]],
                      color = .data[[visit_treatment_var]])) +
    geom_density(aes(weight = weight_per_cell), alpha = 0.3) + 
    theme_minimal() +
    labs(x = "Pseudotime", y = "Density", fill = NULL, color = NULL) +
    theme(
      legend.position = c(0.5, 0.85),
      legend.direction = "horizontal",
      legend.text = element_text(size = 10),
      panel.grid = element_blank(),
      text = element_text(size = 15)
    ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    facet_wrap(~treatment, nrow = 2, ncol = 1)
  
  # Save and upload
  temp_file <- tempfile(fileext = ".png")
  ggsave(filename = temp_file, plot = p, width = 7, height = 5)
  s3$upload_file(temp_file, "attempt", file.path(s3_folder, paste0("attempt_", filename_suffix, "_slingshot_density_trtvisit.png")))
  
  # Run progressionTest
  registerDoSEQ()
  test_result <- progressionTest(sce_object, conditions = df[[visit_treatment_var]])
  print(test_result)
  return(list(plot = p, test_result = test_result))
}

# ===========================================================================
# Function: plot_pseudotime_density_faceted_by_treatment
# ===========================================================================

plot_pseudotime_density_faceted_by_treatment <- function(df,
                                                         pseudotime_var = "slingPseudotime_1",
                                                         visit_col = "visit",
                                                         treatment_col = "treatment",
                                                         visit_treatment_col = "visit_treatment",
                                                         s3_folder = "slingshot",
                                                         filename_suffix = "celltype") {
  # Quantile summary
  raw_summary <- df %>%
    as.data.frame() %>%
    group_by(.data[[treatment_col]], .data[[visit_col]]) %>%
    dplyr::summarise(
      p25 = quantile(.data[[pseudotime_var]], probs = 0.25, na.rm = TRUE),
      p50 = quantile(.data[[pseudotime_var]], probs = 0.5, na.rm = TRUE),
      p75 = quantile(.data[[pseudotime_var]], probs = 0.75, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(visit_treatment = paste(.data[[visit_col]], .data[[treatment_col]]))
  
  # Color palette
  fill_colors <- c("PRE Placebo" = "#fbc4ab", 
                   "POST Placebo" = "#f4978e",
                   "PRE Dapagliflozin" = "#ccd5ae", 
                   "POST Dapagliflozin" = "#828e82")
  
  # Plot
  p <- ggplot(df, aes(x = .data[[pseudotime_var]], fill = .data[[visit_treatment_col]])) +
    geom_density(alpha = 0.5, aes(color = .data[[visit_treatment_col]], weight = weight_per_cell)) +
    facet_wrap(vars(.data[[treatment_col]]), strip.position = "bottom") +
    geom_vline(data = raw_summary, aes(xintercept = p25, color = visit_treatment), linetype = "dashed") +
    geom_vline(data = raw_summary, aes(xintercept = p50, color = visit_treatment), linetype = "dashed") +
    geom_vline(data = raw_summary, aes(xintercept = p75, color = visit_treatment), linetype = "dashed") +
    theme_minimal() +
    labs(x = "Pseudotime", y = "Density", color = NULL, fill = NULL) +
    theme(panel.grid = element_blank(),
          text = element_text(size = 15),
          legend.position = c(0.45, 0.95),
          legend.direction = "horizontal",
          legend.text = element_text(size = 12)) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors)
  
  # Save and upload
  temp_file <- tempfile(fileext = ".png")
  ggsave(filename = temp_file, plot = p, width = 7, height = 5)
  s3$upload_file(temp_file, "attempt", file.path(s3_folder, paste0("attempt_", filename_suffix, "_slingshot_density_trt.png")))
  
  return(p)
}

# ===========================================================================
# Function: plot_delta_percentile_heatmap
# ===========================================================================

plot_delta_percentile_heatmap <- function(df,
                                          percentile_prefix = "rq_tau_",
                                          visit_col = "visit",
                                          treatment_col = "treatment",
                                          s3_folder = "slingshot",
                                          filename_suffix = "celltype",
                                          width = 5,
                                          height = 5) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Prepare data for heatmap
  heatmap_df <- df %>%
    pivot_longer(cols = starts_with(percentile_prefix),
                 names_to = "q",
                 names_prefix = percentile_prefix,
                 values_to = "value") %>%
    pivot_wider(names_from = .data[[visit_col]], values_from = value) %>%
    dplyr::mutate(delta = POST - PRE)
  
  # Plot
  p <- ggplot(heatmap_df, aes(x = .data[[treatment_col]], y = q, fill = delta)) +
    geom_tile() +
    scale_fill_gradient2(low = "#89c2d9", mid = "white", high = "#ee7674", midpoint = 0) +
    labs(x = NULL, y = "Percentile", fill = "\u0394\n(POST-PRE)") +
    theme_minimal() +
    theme(
      legend.title.position = "top",
      legend.title = element_text(hjust = 0.5, size = 10),
      legend.position = c(1.08, 0.5),
      panel.grid = element_blank(),
      text = element_text(size = 15),
      legend.text = element_text(size = 10),
      plot.margin = unit(c(0, 1.8, 0, 0), "cm")
    )
  
  # Save and upload
  temp_file <- tempfile(fileext = ".png")
  ggsave(filename = temp_file, plot = p, width = width, height = height)
  s3$upload_file(temp_file, "attempt", file.path(s3_folder, paste0("attempt_", filename_suffix, "_slingshot_percentile_heatmap.png")))
  
  return(p)
}

# ===========================================================================
# Function: analyze_pseudotime_by_clinvar
# ===========================================================================

library(dplyr)
library(ggplot2)

analyze_pseudotime_by_clinvar <- function(df,
                                          clinical_var,          # unquoted column name of clinical variable
                                          pseudotime_var,        # unquoted column name of pseudotime
                                          subject_id = "subject_id",
                                          visit_col = "visit",
                                          pre_label = "PRE",
                                          post_label = "POST",
                                          bin_probs = 2,
                                          bin_levels_to_compare = c(1, 2),
                                          caption_clinical_var = "",
                                          bucket = "attempt",
                                          celltype_suffix = "",
                                          filesuffix = "",
                                          plot_by_visit_direction = TRUE) {  # New parameter
  theme_transparent <- theme(
    plot.background   = element_rect(fill = "transparent", color = NA),
    panel.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
  clinical_var <- ensym(clinical_var)
  pseudotime_var <- ensym(pseudotime_var)
  clinical_var_chr <- rlang::as_string(clinical_var)
  
  # Step 1: Calculate subject-level delta of the clinical variable
  delta_df <- df %>%
    group_by(.data[[subject_id]]) %>%
    dplyr::summarise(
      delta_value = mean(.data[[clinical_var_chr]][.data[[visit_col]] == post_label], na.rm = TRUE) -
        mean(.data[[clinical_var_chr]][.data[[visit_col]] == pre_label], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Step 2: Join delta back to the original dataframe
  df <- df %>%
    left_join(delta_df, by = subject_id)
  
  # Step 3: Add direction column based on delta
  df <- df %>%
    dplyr::mutate(direction = case_when(
      .data[[clinical_var_chr]] > 0 ~ "+",
      .data[[clinical_var_chr]] < 0 ~ "-",
      .data[[clinical_var_chr]] == 0 ~ "No Change",
      TRUE ~ NA_character_
    ))
  
  # Step 4: Create visit_direction combination
  df <- df %>%
    dplyr::mutate(visit_direction = paste0(.data[[visit_col]], "/", direction))
  
  # Step 5: Bin the clinical variable by quantiles (original binning)
  df_binned <- if (identical(bin_probs, "direction")) {
    df %>%
      dplyr::mutate(clinical_bin = case_when(
        .data[[clinical_var_chr]] > 0 ~ "+",
        .data[[clinical_var_chr]] < 0 ~ "-",
        .data[[clinical_var_chr]] == 0 ~ "No Change"
      )) %>%
      dplyr::mutate(clinical_bin = factor(clinical_bin, levels = c("-", "No Change", "+")))
  } else {
    df %>%
      dplyr::mutate(clinical_bin = cut(.data[[clinical_var_chr]],
                                       breaks = quantile(.data[[clinical_var_chr]], 
                                                         probs = seq(0, 1, 1/bin_probs), na.rm = TRUE),
                                       include.lowest = TRUE)) %>%
      filter(!is.na(clinical_bin))
  }
  
  # Step 6: Plot density of pseudotime by clinical bins (original plot)
  n_bins <- nlevels(droplevels(df_binned$clinical_bin))
  custom_colors <- NULL
  if (n_bins == 2) {
    custom_colors <- c("#457b9d", "#f28482")
    names(custom_colors) <- levels(droplevels(df_binned$clinical_bin))
  }
  
  p_original <- df_binned %>%
    filter(!is.na(clinical_bin)) %>%
    ggplot(aes(x = !!pseudotime_var, fill = clinical_bin)) +
    geom_density(aes(weight = weight_per_cell), alpha = 0.5) +
    theme_minimal() +
    labs(
      y = "Density", x = "Pseudotime", fill = NULL, color = NULL,
      caption = paste0(caption_clinical_var, " across pseudotime")) +
    theme(panel.grid = element_blank(), 
          text = element_text(size = 15),
          legend.position = c(0.5, 0.95),
          legend.direction = "horizontal") +
    theme_transparent + 
    {
      if (!is.null(custom_colors)) {
        list(scale_fill_manual(values = custom_colors))
      } else {
        NULL
      }
    }
  
  # Step 7: New plot - density by visit and direction
  p_visit_direction <- NULL
  if (plot_by_visit_direction) {
    # Define colors for visit/direction combinations
    visit_dir_colors <- c(
      "PRE/+" = "#d6604d",    # Dark red
      "PRE/-" = "#2166ac",    # Dark blue
      "POST/+" = "#f4a582",   # Light red
      "POST/-" = "#92c5de",   # Light blue
      "PRE/No Change" = "#888888",   # Gray
      "POST/No Change" = "#bbbbbb"   # Light gray
    )
    
    # Filter to valid combinations
    df_visit_dir <- df_binned %>%
      filter(!is.na(direction), !is.na(.data[[visit_col]]))
    
    # Get actual levels present in data
    actual_levels <- unique(df_visit_dir$visit_direction)
    colors_to_use <- visit_dir_colors[names(visit_dir_colors) %in% actual_levels]
    
    p_visit_direction <- df_visit_dir %>%
      ggplot(aes(x = !!pseudotime_var, fill = visit_direction)) +
      geom_density(aes(weight = weight_per_cell), alpha = 0.5) +
      theme_minimal() +
      labs(
        y = "Density", 
        x = "Pseudotime", 
        fill = "Visit/Direction", 
        color = NULL,
        caption = paste0("Pseudotime density by visit and ", caption_clinical_var, " direction")
      ) +
      scale_fill_manual(values = colors_to_use) +
      theme(
        panel.grid = element_blank(), 
        text = element_text(size = 14),
        legend.position = "right",
      ) +
      theme_transparent
    
    # Optionally add faceting for clearer visualization
    p_visit_direction_faceted <- df_visit_dir %>%
      ggplot(aes(x = !!pseudotime_var, fill = visit)) +
      geom_density(aes(weight = weight_per_cell), alpha = 0.6) + 
      facet_wrap(~ direction, ncol = 1) +
      theme_minimal() +
      labs(
        y = "Density", 
        x = "Pseudotime", 
        fill = "Direction", 
        color = NULL,
        caption = paste0("Pseudotime density by visit and ", caption_clinical_var, " direction")
      ) +
      scale_fill_manual(values = c("POST" = "#f28482", "PRE" = "#457b9d", "No Change" = "#888888")) +
      theme(
        panel.grid = element_blank(), 
        text = element_text(size = 14),
        legend.position = "top",
        strip.text = element_text(size = 12, face = "bold")
      ) +
      theme_transparent
    
    # Save the visit/direction plot if bucket is specified
    if (!is.null(bucket)) {
      temp_file_vd <- tempfile(fileext = ".png")
      ggsave(p_visit_direction, filename = temp_file_vd, width = 9, height = 6)
      s3$upload_file(temp_file_vd, bucket, 
                     paste0("slingshot/attempt_density_visitdir_", 
                            tolower(celltype_suffix), "_", clinical_var_chr, 
                            filesuffix, "_slingshot.png"))
      
      # Save faceted version too
      temp_file_vd_facet <- tempfile(fileext = ".png")
      ggsave(p_visit_direction_faceted, filename = temp_file_vd_facet, width = 8, height = 8)
      s3$upload_file(temp_file_vd_facet, bucket, 
                     paste0("slingshot/attempt_density_visitdir_faceted_", 
                            tolower(celltype_suffix), "_", clinical_var_chr, 
                            filesuffix, "_slingshot.png"))
    }
    
    # Print both new plots
    print(p_visit_direction)
    print(p_visit_direction_faceted)
  }
  
  # Save original plot if bucket is specified
  if (!is.null(bucket)) {
    temp_file <- tempfile(fileext = ".png")
    ggsave(p_original, filename = temp_file, width = 7, height = 5)
    s3$upload_file(temp_file, bucket, paste0("slingshot/attempt_density_", 
                                             tolower(celltype_suffix), "_", clinical_var_chr, 
                                             filesuffix, "_slingshot.png"))
  }
  
  print(p_original)
  
  # Step 8: KS test between selected bin levels (original functionality)
  bin_levels <- levels(droplevels(df_binned$clinical_bin))
  if (length(bin_levels) >= max(bin_levels_to_compare)) {
    pt1 <- df_binned %>%
      filter(clinical_bin == bin_levels[bin_levels_to_compare[1]]) %>%
      pull(!!pseudotime_var)
    
    pt2 <- df_binned %>%
      filter(clinical_bin == bin_levels[bin_levels_to_compare[2]]) %>%
      pull(!!pseudotime_var)
    
    ks <- ks.test(pt1, pt2)
    print(ks)
  } else {
    warning("Not enough bins to compare the requested levels.")
    ks <- NULL
  }
  
  # Step 9: Additional KS tests for visit/direction comparisons
  ks_visit_direction <- NULL
  if (plot_by_visit_direction) {
    # Compare PRE/+ vs PRE/-
    df_pre_pos <- df_binned %>% 
      filter(.data[[visit_col]] == pre_label, direction == "+")
    df_pre_neg <- df_binned %>% 
      filter(.data[[visit_col]] == pre_label, direction == "-")
    
    if (nrow(df_pre_pos) > 0 && nrow(df_pre_neg) > 0) {
      ks_pre <- ks.test(
        df_pre_pos %>% pull(!!pseudotime_var),
        df_pre_neg %>% pull(!!pseudotime_var)
      )
      cat("\nKS test PRE/+ vs PRE/-:\n")
      print(ks_pre)
    }
    
    # Compare POST/+ vs POST/-
    df_post_pos <- df_binned %>% 
      filter(.data[[visit_col]] == post_label, direction == "+")
    df_post_neg <- df_binned %>% 
      filter(.data[[visit_col]] == post_label, direction == "-")
    
    if (nrow(df_post_pos) > 0 && nrow(df_post_neg) > 0) {
      ks_post <- ks.test(
        df_post_pos %>% pull(!!pseudotime_var),
        df_post_neg %>% pull(!!pseudotime_var)
      )
      cat("\nKS test POST/+ vs POST/-:\n")
      print(ks_post)
    }
    
    # Compare across visits for same direction
    if (nrow(df_pre_pos) > 0 && nrow(df_post_pos) > 0) {
      ks_pos_visits <- ks.test(
        df_pre_pos %>% pull(!!pseudotime_var),
        df_post_pos %>% pull(!!pseudotime_var)
      )
      cat("\nKS test PRE/+ vs POST/+:\n")
      print(ks_pos_visits)
    }
    
    ks_visit_direction <- list(
      pre_comparison = if(exists("ks_pre")) ks_pre else NULL,
      post_comparison = if(exists("ks_post")) ks_post else NULL,
      positive_visits = if(exists("ks_pos_visits")) ks_pos_visits else NULL
    )
  }
  
  return(list(
    plot_original = p_original, 
    plot_visit_direction = p_visit_direction,
    plot_visit_direction_faceted = if(exists("p_visit_direction_faceted")) p_visit_direction_faceted else NULL,
    ks_test = ks,
    ks_visit_direction = ks_visit_direction,
    delta_df = delta_df, 
    df_binned = df_binned
  ))
}

# ===========================================================================
# Function: analyze_pseudotime_by_celltype
# ===========================================================================

library(quantreg)
library(slingshot)
analyze_pseudotime_by_celltype <- function(so, 
                                           celltype_col,    
                                           celltype_groups, 
                                           start_cluster,   
                                           celltype_levels = celltype_groups,
                                           custom_colors,   
                                           suffix, 
                                           n_pcs = 10, 
                                           tau = c(0.25, 0.5, 0.75),
                                           aws_s3,
                                           s3_key) {
  
  start_time <- Sys.time()
  step_time <- start_time
  
  log_step <- function(msg) {
    now <- Sys.time()
    message(sprintf("[%s] %s (elapsed: %.2f sec)", 
                    format(now, "%Y-%m-%d %H:%M:%S"), 
                    msg, as.numeric(difftime(now, step_time, units = "secs"))))
    step_time <<- now
  }
  
  log_step("Step 1: Subsetting Seurat object...")
  cells_to_keep <- rownames(so@meta.data)[so@meta.data[[celltype_col]] %in% celltype_groups]
  so_sub <- subset(so, cells = cells_to_keep)
  
  log_step("Step 2: Running slingshot setup...")
  sling_obj <- slingshot_setup(so_sub, celltype_prefix = toupper(suffix))
  sce <- sling_obj$sce
  
  log_step("Step 3: Inspecting elbow plot...")
  print(sling_obj$elbow_plot)
  
  log_step("Step 4: Running Slingshot...")
  sce_sl <- run_slingshot(sce, sling_obj$pca, 
                          cluster_label  = celltype_col,
                          n_pcs = n_pcs,
                          start_cluster = start_cluster)
  
  log_step("Step 5: Plotting trajectory...")
  res <- plot_slingshot_trajectory(sce_sl = sce_sl,
                                   cluster_label = celltype_col,
                                   celltype_levels = celltype_levels,
                                   custom_colors = custom_colors,
                                   celltype_suffix = suffix)
  print(res$pca_plot)
  print(res$lineage)
  
  log_step("Step 6: Creating pseudotime dataframe...")
  pseudotime_df <- create_pseudotime_df(sce_sl)
  
  log_step("Step 7: Plotting pseudotime violin...")
  plot_pseudotime_violin(df = pseudotime_df, celltype_suffix = suffix)
  
  log_step("Step 8: Plotting 3D trajectory...")
  plot_slingshot_3d(pca_df = res$pca_df,
                    curve_df = res$curve_df,
                    celltype_suffix = suffix)
  
  log_step("Step 9: Preparing data for quantile regression...")
  no_s4 <- setdiff(names(colData(sce_sl)), "slingshot")
  sce_df <- as.data.frame(colData(sce_sl)[, no_s4]) %>%
    dplyr::mutate(visit_treatment = factor(paste(visit, treatment), 
                                           levels = c("PRE Placebo", "POST Placebo",
                                                      "PRE Dapagliflozin", "POST Dapagliflozin"))) %>%
    group_by(subject_id, visit) %>%
    mutate(weight_per_cell = 1/n()) %>%
    ungroup()
  
  sce_sl$weight_per_cell <- sce_df$weight_per_cell
  pseudotime_df$weight_per_cell <- sce_df$weight_per_cell
  
  log_step("Step 10: Running quantile regression...")
  
  rq_fit <- rq(slingPseudotime_1 ~ treatment * visit, 
               tau = tau, 
               data = sce_df)
  
  attr(rq_fit$terms, ".Environment") <- NULL
  attr(rq_fit$formula, ".Environment") <- NULL
  environment(attr(rq_fit$model, "terms")) <- baseenv()
  
  log_step("Step 10: Saving rq_fit...")
  
  temp_file <- tempfile(fileext = ".rds")
  tryCatch({
    R.utils::withTimeout({
      saveRDS(rq_fit, temp_file)
      message(sprintf("Actual file size: %.2f MB", file.size(temp_file) / 1024^2))
    }, timeout = 30, onTimeout = "error")
  }, error = function(e) {
    message("Save timed out - confirms large serialization issue")
  })
  aws_s3$upload_file(temp_file, Bucket = "attempt", 
                     Key = paste0("Results/", suffix, "_rq_fit.rds"))
  if (file.exists(temp_file)) unlink(temp_file)
  
  log_step("Step 10b: Running quantile regression for clinical variables...")
  
  clin_vars <- c(
    "hba1c_delta" = "\u0394HbA1c",
    "weight_delta" = "\u0394Weight", 
    "mgfr_jodal_delta" = "\u0394mGFR",
    "mgfr_jodal_bsa_delta" = "\u0394mGFR (BSA)",
    "tir_delta" = "\u0394TIR",
    "avg_ketones_delta" =  "\u0394Ketones"
  )
  
  for (clin_var in names(clin_vars)) {
    formula_clin <- reformulate(clin_var, response = "slingPseudotime_1")
    rq_fit_clin <- rq(formula_clin, tau = tau, data = sce_df)
    
    # Clean environments for saving
    attr(rq_fit_clin$terms, ".Environment") <- NULL
    attr(rq_fit_clin$formula, ".Environment") <- NULL
    environment(attr(rq_fit_clin$model, "terms")) <- baseenv()
    
    # Save to temp file and upload
    temp_file_clin <- tempfile(fileext = ".rds")
    tryCatch({
      R.utils::withTimeout({
        saveRDS(rq_fit_clin, temp_file_clin)
        message(sprintf("   ✅ Saved rq for %s (%.2f MB)", clin_var, file.size(temp_file_clin) / 1024^2))
      }, timeout = 30, onTimeout = "error")
    }, error = function(e) {
      message(sprintf(" ❌ Save timed out for %s", clin_var))
    })
    
    s3_key_rq <- paste0("Results/", suffix, "_rq_fit_", clin_var, ".rds")
    aws_s3$upload_file(temp_file_clin, Bucket = "attempt", Key = s3_key_rq)
    if (file.exists(temp_file_clin)) unlink(temp_file_clin)
    
    message(sprintf("⬆️  Clinical rq uploaded to s3://%s/%s", "attempt", s3_key_rq))
  }
  log_step("Step 11: Predicting values from quantile regression...")
  rq_grid <- expand.grid(treatment = c("Placebo", "Dapagliflozin"),
                         visit = c("PRE", "POST"))
  rq_pred <- predict(rq_fit, rq_grid)
  tau_names <- paste0("rq_tau_", gsub("\\.", "", as.character(tau * 100)))
  rq_grid <- cbind(rq_grid, setNames(as.data.frame(rq_pred), tau_names))
  
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(rq_grid, temp_file)
  aws_s3$upload_file(temp_file, Bucket = "attempt", 
                     Key = paste0("Results/", suffix, "_rq_predict.rds"))
  unlink(temp_file)
  
  message(sprintf("⬆️  Results uploaded to s3://%s/%s", "attempt", paste0("Results/", suffix, "_rq_predict.rds")))
  
  log_step("Step 12: Plotting density plots...")
  p_pseudotime <- plot_and_test_pseudotime_distribution(pseudotime_df, sce_sl, filename_suffix = suffix)
  p_pseudotime_trt <- plot_pseudotime_density_faceted_by_treatment(sce_df, filename_suffix = suffix)
  p_pseudotime_heatmap <- plot_delta_percentile_heatmap(rq_grid, filename_suffix = suffix)
  
  log_step("Step 13: Running delta clinical variable analyses...")
  
  clin_var_res <- list()
  for (clin_var in names(clin_vars)) {
    message(sprintf("   [%s]   Analyzing %s...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), clin_var))
    res <- analyze_pseudotime_by_clinvar(sce_df, 
                                         !!sym(clin_var), 
                                         slingPseudotime_1,
                                         bin_probs = "direction",
                                         caption_clinical_var = clin_vars[clin_var],
                                         celltype_suffix = suffix,
                                         filesuffix = "delta",
                                         bucket = "attempt")
    clin_var_res[[clin_var]] <- res
  }
  
  log_step("All steps completed successfully!")
  message(sprintf("Total runtime: %.2f sec", as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
  
  return(list(sce = sce_sl, sce_df = sce_df, 
              pseudotime_df = pseudotime_df, 
              rq_fit = rq_fit, rq_predict = rq_grid,
              p_pseudotime = p_pseudotime,
              p_pseudotime_trt = p_pseudotime_trt, 
              p_pseudotime_heatmap = p_pseudotime_heatmap,
              clin_var_res = clin_var_res))
}

# ===========================================================================
# Function: plot_treatment_heatmap
# ===========================================================================
# 

plot_treatment_heatmap <- function(data, 
                                   top_n = 20,
                                   p_col = "p_treatmentDapagliflozin:visitPOST",
                                   logfc_col = "logFC_treatmentDapagliflozin:visitPOST",
                                   fdr_col = "fdr",
                                   heatmap_caption = T,
                                   gene_col = "Gene",
                                   logfc_visit = "logFC_visitPOST",
                                   genes_to_show = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtext)  # Add this for element_markdown
  
  if (!is.null(genes_to_show)) {
    data <- data %>% filter(!!sym(gene_col) %in% genes_to_show)  # Fixed: use gene_col
  }
  
  data <- data %>%
    filter(abs(!!sym(logfc_col)) < 10)
  
  # Identify top positive and negative genes
  top_pos <- data %>%
    dplyr::arrange(!!sym(fdr_col)) %>%
    filter(!!sym(logfc_col) > 0) %>%
    head(top_n) %>%
    pull(!!sym(gene_col))
  
  top_neg <- data %>%
    dplyr::arrange(!!sym(fdr_col)) %>%
    filter(!!sym(logfc_col) < 0) %>%
    head(top_n) %>%
    pull(!!sym(gene_col))
  
  # Prepare long-form data
  heatmap_df <- data %>%
    dplyr::mutate(
      Placebo = !!sym(logfc_visit),
      Dapagliflozin = !!sym(logfc_visit) + !!sym(logfc_col),
      DiD = !!sym(logfc_col)
    ) %>%
    dplyr::select(!!sym(gene_col), Placebo, Dapagliflozin, DiD) %>%
    pivot_longer(cols = c(Placebo, Dapagliflozin, DiD),
                 names_to = "group",
                 values_to = "Slope") %>%
    dplyr::mutate(group = factor(group, levels = c("DiD", "Placebo", "Dapagliflozin")))
  
  # Order genes by DiD slope
  plot_df <- heatmap_df %>%
    filter(!!sym(gene_col) %in% c(top_pos, top_neg)) %>%
    group_by(!!sym(gene_col)) %>%
    dplyr::mutate(DiD_value = Slope[group == "DiD"]) %>%
    ungroup() %>%
    dplyr::arrange(desc(DiD_value)) %>%
    dplyr::mutate(!!sym(gene_col) := reorder(!!sym(gene_col), -DiD_value)) %>%
    dplyr::select(-DiD_value)
  
  abs_max <- max(abs(plot_df$Slope)) * 1.1
  
  gene_levels <- c(top_neg, top_pos)
  
  y_colors <- setNames(
    c(rep("#457b9d", top_n), rep("#f28482", top_n)),
    gene_levels
  )
  
  # Fixed: use gene_col instead of hardcoded "Gene"
  plot_df[[gene_col]] <- factor(plot_df[[gene_col]], levels = rev(gene_levels))
  
  # Heatmap - Fixed to use gene_col parameter
  ggplot(plot_df, aes(x = group, y = !!sym(gene_col), fill = Slope)) +
    geom_tile() +
    scale_fill_gradient2(low = "#89c2d9", mid = "white", high = "#ee7674",
                         midpoint = 0, limits = c(-abs_max, abs_max)) +
    scale_y_discrete(
      position = "right", 
      labels = function(x) {
        mapply(function(label, col) {
          glue::glue("<span style='color:{col}'>{label}</span>")
        }, x, y_colors[x], SIMPLIFY = TRUE)
      }
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = 15),
      legend.title.position = "top",
      legend.title = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(angle = 60, hjust = 1), 
      axis.text.y = element_markdown(face = "bold"),
      plot.caption = element_text(size = 18, hjust = 0.5)
    ) +
    labs(
      x = NULL, y = NULL,
      caption = if(heatmap_caption) {
        expression(
          atop(
            atop("", "Difference-in-Differences: " * beta[3]),
            atop(
              "Placebo slope: (" * beta[0] + beta[1] * ") - " * beta[0] * " = " * beta[1],
              "Dapagliflozin slope: (" * beta[0] + beta[1] + beta[2] + beta[3] *
                ") - (" * beta[0] + beta[2] * ") = " * beta[1] + beta[3]
            )
          )
        )
      }
      else NULL
    )
}

# ===========================================================================
# Function: plot_clinvar_pseudotime_arrows
# ===========================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(rlang)
library(grid)

plot_clinvar_pseudotime_arrows <- function(df,
                                           pseudotime_var = "slingPseudotime_1",
                                           subject_col = "subject_id",
                                           visit_col = "visit",
                                           treatment_col = "treatment",
                                           visit_treatment_col = "visit_treatment",
                                           clinical_var,                # unquoted, e.g. hba1c
                                           clinical_var_label = NULL,   # for axis label
                                           percentile_filter = 50,
                                           percentiles = c(25, 50, 65, 85),
                                           shape_pre = 16,
                                           shape_post = 17,
                                           celltype_suffix = "",
                                           bucket = "attempt",
                                           color_palette = c("PRE Placebo" = "#fbc4ab",
                                                             "POST Placebo" = "#f4978e",
                                                             "PRE Dapagliflozin" = "#ccd5ae",
                                                             "POST Dapagliflozin" = "#828e82",
                                                             "Consistent" = "#160f29",
                                                             "Inconsistent" = "#e5e5e5")) {
  
  clinical_var <- rlang::ensym(clinical_var)
  pseudotime_var_chr <- rlang::as_string(rlang::ensym(pseudotime_var))
  clinical_var_chr <- rlang::as_string(clinical_var)
  
  # Step 1: Pseudotime percentiles
  med <- df %>%
    group_by(.data[[subject_col]], .data[[visit_col]]) %>%
    dplyr::summarise(across(
      .data[[pseudotime_var_chr]],
      list("25" = ~quantile(., 0.25, na.rm = TRUE),
           "50" = ~quantile(., 0.50, na.rm = TRUE),
           "65" = ~quantile(., 0.65, na.rm = TRUE),
           "85" = ~quantile(., 0.85, na.rm = TRUE)),
      .names = "pseudotime_{fn}"
    ), .groups = "drop")
  
  # Step 2: Construct plot_df
  plot_df <- df %>%
    dplyr::select(all_of(c(subject_col, visit_col, treatment_col, clinical_var_chr, visit_treatment_col))) %>%
    distinct(.data[[subject_col]], .data[[visit_col]], .keep_all = TRUE) %>%
    left_join(med, by = c(subject_col, visit_col)) %>%
    pivot_longer(
      cols = starts_with("pseudotime_"),
      names_to = "percentile",
      names_prefix = "pseudotime_",
      values_to = "pseudotime"
    ) %>%
    dplyr::mutate(percentile = as.numeric(percentile))
  
  # Step 3: Format wide for arrow plot
  arrow_df <- plot_df %>%
    dplyr::select(all_of(c(subject_col, treatment_col, clinical_var_chr, "percentile", visit_col, "pseudotime", visit_treatment_col))) %>%
    pivot_wider(names_from = .data[[visit_col]],
                values_from = c(pseudotime, !!clinical_var, !!visit_treatment_col)) %>%
    dplyr::mutate(expectations = case_when(
      pseudotime_POST - pseudotime_PRE > 0 & !!sym(paste0(clinical_var_chr, "_POST")) - !!sym(paste0(clinical_var_chr, "_PRE")) > 0 ~ "Consistent",
      pseudotime_POST - pseudotime_PRE < 0 & !!sym(paste0(clinical_var_chr, "_POST")) - !!sym(paste0(clinical_var_chr, "_PRE")) < 0 ~ "Consistent",
      TRUE ~ "Inconsistent"
    ))
  
  # Step 4: Plot
  p <- arrow_df %>%
    filter(percentile == percentile_filter) %>%
    ggplot() +
    geom_segment(aes(x = pseudotime_PRE,
                     y = !!sym(paste0(clinical_var_chr, "_PRE")),
                     xend = pseudotime_POST,
                     yend = !!sym(paste0(clinical_var_chr, "_POST")),
                     color = expectations),
                 arrow = arrow(length = unit(0.5, "cm")),
                 alpha = 1) +
    geom_point(aes(x = pseudotime_PRE,
                   y = !!sym(paste0(clinical_var_chr, "_PRE")),
                   color = .data[[paste0(visit_treatment_col, "_PRE")]]),
               shape = shape_pre, size = 4, alpha = 0.7) +
    geom_point(aes(x = pseudotime_POST,
                   y = !!sym(paste0(clinical_var_chr, "_POST")),
                   color = .data[[paste0(visit_treatment_col, "_POST")]]),
               shape = shape_post, size = 4, alpha = 0.7) +
    scale_color_manual(values = color_palette) +
    labs(x = "Pseudotime",
         y = clinical_var_label %||% clinical_var_chr,
         color = NULL) +
    theme_minimal(base_size = 15) +
    theme(panel.grid = element_blank())
  
  if (!is.null(bucket)) {
    temp_file <- tempfile(fileext = ".png") # need to create a temporary file
    ggsave(filename = temp_file, width = 7, height = 5)
    s3$upload_file(temp_file, bucket, paste0("slingshot/attempt_arrow_scatter", 
                                             tolower(celltype_suffix), "_", clinical_var, "_slingshot.png"))
  }
  return(p)
}

# ===========================================================================
# Function: plot_clinvar_pseudotime_arrows
# ===========================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(rlang)
library(grid)

plot_clinvar_pseudotime_arrows_delta <- function(df,
                                                 pseudotime_var = "slingPseudotime_1",
                                                 subject_col = "subject_id",
                                                 visit_col = "visit",
                                                 treatment_col = "treatment",
                                                 visit_treatment_col = "visit_treatment",
                                                 clinical_var,                # unquoted, e.g. hba1c
                                                 clinical_var_label = NULL,   # for axis label
                                                 percentile_filter = 50,
                                                 percentiles = c(25, 50, 65, 85),
                                                 shape_pre = 16,
                                                 shape_post = 17,
                                                 bucket = "attempt",
                                                 celltype_suffix = "",
                                                 color_palette = c("PRE Placebo" = "#fbc4ab",
                                                                   "POST Placebo" = "#f4978e",
                                                                   "PRE Dapagliflozin" = "#ccd5ae",
                                                                   "POST Dapagliflozin" = "#828e82",
                                                                   "Healthy state" = "#8d99ae",
                                                                   "Injured state" = "#e5e5e5")) {
  
  clinical_var <- rlang::ensym(clinical_var)
  pseudotime_var_chr <- rlang::as_string(rlang::ensym(pseudotime_var))
  clinical_var_chr <- rlang::as_string(clinical_var)
  
  # Step 1: Pseudotime percentiles
  med <- df %>%
    group_by(.data[[subject_col]], .data[[visit_col]]) %>%
    dplyr::summarise(across(
      .data[[pseudotime_var_chr]],
      list("25" = ~quantile(., 0.25, na.rm = TRUE),
           "50" = ~quantile(., 0.50, na.rm = TRUE),
           "65" = ~quantile(., 0.65, na.rm = TRUE),
           "85" = ~quantile(., 0.85, na.rm = TRUE)),
      .names = "pseudotime_{fn}"
    ), .groups = "drop")
  
  # Step 2: Construct plot_df
  plot_df <- df %>%
    dplyr::select(all_of(c(subject_col, visit_col, treatment_col, clinical_var_chr, visit_treatment_col))) %>%
    distinct(.data[[subject_col]], .data[[visit_col]], .keep_all = TRUE) %>%
    left_join(med, by = c(subject_col, visit_col)) %>%
    pivot_longer(
      cols = starts_with("pseudotime_"),
      names_to = "percentile",
      names_prefix = "pseudotime_",
      values_to = "pseudotime"
    ) %>%
    dplyr::mutate(percentile = as.numeric(percentile))
  
  # Step 3: Format wide for arrow plot
  arrow_df <- plot_df %>%
    dplyr::select(all_of(c(subject_col, treatment_col, clinical_var_chr, "percentile", visit_col, "pseudotime", visit_treatment_col))) %>%
    pivot_wider(names_from = .data[[visit_col]],
                values_from = c(pseudotime, !!clinical_var, !!visit_treatment_col)) %>%
    dplyr::mutate(expectations = case_when(
      pseudotime_POST - pseudotime_PRE < 0 ~ "Healthy state",
      TRUE ~ "Injured state"
    ))
  
  # Step 4: Plot
  p <- arrow_df %>%
    filter(percentile == percentile_filter & 
             !is.na(expectations) & 
             (!is.na(.data[[paste0(visit_treatment_col, "_PRE")]])| !is.na(.data[[paste0(visit_treatment_col, "_POST")]]))) %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#2b2d42") +
    geom_segment(aes(x = pseudotime_PRE,
                     y = !!sym(paste0(clinical_var_chr, "_PRE")),
                     xend = pseudotime_POST,
                     yend = !!sym(paste0(clinical_var_chr, "_POST")),
                     color = expectations),
                 arrow = arrow(length = unit(0.5, "cm")),
                 alpha = 1) +
    geom_point(aes(x = pseudotime_PRE,
                   y = !!sym(paste0(clinical_var_chr, "_PRE")),
                   color = .data[[paste0(visit_treatment_col, "_PRE")]]),
               shape = shape_pre, size = 4, alpha = 0.7) +
    geom_point(aes(x = pseudotime_POST,
                   y = !!sym(paste0(clinical_var_chr, "_POST")),
                   color = .data[[paste0(visit_treatment_col, "_POST")]]),
               shape = shape_post, size = 4, alpha = 0.7) +
    scale_color_manual(values = color_palette) +
    labs(x = "Pseudotime",
         y = clinical_var_label %||% clinical_var_chr,
         color = NULL) +
    theme_minimal(base_size = 15) +
    theme(panel.grid = element_blank())
  
  if (!is.null(bucket)) {
    temp_file <- tempfile(fileext = ".png") # need to create a temporary file
    ggsave(filename = temp_file, width = 7, height = 5)
    s3$upload_file(temp_file, bucket, paste0("slingshot/attempt_arrow_scatter", 
                                             tolower(celltype_suffix), "_", clinical_var, "_slingshot.png"))
  }
  
  return(p)
}



# ===========================================================================
# Function: create_croc_attempt_df
# ===========================================================================

create_croc_attempt_df <- function(attempt_df,
                                   croc_df,
                                   cell_type = "PT",
                                   FC_attempt = "logFC_treatmentDapagliflozin:visitPOST",
                                   FC_croc = "logFC_groupType_1_Diabetes",
                                   p_attempt = "fdr",
                                   p_croc = "p_groupType_1_Diabetes",
                                   attempt_p_cut = 0.05,
                                   croc_p_cut = 0.05,
                                   save_csv = FALSE,
                                   csv_path = NULL) {
  
  # Filter extreme values from croc_df
  croc_df <- croc_df %>% 
    filter(abs(.data[[FC_croc]]) < 10)
  
  # Create the joined dataframe with direction classification
  joined_df <- attempt_df %>%
    dplyr::select(Gene, 
                  logFC_attempt = .data[[FC_attempt]], 
                  p.val_attempt = .data[[p_attempt]]) %>%
    full_join(
      croc_df %>%
        dplyr::select(Gene, 
                      logFC_croc = .data[[FC_croc]], 
                      p.val_croc = .data[[p_croc]]),
      by = "Gene"
    ) %>%
    mutate(
      direction = case_when(
        logFC_attempt < 0 & logFC_croc < 0 & p.val_attempt < attempt_p_cut & p.val_croc < croc_p_cut ~ "non-reversed",
        logFC_attempt > 0 & logFC_croc > 0 & p.val_attempt < attempt_p_cut & p.val_croc < croc_p_cut ~ "non-reversed",
        logFC_attempt < 0 & logFC_croc > 0 & p.val_attempt < attempt_p_cut & p.val_croc < croc_p_cut ~ "reversed",
        logFC_attempt > 0 & logFC_croc < 0 & p.val_attempt < attempt_p_cut & p.val_croc < croc_p_cut ~ "reversed",
        TRUE ~ NA_character_
      ),
      cell_type = cell_type
    ) %>%
    filter(!is.na(direction))
  
  # Save to CSV if requested
  if (save_csv && !is.null(csv_path)) {
    write.csv(joined_df, csv_path, row.names = FALSE)
    message("Saved ", nrow(joined_df), " overlapping genes to: ", csv_path)
  }
  
  return(joined_df)
}

# ===========================================================================
# Function: plot_croc_attempt_volcano
# ===========================================================================

plot_croc_attempt_volcano <- function(attempt_df,
                                      croc_df,
                                      cell_type = "PT",
                                      FC_attempt = "logFC_treatmentDapagliflozin:visitPOST",
                                      FC_croc = "logFC_groupType_1_Diabetes",
                                      p_attempt = "fdr",
                                      p_croc = "p_groupType_1_Diabetes",
                                      top_n = 20,
                                      attempt_p_cut = 0.05,
                                      croc_p_cut = 0.05,
                                      width = 7,
                                      height = 5,
                                      caption = T, 
                                      positive_text = "Upregulated in T1D compared to HC",
                                      negative_text = "Downregulated in T1D compared to HC",
                                      x_axis = "logFC_groupT1D",
                                      y_axis = "-log10(p-value)",
                                      save_path = NULL,
                                      return_df = FALSE) {
  
  # Create the joined dataframe
  joined_df <- create_croc_attempt_df(
    attempt_df = attempt_df,
    croc_df = croc_df,
    cell_type = cell_type,
    FC_attempt = FC_attempt,
    FC_croc = FC_croc,
    p_attempt = p_attempt,
    p_croc = p_croc,
    attempt_p_cut = attempt_p_cut,
    croc_p_cut = croc_p_cut
  )
  
  # Continue with the original plotting code...
  croc_df <- croc_df %>% 
    filter(abs(.data[[FC_croc]]) < 10)
  
  # Subset by significance and direction
  attempt_pos <- attempt_df %>% filter(.data[[FC_attempt]] > 0 & .data[[p_attempt]] < attempt_p_cut)
  attempt_neg <- attempt_df %>% filter(.data[[FC_attempt]] < 0 & .data[[p_attempt]] < attempt_p_cut)
  croc_pos <- croc_df %>% filter(.data[[FC_croc]] > 0 & .data[[p_croc]] < croc_p_cut)
  croc_neg <- croc_df %>% filter(.data[[FC_croc]] < 0 & .data[[p_croc]] < croc_p_cut)
  
  # Matching
  all_pos_match <- croc_pos %>% filter(Gene %in% attempt_pos$Gene)
  all_neg_match <- croc_neg %>% filter(Gene %in% attempt_neg$Gene)
  
  # Top genes for labeling
  top_pos <- croc_pos %>% dplyr::arrange(.data[[p_croc]]) %>% slice_head(n = top_n)
  top_neg <- croc_neg %>% dplyr::arrange(.data[[p_croc]]) %>% slice_head(n = top_n)
  
  top_merged <- rbind(top_pos, top_neg) %>%
    dplyr::mutate(attempt_direction = case_when(Gene %in% attempt_neg$Gene ~ "-",
                                                Gene %in% attempt_pos$Gene ~ "+",
                                                T ~ "NA"),
                  croc_direction = case_when(Gene %in% croc_neg$Gene ~ "-",
                                             Gene %in% croc_pos$Gene ~ "+"),
                  match = case_when(attempt_direction == croc_direction  ~ "nonreversed",
                                    (attempt_direction == "+" & croc_direction == "-") | 
                                      (attempt_direction == "-" & croc_direction == "+") ~ "reversed",
                                    T ~ "no match"),
                  top_color = case_when(Gene %in% croc_pos$Gene ~ "#f28482",
                                        Gene %in% croc_neg$Gene ~ "#457b9d",
                                        TRUE ~ "#ced4da"),
                  reversed_fill = case_when(attempt_direction == "+" & match == "reversed" ~ "+",
                                            attempt_direction == "-" & match == "reversed" ~ "-",
                                            T ~ "nonreversed"),
                  face = case_when(match == "reversed" ~ "bold",
                                   T ~ "plain"))
  
  n_reversed <- length(intersect(attempt_pos$Gene, croc_neg$Gene)) +
    length(intersect(attempt_neg$Gene, croc_pos$Gene))
  n_nonreversed <- nrow(all_pos_match) + nrow(all_neg_match)
  
  # Max and min for annotation arrows
  max_fc <- max(croc_df[[FC_croc]], na.rm = TRUE)
  min_fc <- min(croc_df[[FC_croc]], na.rm = TRUE)
  
  # Get y-axis max for dynamic scaling
  y_max <- max(-log10(croc_df[[p_croc]]), na.rm = TRUE) * 1.1
  
  croc_df <- croc_df %>%
    dplyr::mutate(top_color = case_when(Gene %in% croc_pos$Gene ~ "#f28482",
                                        Gene %in% croc_neg$Gene ~ "#457b9d",
                                        TRUE ~ "#ced4da"))
  
  # Create plot
  p <- ggplot(croc_df, aes(x = .data[[FC_croc]], y = -log10(.data[[p_croc]]))) +
    geom_hline(yintercept = -log10(attempt_p_cut), linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
    geom_point(alpha = 0.3, aes(color = top_color)) +
    geom_label_repel(data = top_merged,
                     aes(label = Gene, fill = reversed_fill, color = top_color, fontface = face),
                     size = 3, max.overlaps = Inf,
                     force = 30, segment.alpha = 0.5, segment.size = 0.4,
                     min.segment.length = 0,
                     segment.color = "#ced4da", seed = 1234, label.size = 0) +
    labs(x = x_axis,
         y = y_axis,
         caption = ifelse(caption, 
                          paste0("\n\nCell type: ", cell_type, " | ",
                                 "N of overlap: ", n_nonreversed + n_reversed, " (",
                                 "reversal: ", n_reversed, ", non-reversal: ", n_nonreversed, ")",
                                 "\n\nUp to top ", top_n, " from each direction are labeled.\n",
                                 "Reversed genes are boxed in fill color corresponding to direction of dapagliflozin effect."),
                          paste0("\n\nCell type: ", cell_type, " | ",
                                 "N of overlap: ", n_nonreversed + n_reversed, " (",
                                 "reversal: ", n_reversed, ", non-reversal: ", n_nonreversed, ")")
         )) +
    scale_size_continuous(range = c(1, 1.3)) +
    scale_fill_manual(values = c("+" = scales::alpha("#ff9996", 0.2),
                                 "-" = scales::alpha("#5c9fc1", 0.2),
                                 "nonreversed" = scales::alpha("white", 0))) +
    scale_color_manual(values = c("#457b9d" = "#457b9d",
                                  "#ced4da" = "#ced4da",
                                  "#f28482" = "#f28482")) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 9, family = "Arial"),
          plot.caption = element_text(size = 8.5, 
                                      hjust = 0.5, margin = margin(t = 8), 
                                      family = "Arial")) +
    guides(color = "none", size = "none", fill = "none") +
    annotate("segment", 
             x=max_fc/8, 
             xend=(max_fc*7)/8, 
             y=-y_max * 0.1,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(max_fc/8, (max_fc*7)/8)), 
             y=-y_max * 0.14, 
             label=positive_text,
             size=3, color="#343a40") +
    annotate("segment", 
             x=min_fc/8, 
             xend=(min_fc*7)/8, 
             y=-y_max * 0.1,
             col="darkgrey", arrow=arrow(length=unit(0.2, "cm"))) +
    annotate("text", 
             x=mean(c(min_fc/8, (min_fc*7)/8)), 
             y=-y_max * 0.14, 
             label=negative_text,
             size=3, color="#343a40") +
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(ylim = c(0, y_max), clip="off")
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = width, height = height)
  }
  
  if (return_df) {
    return(list(plot = p, data = joined_df))
  } else {
    return(p)
  }
}
# ===========================================================================
# Function: soma_corr
# ===========================================================================

soma_corr <- function(
    clinical_raw,          # raw clinical data (long format)
    soma_raw,              # raw soma data (long format)
    did_df,                # DiD results
    analytes_info,         # metadata (AptName + Target)
    clinical_var,          # base clinical variable name (e.g., "mgfr_si")
    feature_suffix,        # suffix for aptamer features
    top_n = 20,         # number of top labels to show
    proper_name,
    pre_visit = 0,
    post_visit = 16,
    label = "Target"
) {
  library(psych)
  library(dplyr)
  library(stringr)
  library(tibble)
  
  # --- Helper for sequence detection ---
  is_seq <- function(.x, feature_prefix = "seq", feature_suffix_local = feature_suffix) {
    # Build regex pattern depending on suffix
    pattern <- if (feature_suffix_local == "") {
      paste0("^", feature_prefix, ".*[0-9]+$")
    } else {
      paste0("^", feature_prefix, ".*", feature_suffix_local, "$")
    }
    
    grepl(pattern, .x)
  }
  
  # --- Filter relevant visits ---
  clinical <- clinical_raw %>%
    filter(visit %in% c(pre_visit, post_visit))
  clinical$record_id <- clinical$subject_id
  
  # --- Reshape clinical data wide ---
  keep_vars <- c("record_id", "visit", clinical_var)
  clinical_wide <- clinical[, keep_vars] %>%
    pivot_wider(names_from = visit,
                values_from = clinical_var,
                names_prefix = paste0(clinical_var, "_"))
  
  # dynamically compute delta
  col_pre  <- paste0(clinical_var, "_", pre_visit)
  col_post <- paste0(clinical_var, "_", post_visit)
  clinical_wide[[paste0("delta_", clinical_var)]] <- 
    clinical_wide[[col_post]] - clinical_wide[[col_pre]]
  
  # --- Reshape SOMA data wide and compute deltas ---
  seq_features <- names(soma_raw)[is_seq(names(soma_raw))]
  soma_keep <- soma_raw %>%
    dplyr::select(record_id, visit, all_of(seq_features))
  
  soma_wide <- reshape(soma_keep,
                       timevar = "visit",
                       idvar   = "record_id",
                       direction = "wide")
  
  for (feat in seq_features) {
    soma_wide[[paste0("delta_", feat)]] <-
      soma_wide[[paste0(feat, ".POST")]] - soma_wide[[paste0(feat, ".PRE")]]
  }
  
  # ---- filter DiD up/down ----
  did_up <- did_df %>% filter(logFC > 0 & adj.P.Val < 0.05)
  did_down <- did_df %>% filter(logFC < 0 & adj.P.Val < 0.05)
  
  # ---- correlation matrix ----
  x <- clinical_wide %>% 
    filter(record_id %in% soma_wide$record_id) %>%
    dplyr::select(all_of(paste0("delta_", clinical_var)))
  y <- soma_wide %>% dplyr::select(starts_with("delta_"))
  
  corr_res <- corr.test(x = x, y = y, method = "spearman", adjust = "none")
  corr_mat <- data.frame(
    spearman = t(corr_res$r),
    p.value  = t(corr_res$p.adj)
  )
  names(corr_mat)[1] <- paste0("spearman_delta_", clinical_var)
  names(corr_mat)[2] <- "p.value"
  
  corr_mat$AptName <- rownames(corr_mat) %>% str_remove("^delta_")
  
  corr_mat <- corr_mat %>%
    left_join(analytes_info %>% 
                dplyr::mutate(AptName = paste0(AptName, feature_suffix)), by = "AptName") %>%
    dplyr::mutate(Target_apt = paste0(Target, " (", AptName, ")")) %>%
    column_to_rownames("Target_apt") %>%
    dplyr::arrange(p.value)
  
  # ---- subset for DiD up/down ----
  corr_did_up <- corr_mat %>% filter(AptName %in% did_up$AptName)
  corr_did_down <- corr_mat %>% filter(AptName %in% did_down$AptName)
  
  # ---- volcano plot ----
  p <-  soma_plot_corr_volcano(
    data = corr_mat,
    FC = paste0("spearman_delta_", clinical_var),
    p.value = "p.value",
    labels = label,
    pos_did = corr_did_up,
    neg_did = corr_did_down,
    cutoff = 0.05,
    top_n = top_n,
    x.lab = paste0(proper_name, " Spearman Correlation"))
  
  
  list(
    correlation_table = corr_mat,
    plot_all   = p
  )
}

# ===========================================================================
# Function: soma_plot_corr_volcano
# ===========================================================================

soma_plot_corr_volcano <- function (data, 
                                    FC = "FC",           # Column name as string
                                    p.value = "p.value", # Column name as string
                                    labels = "labels",   # Column name as string
                                    identify = T, 
                                    identify_manual = NULL, 
                                    pt.size = 3, 
                                    text.size = 4, 
                                    cutoff = 0.05/nrow(data), 
                                    sig_pos_lab = "Significant (+)",
                                    sig_neg_lab = "Significant (-)", 
                                    ns_lab = "Non-Significant", 
                                    sig_pos_col = "#f28482", 
                                    sig_neg_col = "#457b9d", 
                                    ns_col = "#dad7cd",
                                    main = NULL, 
                                    x.lab = NULL, 
                                    overlaps = 10, 
                                    top_n = 20,
                                    pos_did,
                                    neg_did,
                                    ...) 
{
  # Convert column names to symbols for use with ggplot
  .fc_sym <- sym(FC)
  .p_sym <- sym(p.value)
  .label_sym <- sym(labels)
  
  # Check if FC values are log-transformed
  if (all(data[[FC]] >= 0)) {
    warning("It appears you are not passing log2-transformed ", 
            "fold-change values. Please check.", call. = FALSE)
  }
  
  if (is.null(main)) {
    main <- "Volcano Plot of Significant Changes"
  }
  if (is.null(x.lab)) {
    x.lab <- bquote(italic(log)[2] ~ (Fold ~ Change))
  }
  y.lab <- bquote(-italic(log)[10] ~ (p - value))
  
  # Create plot dataframe with simplified categorization
  plot_df <- data %>%
    dplyr::mutate(
      # Main categorization based on significance and direction
      significance_status = case_when(
        -log10(.data[[p.value]]) >= -log10(cutoff) & .data[[FC]] > 0 ~ sig_pos_lab,
        -log10(.data[[p.value]]) >= -log10(cutoff) & .data[[FC]] < 0 ~ sig_neg_lab,
        TRUE ~ ns_lab
      ),
      # Track if it's in DiD models
      in_pos_did = AptName %in% pos_did$AptName,
      in_neg_did = AptName %in% neg_did$AptName,
      DiD_status = case_when(in_pos_did | in_neg_did ~ "DiD",
                             T ~ "not-DiD"),
      # Create label for DiD entries
      did_label = case_when(
        in_pos_did ~ "DiD +",
        in_neg_did ~ "DiD -",
        T ~ "not-DiD"
      ),
      
      
      # Flag for significant entries
      is_significant = -log10(.data[[p.value]]) >= -log10(cutoff),
      sig_face = case_when(is_significant ~ "bold",
                           T ~ "plain"),
      did_face = case_when(DiD_status == "DiD" ~ "plain",
                           T ~ "bold"),
      combined_status = case_when(!is.na(significance_status) & !is.na(did_label) ~ paste0(significance_status," & ", did_label),
                                  !is.na(significance_status) ~ significance_status,
                                  !is.na(did_label) ~ did_label))
  
  # Set up colors
  sig_cols <- setNames(
    c(sig_pos_col, sig_neg_col, ns_col),
    c(sig_pos_lab, sig_neg_lab, ns_lab)
  )
  
  combined_cols <- setNames(
    c(sig_pos_col, sig_pos_col, ns_col,
      sig_neg_col, sig_neg_col,
      sig_pos_col, sig_neg_col,
      sig_pos_col, sig_neg_col),
    c(paste0(sig_pos_lab, " & ", "DiD (+)"), paste0(sig_neg_lab, " & ", "DiD (+)"), ns_lab,
      paste0(sig_pos_lab, " & ", "DiD (-)"), paste0(sig_neg_lab, " & ", "DiD (-)"),
      "DiD (+)", "DiD (-)",
      sig_pos_lab, sig_neg_lab)
  )
  
  did_cols <- setNames(
    c(fill_alpha(sig_pos_col, 0.3), fill_alpha(sig_neg_col, 0.3), fill_alpha("white", 0)),
    c("DiD +", "DiD -", "not-DiD")
  )
  plot_df$did_label = factor(plot_df$did_label, levels = c("DiD +", "DiD -", "not-DiD"))
  
  n_sig = nrow(plot_df %>% filter(is_significant))
  n_pos_sig = nrow(plot_df %>% filter(significance_status == "Significant (+)"))
  n_neg_sig = nrow(plot_df %>% filter(significance_status == "Significant (-)"))
  n_pos_did = nrow(plot_df %>% filter(is_significant & did_label == "DiD +"))
  n_neg_did = nrow(plot_df %>% filter(is_significant & did_label == "DiD -"))
  # Create base plot
  p <- ggplot(plot_df, aes(x = !!.fc_sym, y = -log10(!!.p_sym), color = significance_status)) + 
    geom_hline(yintercept = -log10(cutoff), 
               color = "grey50", linetype = "dashed") + 
    theme_minimal() + 
    theme(panel.grid = element_blank(),
          legend.position = c(0.5, 0.9),
          legend.direction = "horizontal",
          plot.caption = element_text(hjust = 0.5, size = 10)) + 
    geom_point(alpha = 0.5, size = pt.size) + 
    geom_label_repel(data = plot_df %>% 
                       dplyr::arrange(p.value) %>%
                       head(top_n) %>%
                       dplyr::filter(is_significant | DiD_status == "DiD"), 
                     aes(label = !!.label_sym,
                         fill = did_label,
                         fontface = sig_face),  
                     label.size = 0,
                     segment.color = "grey",
                     size = text.size, 
                     color = "black", 
                     max.overlaps = Inf) +
    scale_fill_manual(values = did_cols, name = "") +
    scale_color_manual(values = sig_cols, name = "") + 
    guides(color = "none") +
    labs(x = x.lab, y = y.lab,
         caption = paste0("\nSignificant findings:\n",
                          "Correlations n: ", n_sig,
                          " | Positive correlations n: ", n_pos_sig,
                          " | Negative correlations n: ", n_neg_sig,
                          "\nCorrelations & positive DiD n: ", n_pos_did,
                          " | Correlations & negative DiD n: ", n_neg_did,
                          "\n\nUp to top ", top_n, " proteins are labeled."))
  p
  
  return(p)
}

# ===========================================================================
# Function: run_limma_proteomics
# ===========================================================================


# Proteomics
library(dplyr)
library(tibble)
library(tidyr)
library(limma)

run_limma_proteomics <- function(
    data, 
    analyte_info, 
    treatment = "Placebo",
    visit = "PRE",
    visit_var = "visit",
    treatment_var = "treatment_arm",
    feature_prefix = "seq",
    feature_suffix = "",
    creatinine = NULL,
    model_type = c("within_treatment", "within_visit", "interaction", "interaction_random")
) {
  model_type <- match.arg(model_type)
  
  # --- Ensure factors have consistent levels ---
  data[[treatment_var]] <- factor(data[[treatment_var]], levels = c("Placebo", "Dapagliflozin"))
  data[[visit_var]] <- factor(data[[visit_var]], levels = c("PRE", "POST"))
  
  # --- Conditional filtering ---
  data_sub <- switch(model_type,
                     within_treatment = dplyr::filter(data, !!sym(treatment_var) == treatment),
                     within_visit     = dplyr::filter(data, !!sym(visit_var) == visit),
                     interaction      = data,
                     interaction_random = data
  ) %>%
    dplyr::arrange(record_id)
  
  # --- Check if creatinine variable exists when specified ---
  if (!is.null(creatinine) && !creatinine %in% names(data_sub)) {
    stop(paste("Creatinine variable", creatinine, "not found in data"))
  }
  
  # --- Select proteomics features only ---
  if (feature_suffix == "") {
    pattern <- paste0("^", feature_prefix, ".*[0-9]+$")
  } else {
    pattern <- paste0("^", feature_prefix, ".*", feature_suffix, "$")
  }
  
  features <- data_sub %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "record_id_visit") %>%
    dplyr::select(dplyr::matches(pattern))
  
  # --- Log2 transform ---
  y <- log2(t(features))
  
  if (model_type == "within_treatment") {
    # Build formula with or without creatinine
    if (!is.null(creatinine)) {
      design_formula <- as.formula(paste0("~0 + ", visit_var, " + ", creatinine))
    } else {
      design_formula <- as.formula(paste0("~0 + ", visit_var))
    }
    
    design_mat <- model.matrix(design_formula, data = data_sub)
    
    # Clean column names and create contrast
    colnames(design_mat) <- gsub(visit_var, "", colnames(design_mat))
    contrast <- limma::makeContrasts("POST - PRE", levels = design_mat)
    coef_name <- 1
    
    fit <- limma::lmFit(y, design_mat)
    fit <- limma::contrasts.fit(fit, contrasts = contrast)
    fit <- limma::eBayes(fit)
    res <- limma::topTable(fit, coef = coef_name, number = nrow(y), sort.by = "p")
    res <- res %>%
      rownames_to_column("feature")
    
  } else if (model_type == "within_visit") {
    # Build formula with or without creatinine
    if (!is.null(creatinine)) {
      design_formula <- as.formula(paste0("~0 + ", treatment_var, " + ", creatinine))
    } else {
      design_formula <- as.formula(paste0("~0 + ", treatment_var))
    }
    
    design_mat <- model.matrix(design_formula, data = data_sub)
    
    # Clean column names and create contrast
    colnames(design_mat) <- gsub(treatment_var, "", colnames(design_mat))
    contrast <- limma::makeContrasts("Dapagliflozin - Placebo", levels = design_mat)
    coef_name <- 1
    
    fit <- limma::lmFit(y, design_mat)
    fit <- limma::contrasts.fit(fit, contrasts = contrast)
    fit <- limma::eBayes(fit)
    res <- limma::topTable(fit, coef = coef_name, number = nrow(y), sort.by = "p")
    res <- res %>%
      rownames_to_column("feature")
    
  } else if (model_type == "interaction") {
    # Build formula with or without creatinine
    if (!is.null(creatinine)) {
      design_formula <- as.formula(paste0("~ ", treatment_var, " * ", visit_var, " + ", creatinine))
    } else {
      design_formula <- as.formula(paste0("~ ", treatment_var, " * ", visit_var))
    }
    
    design_mat <- model.matrix(design_formula, data = data_sub)
    coef_name <- paste0(treatment_var, "Dapagliflozin:", visit_var, "POST")
    
    fit <- limma::lmFit(y, design_mat)
    fit <- limma::eBayes(fit)
    res <- limma::topTable(fit, coef = coef_name, number = nrow(y), sort.by = "p")
    res <- res %>%
      rownames_to_column("feature")
    
  } else if (model_type == "interaction_random") {
    # --- Mixed model for each analyte ---
    if (!"record_id" %in% names(data_sub)) {
      stop("For interaction_random, data must have a 'subject_id' column")
    }
    suppressPackageStartupMessages({
      library(lme4)
      library(lmerTest)
      library(broom.mixed)
    })
    
    coef_name <- paste0(treatment_var, "Dapagliflozin:", visit_var, "POST")
    
    singular_count <- 0  # initialize counter
    
    results_list <- lapply(colnames(features), function(feat) {
      df_tmp <- data_sub %>%
        dplyr::select(record_id, all_of(visit_var), all_of(treatment_var), all_of(feat), 
                      if (!is.null(creatinine)) all_of(creatinine) else NULL) %>%
        dplyr::rename(y = !!feat) %>%
        dplyr::mutate(y = log2(y))
      
      # Build formula with or without creatinine
      if (!is.null(creatinine)) {
        form <- as.formula(paste("y ~", treatment_var, "*", visit_var, "+", creatinine, "+ (1|record_id)"))
      } else {
        form <- as.formula(paste("y ~", treatment_var, "*", visit_var, "+ (1|record_id)"))
      }
      
      # Suppress warnings from singular fit
      suppressWarnings({
        fit <- lmer(form, data = df_tmp)
      })
      
      # Count if singular
      if (lme4::isSingular(fit, tol = 1e-4)) {
        singular_count <<- singular_count + 1
      }
      
      # Get all fixed effects for res_full
      full_results <- broom.mixed::tidy(fit, effects = "fixed") %>%
        dplyr::mutate(feature = feat)
      
      # Return a list with both filtered and full results
      list(
        filtered = full_results %>% dplyr::filter(term == coef_name),
        full = full_results
      )
    })
    
    # Extract filtered results for res
    res <- dplyr::bind_rows(lapply(results_list, function(x) x$filtered)) %>%
      dplyr::rename(logFC = estimate, P.Value = p.value) %>%
      dplyr::arrange(P.Value) %>%
      dplyr::mutate(adj.P.Val = p.adjust(P.Value, method = "fdr"))
    
    # Extract full results for res_full
    res_full <- dplyr::bind_rows(lapply(results_list, function(x) x$full)) %>%
      dplyr::rename(logFC = estimate, P.Value = p.value) %>%
      dplyr::arrange(feature, term) %>%
      dplyr::group_by(term) %>%
      dplyr::mutate(adj.P.Val = p.adjust(P.Value, method = "fdr")) %>%
      dplyr::ungroup()
    
    res_full <- res_full %>%
      dplyr::rename(AptName = feature) %>%
      dplyr::left_join(
        analyte_info %>% 
          dplyr::mutate(AptName = paste0(AptName, feature_suffix)),
        by = "AptName"
      )
    
    # Output singular count
    message("Number of singular fits: ", singular_count, " out of ", ncol(features),
            " (", round(100*singular_count/ncol(features), 1), "%)")
  }
  
  # --- Add annotation ---
  res_save <- res %>%
    dplyr::rename(AptName = feature) %>%
    dplyr::left_join(
      analyte_info %>% 
        dplyr::mutate(AptName = paste0(AptName, feature_suffix)),
      by = "AptName"
    )
  
  return(list(
    fit = if (model_type == "interaction_random") NULL else fit,
    results = res,
    results_full = if (model_type == "interaction_random") res_full else NULL, 
    results_annotated = res_save
  ))
}

# ===========================================================================
# Function: run_fgsea_analysis
# ===========================================================================

run_fgsea_analysis <- function(bg_path = file.path(root_path, "GSEA/"),
                               results_annotated,
                               stat_col = "t",
                               gene_col = "EntrezGeneSymbol",
                               minSize_kegg = 3,
                               maxSize_kegg = 500,
                               minSize_reactome = 3,
                               maxSize_reactome = 500,
                               minSize_go = 5,
                               maxSize_go = 500,
                               minSize_full = 5,
                               maxSize_full = 500,
                               minSize_hallmark = 5,
                               maxSize_hallmark = 500,
                               nPermSimple = 10000,
                               nproc = 1,
                               seed = 1234,
                               references = c("kegg_legacy", "reactome", "go", "full", "hallmark")) {
  
  # --- Prepare GMT files ---
  gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
  
  # Initialize pathway variables
  kegg_legacy <- NULL
  reactome <- NULL
  go <- NULL
  full <- NULL
  hallmark <- NULL
  
  # Only prepare GMT files that are in references
  if ("kegg_legacy" %in% references) {
    kegg_legacy <- prepare_gmt(gmt_files[1], unique(results_annotated[[gene_col]]), savefile = FALSE)
  }
  
  if ("reactome" %in% references) {
    reactome <- prepare_gmt(gmt_files[3], unique(results_annotated[[gene_col]]), savefile = FALSE)
  }
  
  if ("go" %in% references) {
    go <- prepare_gmt(gmt_files[4], unique(results_annotated[[gene_col]]), savefile = FALSE)
  }
  
  if ("full" %in% references) {
    full <- prepare_gmt(gmt_files[5], unique(results_annotated[[gene_col]]), savefile = FALSE)
  }
  
  if ("hallmark" %in% references) {
    hallmark <- prepare_gmt(gmt_files[6], unique(results_annotated[[gene_col]]), savefile = FALSE)
  }
  
  
  # --- Rank genes by specified statistic ---
  if (!stat_col %in% names(results_annotated)) {
    stop(paste("Column", stat_col, "not found in results_annotated"))
  }
  
  rankings <- results_annotated[[stat_col]]
  names(rankings) <- results_annotated[[gene_col]]
  rankings <- sort(rankings, decreasing = TRUE)
  
  # --- Run FGSEA ---
  set.seed(seed)
  
  # Initialize result variables as NULL
  kegg_res <- NULL
  reactome_res <- NULL
  go_res <- NULL
  full_res <- NULL
  hallmark_res <- NULL
  
  if ("kegg_legacy" %in% references) {
    kegg_res <- fgsea(pathways = kegg_legacy,
                      stats = rankings,
                      scoreType = 'std', 
                      minSize = minSize_kegg,
                      maxSize = maxSize_kegg,
                      nproc = nproc,
                      nPermSimple = nPermSimple)
  }
  
  if ("reactome" %in% references) {
    reactome_res <- fgsea(pathways = reactome,
                          stats = rankings,
                          scoreType = 'std', 
                          minSize = minSize_reactome,
                          maxSize = maxSize_reactome,
                          nproc = nproc,
                          nPermSimple = nPermSimple)
  }
  
  if ("go" %in% references) {
    go_res <- fgsea(pathways = go,
                    stats = rankings,
                    scoreType = "std",
                    minSize = minSize_go,
                    maxSize = maxSize_go,
                    nPermSimple = nPermSimple,
                    nproc = nproc)
  }
  
  if ("full" %in% references) {
    full_res <- fgsea(pathways = full,
                      stats = rankings,
                      scoreType = "std",
                      minSize = minSize_full,
                      maxSize = maxSize_full,
                      nPermSimple = nPermSimple,
                      nproc = nproc)
  }
  
  if ("hallmark" %in% references) {
    hallmark_res <- fgsea(pathways = hallmark,
                          stats = rankings,
                          scoreType = "std",
                          minSize = minSize_hallmark,
                          maxSize = maxSize_hallmark,
                          nPermSimple = nPermSimple,
                          nproc = nproc)
  }
  
  # --- Build summary dataframe dynamically ---
  summary_list <- list()
  
  if ("kegg_legacy" %in% references && !is.null(kegg_res)) {
    summary_list[["KEGG Legacy"]] <- c(
      sum(kegg_res$padj < 0.05, na.rm = TRUE),
      sum(kegg_res$pval < 0.05, na.rm = TRUE)
    )
  }
  
  if ("reactome" %in% references && !is.null(reactome_res)) {
    summary_list[["REACTOME"]] <- c(
      sum(reactome_res$padj < 0.05, na.rm = TRUE),
      sum(reactome_res$pval < 0.05, na.rm = TRUE)
    )
  }
  
  if ("go" %in% references && !is.null(go_res)) {
    summary_list[["GO"]] <- c(
      sum(go_res$padj < 0.05, na.rm = TRUE),
      sum(go_res$pval < 0.05, na.rm = TRUE)
    )
  }
  
  if ("full" %in% references && !is.null(full_res)) {
    summary_list[["FULL"]] <- c(
      sum(full_res$padj < 0.05, na.rm = TRUE),
      sum(full_res$pval < 0.05, na.rm = TRUE)
    )
  }
  
  if ("hallmark" %in% references && !is.null(hallmark_res)) {
    summary_list[["HALLMARK"]] <- c(
      sum(hallmark_res$padj < 0.05, na.rm = TRUE),
      sum(hallmark_res$pval < 0.05, na.rm = TRUE)
    )
  }
  
  # Convert list to dataframe
  if (length(summary_list) > 0) {
    summary_df <- as.data.frame(summary_list)
    rownames(summary_df) <- c("adj.pval", "p.val")
  } else {
    summary_df <- data.frame()
  }
  
  # --- Return results ---
  return(list(
    summary = summary_df,
    kegg = if("kegg_legacy" %in% references) kegg_res else NULL,
    reactome = if("reactome" %in% references) reactome_res else NULL,
    go = if("go" %in% references) go_res else NULL,
    full = if("full" %in% references) full_res else NULL,
    hallmark = if("hallmark" %in% references) hallmark_res else NULL
  ))
}

# ===========================================================================
# Function: process_enrichr_data
# ===========================================================================

# Function to process enrichr results
process_enrichr_data <- function(enrichr_result, response_type, top_n = 20) {
  # Extract the Reactome 2022 results (or whichever database you used)
  data <- enrichr_result[["Reactome_Pathways_2024"]] # Adjust database name as needed
  
  # Clean and process data
  data_clean <- data %>%
    filter(Adjusted.P.value < 0.05) %>%  # Filter significant pathways
    head(top_n) %>%
    dplyr::mutate(
      neg_log_p = -log10(P.value),
      neg_log_adj_p = -log10(Adjusted.P.value),
      gene_count = as.numeric(str_extract(Overlap, "\\d+")),
      response_type = response_type,
      # Create simplified pathway names for plotting
      pathway_short = str_trunc(Term, 40),
      # Assign categories based on pathway names
      category = case_when(
        str_detect(Term, "Immune|Neutrophil|Cytokine|Interferon|Interleukin") ~ "Immune",
        str_detect(Term, "Metabolism|Metabolic|Amino Acid|Fatty Acid|Glucose") ~ "Metabolism", 
        str_detect(Term, "Signal|Signaling|Receptor|Tyrosine|Growth Factor") ~ "Signaling",
        str_detect(Term, "Proteasome|Autophagy|Ubiquitin|Degradation|Quality") ~ "Quality Control",
        str_detect(Term, "Stress|Oxidative|Response|NFE2L2|KEAP1") ~ "Stress Response",
        str_detect(Term, "Mitochondrial|Respiratory|Electron|ATP") ~ "Mitochondrial",
        str_detect(Term, "Apoptosis|Cell Death|Programmed") ~ "Cell Death",
        str_detect(Term, "Transport|Trafficking|Vesicle|Membrane") ~ "Transport",
        TRUE ~ "Other"
      )
    )
  
  return(data_clean)
}

# ===========================================================================
# Function: prepare_pathway_data
# ===========================================================================

# Process enrichr results
prepare_pathway_data <- function(negative_paths_df, positive_paths_df, discordant_paths_df) {
  # Process each dataset
  negative_paths <- process_enrichr_data(negative_paths_df, "Negative")
  positive_paths <- process_enrichr_data(positive_paths_df, "Positive") 
  discordant_paths <- process_enrichr_data(discordant_paths_df, "Discordant")
  
  return(list(
    negative = negative_paths,
    positive = positive_paths,
    discordant = discordant_paths
  ))
}

# ===========================================================================
# Function: perform_concordance_analysis
# ===========================================================================
# Generalized Concordance Analysis Function
# Works for any cell type and protein source (urine/plasma)

library(ggplot2)
library(dplyr)
library(ggrepel)
library(enrichR)
library(stringr)

# Main function to perform concordance analysis
perform_concordance_analysis <- function(
    protein_data,             # The actual protein dataframe (not wrapped in results_annotated)
    transcript_data,          # Single-cell transcript data 
    celltype_name,            # Name for labeling (e.g., "PT", "LOH", etc.)
    protein_source = "urine", # "urine" or "plasma" 
    fdr_cutoff = 0.05,        # FDR cutoff for transcript significance
    adj_p_cutoff = 0.05,      # Adjusted p-value cutoff for proteins
    logfc_col = "logFC",      # Column name for protein logFC
    adj_p_col = "adj.P.Val",  # Column name for adjusted p-value
    gene_col = "EntrezGeneSymbol", # Column name for gene symbols
    transcript_logfc_col = "logFC_treatmentDapagliflozin:visitPOST", # Transcript logFC column
    save_plots = TRUE,        # Whether to save plots
    output_dir = ".",         # Output directory for plots
    databases = c("Reactome_Pathways_2024"), # Enrichr databases
    max_overlaps = Inf,       # Max overlaps for text labels
    volcano_width = 7,          # Plot width for saving
    volcano_height = 7,          # Plot height for saving
    bubble_width = 15,
    bubble_height = 7
) {
  
  # Validate inputs
  if (!protein_source %in% c("urine", "plasma")) {
    stop("protein_source must be either 'urine' or 'plasma'")
  }
  
  # Check if required columns exist
  required_protein_cols <- c(logfc_col, adj_p_col, gene_col)
  missing_cols <- required_protein_cols[!required_protein_cols %in% colnames(protein_data)]
  if (length(missing_cols) > 0) {
    stop("Missing columns in protein data: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check transcript data columns
  if (!transcript_logfc_col %in% colnames(transcript_data)) {
    cat("Available columns in transcript data:", paste(colnames(transcript_data), collapse = ", "), "\n")
    stop("Missing transcript logFC column: ", transcript_logfc_col)
  }
  
  # Add concordance classification if not already present
  if (!"conc_class" %in% colnames(protein_data)) {
    
    # Get significant transcript genes with direction
    sig_transcript_up <- subset(transcript_data, fdr < fdr_cutoff & transcript_data[transcript_logfc_col] > 0)$Gene
    sig_transcript_down <- subset(transcript_data, fdr < fdr_cutoff & transcript_data[transcript_logfc_col] < 0)$Gene
    
    protein_data <- protein_data %>%
      dplyr::mutate(
        label = case_when(
          !!sym(gene_col) %in% subset(transcript_data, fdr < fdr_cutoff)$Gene & 
            !!sym(adj_p_col) < adj_p_cutoff ~ !!sym(gene_col),
          TRUE ~ ""
        ),
        transcript_direction = case_when(
          !!sym(gene_col) %in% sig_transcript_up ~ "+",    # Only significant UP genes
          !!sym(gene_col) %in% sig_transcript_down ~ "-",  # Only significant DOWN genes  
          TRUE ~ "No transcript data"  # No significant transcript data
        ),
        # concordance relative to protein logFC
        conc_class = case_when(
          !!sym(logfc_col) >  0 & transcript_direction == "+" ~ "Concordant (+/+)",
          !!sym(logfc_col) <  0 & transcript_direction == "-" ~ "Concordant (-/-)",
          !!sym(logfc_col) >  0 & transcript_direction == "-" ~ "Discordant (+/-)",
          !!sym(logfc_col) <  0 & transcript_direction == "+" ~ "Discordant (-/+)",
          TRUE ~ "No transcript call"
        )
      )
  } else {
    # If conc_class exists, make sure we have label for plotting
    if (!"label" %in% colnames(protein_data)) {
      protein_data <- protein_data %>%
        dplyr::mutate(
          label = case_when(
            !!sym(gene_col) %in% subset(transcript_data, fdr < fdr_cutoff)$Gene & 
              !!sym(adj_p_col) < adj_p_cutoff ~ !!sym(gene_col),
            TRUE ~ ""
          )
        )
    }
    
    # Update transcript_direction if it exists but uses the old " " format or wrong logic
    if ("transcript_direction" %in% colnames(protein_data)) {
      
      # Get significant transcript genes with direction  
      sig_transcript_up <- subset(transcript_data, fdr < fdr_cutoff & transcript_data[transcript_logfc_col] > 0)$Gene
      sig_transcript_down <- subset(transcript_data, fdr < fdr_cutoff & transcript_data[transcript_logfc_col] < 0)$Gene
      
      protein_data <- protein_data %>%
        dplyr::mutate(
          transcript_direction = case_when(
            !!sym(gene_col) %in% sig_transcript_up ~ "+",
            !!sym(gene_col) %in% sig_transcript_down ~ "-", 
            TRUE ~ "No transcript data"
          ),
          # Recalculate concordance with corrected transcript direction
          conc_class = case_when(
            !!sym(logfc_col) >  0 & transcript_direction == "+" ~ "Concordant (+/+)",
            !!sym(logfc_col) <  0 & transcript_direction == "-" ~ "Concordant (-/-)",
            !!sym(logfc_col) >  0 & transcript_direction == "-" ~ "Discordant (+/-)",
            !!sym(logfc_col) <  0 & transcript_direction == "+" ~ "Discordant (-/+)",
            TRUE ~ "No transcript call"
          )
        )
    }
  }
  
  # Calculate counts
  counts <- protein_data %>%
    dplyr::summarise(
      n_pos_conc = sum(conc_class == "Concordant (+/+)", na.rm = TRUE),
      n_neg_conc = sum(conc_class == "Concordant (-/-)", na.rm = TRUE),
      n_disconc  = sum(conc_class %in% c("Discordant (+/-)", "Discordant (-/+)"), na.rm = TRUE),
      n_no_tx    = sum(conc_class == "No transcript call", na.rm = TRUE)
    )
  
  # Create caption text
  cap_txt <- sprintf(
    "Cell type: %s | Protein source: %s\n+ concordant: %d | - concordant: %d\ndiscordant: %d | no transcript call: %d",
    celltype_name, protein_source, counts$n_pos_conc, counts$n_neg_conc, 
    counts$n_disconc, counts$n_no_tx
  )
  
  # Create volcano plot with directional annotations
  # Max and min for annotation arrows
  max_fc <- max(protein_data[[logfc_col]], na.rm = TRUE)
  min_fc <- min(protein_data[[logfc_col]], na.rm = TRUE)
  
  # Get y-axis max for dynamic scaling
  y_max <- max(-log10(protein_data[[adj_p_col]]), na.rm = TRUE)
  
  # Annotation parameters
  arrow_padding <- 0.05  # How far below x-axis to put arrows
  arrow_text_padding <- 0.08  # How far below arrows to put text
  positive_text <- paste0("Higher in ", protein_source, " after SGLT2i")
  negative_text <- paste0("Lower in ", protein_source, " after SGLT2i")
  
  volcano_plot <- protein_data %>%
    ggplot(aes(x = !!sym(logfc_col), y = -log10(!!sym(adj_p_col)), color = transcript_direction)) +
    geom_hline(yintercept = -log10(adj_p_cutoff), linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
    geom_point(alpha = 0.4) +
    geom_text_repel(
      aes(label = label),
      size = 3, max.overlaps = max_overlaps, force = 30,
      segment.alpha = 0.5, segment.size = 0.4,
      min.segment.length = 0, segment.color = "#ced4da",
      seed = 1234, label.size = 0
    ) +
    annotate("segment", 
             x = max_fc/8, 
             xend = (max_fc*7)/8, 
             y = -y_max * arrow_padding,
             col = "darkgrey", arrow = arrow(length = unit(0.2, "cm"))) +
    annotate("text", 
             x = mean(c(max_fc/8, (max_fc*7)/8)), 
             y = -y_max * arrow_text_padding, 
             label = positive_text,
             size = 3, color = "#343a40") +
    annotate("segment", 
             x = min_fc/8, 
             xend = (min_fc*7)/8, 
             y = -y_max * arrow_padding,
             col = "darkgrey", arrow = arrow(length = unit(0.2, "cm"))) +
    annotate("text", 
             x = mean(c(min_fc/8, (min_fc*7)/8)), 
             y = -y_max * arrow_text_padding, 
             label = negative_text,
             size = 3, color = "#343a40") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = 10),
      title = element_text(size = 11),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.spacing.x = unit(0.05, 'cm'),
      plot.margin = margin(t = 10, r = 20, b = 25, l = 20),
      axis.title.x = element_text(margin = margin(t = 32)),
      plot.caption = element_text(size = 8.5, hjust = 0.5, margin = margin(t = 15)),
      legend.margin = margin(t = 5, b = 5),
      legend.title = element_blank()
    ) +
    scale_color_manual(
      values = c("+" = "#f28482", "-" = "#457b9d", "No transcript data" = "#edede9"),
      labels = c("+" = "Transcript ↑", "-" = "Transcript ↓", "No transcript data" = "No transcript data")
    ) +
    labs(x = paste0("logFC (", protein_source, " proteins)"), 
         y = "-log10(adj p-value)",
         caption = cap_txt)
  
  # Create concordant and discordant gene lists (STRATIFIED)
  neg_concordant <- protein_data %>%
    filter(!!sym(logfc_col) < 0) %>% pull(!!sym(gene_col))
  pos_concordant <- protein_data %>%
    filter(!!sym(logfc_col) > 0) %>% pull(!!sym(gene_col))
  
  transcript_neg <- transcript_data %>%
    filter(transcript_data[transcript_logfc_col] < 0) %>% 
    pull(Gene)
  transcript_pos <- transcript_data %>%
    filter(transcript_data[transcript_logfc_col] > 0) %>% 
    pull(Gene)
  
  # Create stratified concordant and discordant gene lists
  neg_concordant_genes <- neg_concordant[neg_concordant %in% transcript_neg]
  pos_concordant_genes <- pos_concordant[pos_concordant %in% transcript_pos]
  discordant_pos_neg <- pos_concordant[pos_concordant %in% transcript_neg]  # Protein ↑, Transcript ↓
  discordant_neg_pos <- neg_concordant[neg_concordant %in% transcript_pos]  # Protein ↓, Transcript ↑
  
  # Remove empty strings and NA values
  neg_concordant_genes <- neg_concordant_genes[neg_concordant_genes != "" & !is.na(neg_concordant_genes)]
  pos_concordant_genes <- pos_concordant_genes[pos_concordant_genes != "" & !is.na(pos_concordant_genes)]
  discordant_pos_neg <- discordant_pos_neg[discordant_pos_neg != "" & !is.na(discordant_pos_neg)]
  discordant_neg_pos <- discordant_neg_pos[discordant_neg_pos != "" & !is.na(discordant_neg_pos)]
  
  # Perform pathway enrichment analysis (STRATIFIED)
  enrichment_results <- list()
  
  if (length(neg_concordant_genes) > 2) {
    cat("Running enrichment for", length(neg_concordant_genes), "negative concordant genes...\n")
    enrichment_results$negative <- enrichr(neg_concordant_genes, databases)
  }
  if (length(pos_concordant_genes) > 2) {
    cat("Running enrichment for", length(pos_concordant_genes), "positive concordant genes...\n")
    enrichment_results$positive <- enrichr(pos_concordant_genes, databases)
  }
  if (length(discordant_pos_neg) > 2) {
    cat("Running enrichment for", length(discordant_pos_neg), "discordant (+/-) genes...\n")
    enrichment_results$discordant_damage <- enrichr(discordant_pos_neg, databases)
  }
  if (length(discordant_neg_pos) > 2) {
    cat("Running enrichment for", length(discordant_neg_pos), "discordant (-/+) genes...\n")
    enrichment_results$discordant_protective <- enrichr(discordant_neg_pos, databases)
  }
  
  # Create bubble plot if we have enrichment results
  bubble_plot <- NULL
  if (length(enrichment_results) > 0) {
    
    tryCatch({
      pathway_data <- prepare_pathway_data_generic(enrichment_results, databases[1])
      
      if (nrow(pathway_data) > 0) {
        bubble_plot <- pathway_data %>%
          ggplot(aes(x = Odds.Ratio, y = -log10(Adjusted.P.value), color = response_type)) +
          geom_point(aes(size = gene_count), alpha = 0.6) +
          geom_text_repel(
            aes(label = pathway_short), 
            max.overlaps = max_overlaps,
            force = 15, 
            segment.alpha = 0.3, 
            segment.size = 0.3,
            size = 3
          ) +
          theme_minimal() +
          scale_size_continuous(range = c(2, 8), name = "Gene Count") +
          theme(
            text = element_text(size = 10),
            panel.grid = element_blank(),
            legend.position = "right"
          ) +
          scale_x_log10() +
          scale_y_log10() +
          scale_color_manual(values = c(
            "Positive" = "#f28482",
            "Negative" = "#457b9d",
            "Discordant (+/-)" = "#adc178",     # Orange-red for damage
            "Discordant (-/+)" = "#3a5a40",     # Teal for protection
            "Discordant" = "#3a5a40"            # Keep old green for compatibility
          )) +
          labs(
            color = "Concordance",
            x = "Odds Ratio",
            y = "-log10(Adjusted P-value)",
            caption = cap_txt
          )
      }
    }, error = function(e) {
      warning("Could not create bubble plot: ", e$message)
    })
  }
  
  # Save plots if requested
  if (save_plots && !is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    volcano_filename <- file.path(output_dir, paste0(tolower(celltype_name), "_", 
                                                     protein_source, "_volcano.png"))
    ggsave(volcano_filename, volcano_plot, width = volcano_width, height = volcano_height)
    cat("Volcano plot saved:", volcano_filename, "\n")
    
    if (!is.null(bubble_plot)) {
      bubble_filename <- file.path(output_dir, paste0(tolower(celltype_name), "_", 
                                                      protein_source, "_pathways.png"))
      ggsave(bubble_filename, bubble_plot, width = bubble_width, height = bubble_height)
      cat("Pathway bubble plot saved:", bubble_filename, "\n")
    }
  }
  
  # Return results
  results <- list(
    counts = counts,
    gene_lists = list(
      negative_concordant = neg_concordant_genes,
      positive_concordant = pos_concordant_genes,
      discordant_damage = discordant_pos_neg,      # Protein ↑, Transcript ↓
      discordant_protective = discordant_neg_pos   # Protein ↓, Transcript ↑
    ),
    enrichment_results = enrichment_results,
    annotated_data = protein_data,
    plots = list(
      volcano = volcano_plot,
      bubble = bubble_plot
    ),
    summary = list(
      celltype = celltype_name,
      protein_source = protein_source,
      total_proteins = nrow(protein_data),
      significant_proteins = sum(protein_data[[adj_p_col]] < adj_p_cutoff, na.rm = TRUE),
      overlapping_genes = length(neg_concordant_genes) + length(pos_concordant_genes) + 
        length(discordant_pos_neg) + length(discordant_neg_pos)
    )
  )
  
  return(results)
}

# Helper function to process pathway data generically (UPDATED FOR STRATIFIED DISCORDANCE)
prepare_pathway_data_generic <- function(enrichment_results, database_name) {
  
  pathway_list <- list()
  
  for (response_type in names(enrichment_results)) {
    if (database_name %in% names(enrichment_results[[response_type]])) {
      data <- enrichment_results[[response_type]][[database_name]]
      
      if (nrow(data) > 0) {
        processed_data <- data %>%
          filter(Adjusted.P.value < 0.05) %>%
          head(20) %>%
          dplyr::mutate(
            response_type = case_when(
              response_type == "negative" ~ "Negative",
              response_type == "positive" ~ "Positive", 
              response_type == "discordant_damage" ~ "Discordant (+/-)",
              response_type == "discordant_protective" ~ "Discordant (-/+)",
              response_type == "discordant" ~ "Discordant"  # Keep old version for compatibility
            ),
            gene_count = as.numeric(stringr::str_extract(Overlap, "\\d+")),
            pathway_short = stringr::str_trunc(Term, 40)
          )
        
        pathway_list[[response_type]] <- processed_data
      }
    }
  }
  
  if (length(pathway_list) > 0) {
    return(do.call(rbind, pathway_list))
  } else {
    return(data.frame())
  }
}

# Batch analysis function for multiple cell types with single protein dataset
batch_concordance_analysis <- function(
    protein_data,             # Single protein dataframe (same urine/plasma for all cell types)
    transcript_data_list,     # Named list of transcript data (e.g., list(PT = pt_kpmp, LOH = loh_kpmp))
    protein_source = "urine", # "urine" or "plasma"
    output_dir = ".",         # Output directory
    ...                       # Additional arguments passed to perform_concordance_analysis
) {
  
  results_list <- list()
  
  for (celltype in names(transcript_data_list)) {
    cat("Processing", celltype, "vs", protein_source, "proteomics...\n")
    
    # Create cell-type specific output directory
    celltype_dir <- file.path(output_dir, celltype)
    if (!dir.exists(celltype_dir)) {
      dir.create(celltype_dir, recursive = TRUE)
    }
    
    results_list[[celltype]] <- perform_concordance_analysis(
      protein_data = protein_data,                          # Same protein data for all
      transcript_data = transcript_data_list[[celltype]],   # Different transcript data
      celltype_name = celltype,
      protein_source = protein_source,
      output_dir = celltype_dir,
      ...
    )
  }
  
  return(results_list)
}

# ===========================================================================
# Function: vertical_upset
# ===========================================================================

# install.packages(c("dplyr","tidyr","purrr","ggplot2","cowplot"))
# install.packages(c("dplyr","tidyr","purrr","ggplot2","cowplot"))
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(cowplot)

vertical_upset <- function(df, sets, top_n = Inf, min_size = 1,
                           order_sets = sets,
                           bar_fill = "grey60",
                           show_set_sizes = TRUE,
                           set_bar_fill = "grey50",
                           connect_lines = TRUE,
                           line_color = "grey60",
                           line_width = 0.6,
                           scale_type = "linear",  # "linear", "log10", "sqrt", or "break"
                           break_point = NULL,     # For scale_type = "break"
                           break_ratio = 0.5,      # Proportion of plot for lower range
                           set_colors = NULL,      # Custom colors for each set
                           multi_bar_fill = "#7A9A9F",  # Color for bars with 3+ intersections
                           multi_threshold = 3) {  # Threshold for coloring bars
  stopifnot(all(sets %in% names(df)))
  
  # Define muted color palette if not provided
  if (is.null(set_colors)) {
    # Color palette from the provided image
    muted_palette <- c(
      "#3B5A5B",  # Dark teal/grey
      "#5A8A80",  # Medium teal
      "#7ABAA2",  # Light teal
      "#A8B88A",  # Sage green
      "#E6C86E",  # Golden yellow
      "#E89B5C",  # Light orange
      "#D16558",  # Coral red
      "#8C6BB1"   # Muted plum
    )
    set_colors <- setNames(muted_palette[1:length(sets)], sets)
  } else {
    # Ensure all sets have colors
    if (!all(sets %in% names(set_colors))) {
      stop("set_colors must include all sets")
    }
  }
  
  # Ensure 0/1 membership matrix
  m <- df %>%
    dplyr::select(all_of(sets)) %>%
    dplyr::mutate(across(everything(), ~ as.integer(.x > 0)))
  
  # Map each row to a combination string (e.g., "A & C"), drop rows with no set membership
  row_combos <- m %>%
    dplyr::mutate(.row = dplyr::row_number()) %>%
    tidyr::pivot_longer(all_of(sets), names_to = "set", values_to = "present") %>%
    dplyr::group_by(.row) %>%
    dplyr::summarise(
      members = list(set[present == 1]), 
      n_members = sum(present == 1),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_members > 1) %>%  # Filter out single-set memberships BEFORE creating combinations
    dplyr::mutate(
      combination = purrr::map_chr(
        members,
        ~ paste(.x, collapse = " & ")
      )
    )
  
  # Intersection sizes (counts of identical combinations)
  combos <- row_combos %>%
    dplyr::count(combination, name = "size", sort = TRUE) %>%
    dplyr::filter(size >= min_size)
  
  if (!is.finite(top_n)) top_n <- nrow(combos)
  combos <- combos %>% dplyr::slice_head(n = top_n)
  
  # IMPORTANT: Ensure we only keep combinations with size > 0
  combos <- combos %>% dplyr::filter(size > 0)
  
  # Keep only rows belonging to the kept combinations (for consistent set sizes)
  rows_kept <- row_combos %>%
    dplyr::filter(combination %in% combos$combination) %>%
    dplyr::pull(.row)
  m_kept <- m[rows_kept, , drop = FALSE]
  
  # Set sizes among the kept rows (so top bars match the currently displayed intersections)
  set_sizes <- tibble::tibble(
    set = factor(sets, levels = order_sets),
    size = colSums(m_kept[, sets, drop = FALSE] > 0)
  )
  
  # Dot-matrix long form for kept combinations only
  mat_long <- tibble::tibble(combination = combos$combination) %>%
    tidyr::crossing(set = sets) %>%
    dplyr::left_join(
      combos %>%
        dplyr::mutate(members = strsplit(combination, " & ", fixed = TRUE)) %>%
        dplyr::select(combination, members),
      by = "combination"
    ) %>%
    dplyr::mutate(present = purrr::map2_lgl(members, set, ~ .y %in% .x)) %>%
    dplyr::select(-members)
  
  # Factor orders - REVERSED for highest on top
  comb_levels <- rev(combos$combination)  # Reverse order for highest first
  
  # Add number of intersections to combos
  combos <- combos %>%
    dplyr::mutate(
      n_intersections = stringr::str_count(combination, " & ") + 1,
      bar_color = ifelse(n_intersections >= multi_threshold, multi_bar_fill, bar_fill)
    )
  
  # Filter mat_long to only include combinations that exist in combos
  mat_long <- mat_long %>%
    dplyr::filter(combination %in% combos$combination) %>%
    dplyr::mutate(
      combination = factor(combination, levels = comb_levels),
      set         = factor(set, levels = order_sets)
    )
  
  combos <- combos %>%
    dplyr::mutate(combination = factor(combination, levels = comb_levels))
  
  # Optional lines connecting dots: compute numeric coords AFTER factoring
  lines_df <- NULL
  if (connect_lines) {
    lines_df <- mat_long %>%
      dplyr::filter(present) %>%
      dplyr::mutate(
        xn = as.integer(set),            # column index (set position)
        yn = as.integer(combination)     # row index (combination order)
      ) %>%
      dplyr::arrange(combination, xn)
  }
  
  # --- Panels ---
  
  # Top set-size bars (aligned with sets) - now with matching colors
  p_setsizes <-
    ggplot(set_sizes, aes(x = set, y = size, fill = set)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = set_colors, guide = "none") +
    scale_x_discrete(position = "top", drop = FALSE) +
    labs(x = NULL, y = "Set size") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(2, 5, 0, 5),  # Reduced top margin
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    )
  
  # Left: dot matrix with colored dots
  p_dots <-
    ggplot(mat_long, aes(x = set, y = combination)) +
    # faint grid dots for alignment
    geom_point(shape = 16, size = 1.8, color = "grey90") +
    # connecting lines behind filled dots
    { if (!is.null(lines_df))
      geom_path(
        data = lines_df,
        aes(x = xn, y = yn, group = combination),
        inherit.aes = FALSE, linewidth = line_width,
        color = line_color, lineend = "round"
      )
      else NULL } +
    # filled dots where present - now colored by set
    geom_point(data = dplyr::filter(mat_long, present), 
               aes(color = set), shape = 16, size = 2.8) +
    scale_color_manual(values = set_colors, guide = "none") +
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_y_discrete(limits = comb_levels, drop = FALSE,
                     expand = expansion(mult = c(0.02, 0.02))) +  # Small expansion for visual padding
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9, angle = 60, hjust = 0,
                                 color = set_colors[order_sets],
                                 face = "bold"),  # Color labels to match
      plot.margin = margin(0, 5, 5, 5),  # Adjusted margin
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    )
  
  # Right: horizontal bars
  # Double-check we have no zeros
  if(any(combos$size == 0)) {
    warning("Found combinations with size 0 that should have been filtered")
    combos <- combos %>% dplyr::filter(size > 0)
  }
  
  p_bars <- ggplot(combos, aes(y = combination, x = size, fill = bar_color)) +
    geom_col(width = 0.7) +
    scale_fill_identity() +  # Use the actual color values
    scale_y_discrete(limits = comb_levels, drop = FALSE,
                     expand = expansion(mult = c(0.02, 0.02))) +  # Same expansion as dot plot
    scale_x_continuous(position = "top") +  # Move x-axis to top
    labs(x = "Intersection size", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(2, 5, 5, 0),  # Reduced top margin
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    )
  
  # Apply scale transformations based on scale_type
  if (scale_type == "log10") {
    # Log scale (adding 1 to handle zeros)
    p_bars <- p_bars + 
      scale_x_continuous(
        trans = scales::pseudo_log_trans(base = 10),
        breaks = c(0, 1, 10, 100, 1000),
        labels = scales::comma,
        position = "top"  # Keep position on top
      )
  } else if (scale_type == "sqrt") {
    # Square root scale
    p_bars <- p_bars + 
      scale_x_continuous(
        trans = "sqrt",
        breaks = scales::pretty_breaks(n = 5),
        labels = scales::comma,
        position = "top"  # Keep position on top
      )
  } else if (scale_type == "break" && !is.null(break_point)) {
    # Custom broken axis
    # This is more complex - we'll use a transformation that compresses the upper range
    compress_trans <- function(break_point, break_ratio) {
      trans <- function(x) {
        ifelse(x <= break_point, 
               x * break_ratio / break_point,
               break_ratio + (x - break_point) * (1 - break_ratio) / (max(x) - break_point))
      }
      
      inv <- function(x) {
        max_val <- max(combos$size)
        ifelse(x <= break_ratio,
               x * break_point / break_ratio,
               break_point + (x - break_ratio) * (max_val - break_point) / (1 - break_ratio))
      }
      
      scales::trans_new("compress", trans, inv)
    }
    
    p_bars <- p_bars + 
      scale_x_continuous(
        trans = compress_trans(break_point, break_ratio),
        breaks = c(0, break_point/2, break_point, 
                   break_point + (max(combos$size) - break_point)/2, 
                   max(combos$size)),
        labels = scales::comma,
        position = "top"  # Keep position on top
      ) +
      # Add visual indicator of break
      annotate("segment", 
               x = break_point, xend = break_point,
               y = 0.5, yend = length(comb_levels) + 0.5,
               linetype = "dashed", color = "grey70", size = 0.5)
  } else {
    # Default linear scale with optimized limits
    p_bars <- p_bars + 
      scale_x_continuous(
        limits = c(0, max(combos$size) * 1.05),  # 5% padding
        expand = c(0, 0),
        labels = scales::comma,
        position = "top"  # Keep position on top
      )
  }
  
  # Create empty plot for alignment if showing set sizes
  if (show_set_sizes) {
    # Try patchwork instead of cowplot for better control
    if (requireNamespace("patchwork", quietly = TRUE)) {
      library(patchwork)
      
      p_empty <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "transparent", color = NA),
                                                 panel.background = element_rect(fill = "transparent", color = NA),
                                                 legend.background = element_rect(fill = "transparent", color = NA))
      
      final <- (p_setsizes | p_empty) / (p_dots | p_bars) +
        plot_layout(widths = c(1, 1.6), heights = c(1, 5))
    } else {
      # Fallback to cowplot with adjusted parameters
      p_empty <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "transparent", color = NA),
                                                 panel.background = element_rect(fill = "transparent", color = NA),
                                                 legend.background = element_rect(fill = "transparent", color = NA))
      
      # Try assembling without align parameter first
      top_row <- cowplot::plot_grid(p_setsizes, p_empty, ncol = 2, rel_widths = c(1, 1.6))
      bottom_row <- cowplot::plot_grid(p_dots, p_bars, ncol = 2, rel_widths = c(1, 1.6), align = "h")
      
      final <- cowplot::plot_grid(
        top_row, bottom_row,
        ncol = 1,
        rel_heights = c(1, 5)  # Much smaller top section
      )
    }
  } else {
    # Without set sizes, just align the two main plots
    final <- cowplot::plot_grid(
      p_dots, p_bars,
      ncol = 2, 
      rel_widths = c(1, 1.6), 
      align = "h", 
      axis = "tb"
    )
  }
  
  final <- final + 
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )
  
  return(final)
}

# Example usage with your data:
# vertical_upset(
#   upset_combined_dat,
#   sets = c("PT", "TAL", "EC", "IC", "IMMUNE", "VSMC_P_FIB", "POD"),
#   show_set_sizes = TRUE,
#   connect_lines  = TRUE, 
#   min_size = 3
# )

# Example with custom colors and multi-intersection highlighting:
# my_colors <- c(
#   "PT" = "#8B7D6B",
#   "TAL" = "#6B8E9F", 
#   "EC" = "#8FA68E",
#   "IC" = "#B8A9C9",
#   "IMMUNE" = "#D4A76A",
#   "VSMC_P_FIB" = "#9B6B6B",
#   "POD" = "#7A9A9F"
# )
# 
# vertical_upset(
#   upset_combined_dat,
#   sets = c("PT", "TAL", "EC", "IC", "IMMUNE", "VSMC_P_FIB", "POD"),
#   show_set_sizes = TRUE,
#   connect_lines  = TRUE, 
#   min_size = 3,
#   set_colors = my_colors,
#   multi_bar_fill = "#4A90A4",  # Color for 3+ intersections
#   multi_threshold = 3          # Highlight bars with 3 or more cell types
# )


# ===========================================================================
# Function: shorten_pathway_names
# ===========================================================================

shorten_pathway_names <- function(pathway_names, max_length = 40, aggressive = FALSE) {
  
  # Define replacement patterns
  # [Keep all your existing replacements as they are]
  word_replacements <- c(
    "Hallmark " = "",
    "Role of " = "",
    "Function of " = "",
    "Eukaryotic " = "",
    "in the" = "in",
    "by the" = "by",
    "of the" = "of",
    "Phosphorylation" = "Phosph",
    "Metabolism" = "Metab",
    "Extracellular" = "EC",
    "Intracellular" = "IC",
    "Dysfunction" = "Adaptation",
    "Alternative" = "Alt",
    "Dependent" = "Dep",
    "Independent" = "Indep",
    "Associated" = "Assoc",
    "Inflammatory" = "Inflam",
    "Inflammation" = "Inflam",
    "Lymphocyte" = "Lymph",
    "Lymphocytes" = "Lymph",
    "Fibroblast" = "Fib",
    "Fibroblasts" = "Fib",
    "Endothelial" = "Endo",
    "Epithelial" = "Epi",
    "Transition" = "Trans",
    "Rheumatoid Arthritis" = "RA",
    "Multiple Sclerosis" = "MS",
    "Alzheimer's Disease" = "Alzheimer's",
    "Parkinson's Disease" = "Parkinson's",
    "Respiratory" = "Resp",
    "Gastrointestinal" = "GI",
    "Central Nervous System" = "CNS",
    "Oxidative" = "Ox",
    "between" = "b/w",
    "through" = "via",
    "and" = "&"
  )
  
  # More aggressive replacements
  aggressive_replacements <- c(
    " the " = " ",
    " in " = " ",
    " of " = " ",
    " to " = "→",
    " from " = "←"
  )
  
  # Cell/molecule specific replacements
  molecule_replacements <- c(
    "Natural Killer" = "NK",
    "Tumor Necrosis Factor" = "TNF",
    "Transforming Growth Factor" = "TGF",
    "Platelet Derived Growth Factor" = "PDGF",
    "Epidermal Growth Factor" = "EGF",
    "Vascular Endothelial Growth Factor" = "VEGF",
    "Insulin-like Growth Factor" = "IGF",
    "Interferon" = "IFN",
    "Interleukin" = "IL",
    "G-Protein Coupled" = "GPCR",
    "G Protein Coupled" = "GPCR",
    "Activator Protein" = "AP",
    "Electron Transport" = "e- Transport",
    "Nitric Oxide" = "NO",
    "Reactive Oxygen Species" = "ROS",
    "Amino Acid" = "AA",
    "Fatty Acid" = "FA"
  )
  
  # Function to apply replacements
  apply_replacements <- function(text, replacements) {
    for (pattern in names(replacements)) {
      text <- gsub(pattern, replacements[pattern], text, ignore.case = FALSE)
    }
    return(text)
  }
  
  # Track which pathways got truncated at Step 6
  step6_truncated <- rep(FALSE, length(pathway_names))
  
  # Process each pathway name
  shortened <- sapply(seq_along(pathway_names), function(i) {
    pathway <- pathway_names[i]
    if (is.na(pathway)) return(NA)
    
    # Step 1: Apply molecule/cell replacements first
    result <- apply_replacements(pathway, molecule_replacements)
    
    # Step 2: Apply standard word replacements
    result <- apply_replacements(result, word_replacements)
    
    # Step 3: If still too long and aggressive mode is on
    if (nchar(result) > max_length && aggressive) {
      result <- apply_replacements(result, aggressive_replacements)
    }
    
    # Step 4: Remove parenthetical content if still too long
    if (nchar(result) > max_length) {
      result <- gsub(" \\([^)]+\\)", "", result)
    }
    
    # # Step 5: If still too long, abbreviate remaining long words
    # if (nchar(result) > max_length) {
    #   words <- strsplit(result, " ")[[1]]
    #   # Abbreviate words longer than 8 characters that aren't already abbreviated
    #   words <- sapply(words, function(word) {
    #     if (nchar(word) > 8 && !grepl("[A-Z]{2,}", word)) {
    #       # Keep first 3-4 letters and last letter
    #       if (nchar(word) > 10) {
    #         paste0(substr(word, 1, 3), ".", substr(word, nchar(word), nchar(word)))
    #       } else {
    #         paste0(substr(word, 1, 4), ".")
    #       }
    #     } else {
    #       word
    #     }
    #   })
    #   result <- paste(words, collapse = " ")
    # }
    
    # Step 6: Final truncation if needed
    if (nchar(result) > max_length) {
      step6_truncated[i] <<- TRUE  # Mark this pathway as truncated at Step 6
      result <- paste0(substr(result, 1, max_length - 3), "...")
    }
    
    return(result)
  })
  
  # Return a list with both the shortened names and the Step 6 truncation info
  return(list(
    shortened = shortened,
    step6_truncated = step6_truncated
  ))
}

# Helper function to create a lookup table for manual review
create_pathway_lookup <- function(original, shortened) {
  data.frame(
    original = original,
    shortened = shortened,
    length_original = nchar(original),
    length_shortened = nchar(shortened),
    reduction = round((1 - nchar(shortened)/nchar(original)) * 100, 1),
    stringsAsFactors = FALSE
  )
}

# Extract the legend from one of the plots (e.g., p2)
get_legend <- function(myplot) {
  tmp <- ggplot_gtable(ggplot_build(myplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# ===========================================================================
# Function: concordance_by_celltype
# ===========================================================================

concordance_by_celltype <- function(
    biofluid_df,
    scrna_df,
    biofluid_label   = "urine",   # could be "plasma"
    urine_gene_col   = Gene,
    urine_logFC_col  = `logFC_log2(protein)`,
    urine_p_col      = `p_log2(protein)`,
    scrna_gene_col   = Gene,
    scrna_logFC_col  = `logFC_treatmentDapagliflozin:visitPOST`,
    scrna_p_col      = `p_treatmentDapagliflozin:visitPOST`,
    p_thresh         = 0.05,
    celltype         = "PT",
    out_dir          = ".",
    include_not_in_both = TRUE
) {
  stopifnot(requireNamespace("dplyr", quietly = TRUE),
            requireNamespace("tidyr", quietly = TRUE),
            requireNamespace("readr", quietly = TRUE),
            requireNamespace("rlang", quietly = TRUE))
  
  library(dplyr); library(tidyr); library(rlang)
  
  # Tidy-eval symbols
  bg   <- enquo(urine_gene_col)
  blfc <- enquo(urine_logFC_col)
  bp   <- enquo(urine_p_col)
  sg   <- enquo(scrna_gene_col)
  slfc <- enquo(scrna_logFC_col)
  sp   <- enquo(scrna_p_col)
  
  standardize <- function(df, gene_q, logfc_q, p_q, source_label) {
    df %>%
      filter(!is.na(!!p_q), !is.na(!!logfc_q), !!p_q < p_thresh) %>%
      mutate(source = source_label) %>%
      transmute(
        source,
        Gene  = as.character(!!gene_q),
        logFC = as.numeric(!!logfc_q),
        p     = as.numeric(!!p_q)
      )
  }
  
  biofluid_std <- standardize(biofluid_df, bg, blfc, bp, biofluid_label)
  scrna_std    <- standardize(scrna_df, sg, slfc, sp, "scRNA")
  
  concordance_tbl <- bind_rows(biofluid_std, scrna_std) %>%
    mutate(.dir = case_when(
      logFC > 0 ~ "pos",
      logFC < 0 ~ "neg",
      TRUE ~ NA_character_
    )) %>%
    group_by(Gene, source) %>%
    summarise(
      n_entries = n(),
      n_pos = sum(.dir == "pos", na.rm = TRUE),
      n_neg = sum(.dir == "neg", na.rm = TRUE),
      source_dir = case_when(
        n_pos > 0 & n_neg > 0 ~ "mixed",
        n_pos > 0 ~ "pos",
        n_neg > 0 ~ "neg",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    ) %>%
    dplyr::select(Gene, source, source_dir, n_entries) %>%
    pivot_wider(
      names_from = source,
      values_from = c(source_dir, n_entries),
      names_sep = "__"
    ) %>%
    mutate(
      concordance = case_when(
        is.na(!!sym(paste0("source_dir__", biofluid_label))) |
          is.na(source_dir__scRNA) ~ "not_in_both",
        !!sym(paste0("source_dir__", biofluid_label)) == "mixed" |
          source_dir__scRNA == "mixed" ~ "mixed",
        !!sym(paste0("source_dir__", biofluid_label)) == "pos" &
          source_dir__scRNA == "pos" ~ "positive",
        !!sym(paste0("source_dir__", biofluid_label)) == "neg" &
          source_dir__scRNA == "neg" ~ "negative",
        TRUE ~ "non-concordant"
      )
    ) %>%
    { if (!include_not_in_both) filter(., concordance != "not_in_both") else . } %>%
    arrange(match(concordance, c("mixed","positive","negative","non-concordant","not_in_both")),
            Gene)
  
  out_file <- file.path(out_dir,
                        paste0(tolower(celltype), "_", biofluid_label, "_concordance.csv"))
  readr::write_csv(concordance_tbl, out_file)
  
  print(table(concordance_tbl$concordance, useNA = "ifany"))
  message("Saved: ", out_file)
  
  return(concordance_tbl)
}


# ===========================================================================
# Function: plot_gsea_results
# ===========================================================================

plot_gsea_results <- function(gsea_list, 
                              cell_name,
                              reference = "hallmark",
                              top_n = 20,
                              max_pathway_length = 45,
                              caption_width = 70,
                              p_threshold = 0.05,
                              min_x_limit = 5,
                              low_color = "#89c2d9",
                              mid_color = "white", 
                              high_color = "#ee7674",
                              show_truncated_in_caption = TRUE) {
  
  # Extract the specified reference results
  if (!reference %in% names(gsea_list)) {
    stop(paste0("Reference '", reference, "' not found in gsea_list. ",
                "Available references: ", paste(names(gsea_list), collapse = ", ")))
  }
  
  if (reference == "go") {
    gsea_list[[reference]] <- gsea_list[[reference]] %>%
      filter(grepl("^GOBP_", pathway))
  }
  # Extract and prepare top pathways from the specified reference
  top_pathways <- gsea_list[[reference]] %>%
    arrange(pval) %>%
    head(top_n)
  
  # Shorten pathway names
  shorten_result <- shorten_pathway_names(top_pathways$pathway, max_length = max_pathway_length)
  
  top_pathways <- top_pathways %>%
    mutate(
      was_step6_truncated = shorten_result$step6_truncated,
      shortened_pathway = shorten_result$shortened,
      clean_pathway = clean_pathway_names(shortened_pathway),
      neg_log_p = -log10(pval)
    )
  
  # Create caption
  base_caption <- paste0("Cell type: ", gsub("_", " ", cell_name), 
                         " | Reference: ", toupper(reference))
  
  if (show_truncated_in_caption) {
    # Get truncated pathways for caption
    truncated_pathways <- top_pathways %>%
      filter(was_step6_truncated) %>%
      dplyr::mutate(
        full_clean = clean_pathway_names(pathway),
        full_clean_wrapped = str_wrap(full_clean, width = caption_width, indent = 0, exdent = 4)
      ) %>%
      dplyr::select(clean_pathway, full_clean_wrapped)
    
    if (nrow(truncated_pathways) > 0) {
      truncated_text <- paste(
        apply(truncated_pathways, 1, function(x) {
          paste0(x["full_clean_wrapped"])
        }), 
        collapse = "\n"
      )
      full_caption <- paste0(base_caption, "\n\nTruncated pathways:\n", truncated_text)
    } else {
      full_caption <- base_caption
    }
  } else {
    full_caption <- base_caption
  }
  
  # Set factor levels for plotting order
  top_pathways$clean_pathway <- factor(top_pathways$clean_pathway, 
                                       levels = rev(top_pathways$clean_pathway))
  
  # Determine x-axis limits
  if (!is.null(min_x_limit)) {
    actual_max <- max(top_pathways$neg_log_p, na.rm = TRUE)
    x_upper_limit <- max(actual_max, min_x_limit)
  } else {
    x_upper_limit <- NULL  # Let ggplot2 determine automatically
  }
  
  # Create plot
  p <- top_pathways %>%
    ggplot(aes(y = clean_pathway, x = neg_log_p, fill = NES)) +
    geom_col(width = 0.9) + 
    geom_vline(xintercept = -log10(p_threshold), linetype = "dashed", color = "#aaaaaa") +
    geom_text(aes(label = clean_pathway), 
              x = -log10(p_threshold) + 0.1, hjust = 0, 
              fontface = "bold", family = "Arial",
              color = "#2b2b2b") +
    scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, 
                         midpoint = 0,
                         guide = guide_colorbar(barheight = 0.4, barwidth = 8)) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top",
          axis.text.y = element_blank(),
          legend.title = element_text(vjust = 0.8),
          plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.2)) +
    labs(y = NULL, 
         x = "-log(p-value)", 
         fill = "NES",
         caption = full_caption,
         title = cell_name)
  
  # Apply x-axis limits
  if (!is.null(min_x_limit)) {
    p <- p + scale_x_continuous(limits = c(0, x_upper_limit), expand = c(0, 0))
  } else {
    p <- p + scale_x_continuous(expand = expansion(mult = c(0, 0)))
  }
  
  return(p)
}

# ===========================================================================
# Function: calculate_celltype_proportions
# ===========================================================================

calculate_celltype_proportions <- function(data, bin_width = 10, group_by = "KPMP_celltype") {
  
  # Create pseudotime bins
  data_with_bins <- data %>%
    dplyr::mutate(
      bin_start = floor(slingPseudotime_1 / bin_width) * bin_width,
      bin_end = bin_start + bin_width,
      bin_label = paste0(bin_start, "-", bin_end)
    )
  
  # Calculate proportions for each bin
  proportions <- data_with_bins %>%
    group_by(bin_label, bin_start, !!sym(group_by)) %>%
    dplyr::summarise(n = n(), .groups = "drop_last") %>%
    dplyr::mutate(
      total_in_bin = sum(n),
      proportion = (n / total_in_bin) * 100
    ) %>%
    ungroup()
  
  # Ensure all groups are represented in each bin
  all_bins <- unique(proportions$bin_label)
  all_groups <- unique(data[[group_by]])
  
  complete_proportions <- proportions %>%
    complete(
      bin_label = all_bins,
      !!sym(group_by) := all_groups,
      fill = list(n = 0, proportion = 0)
    ) %>%
    group_by(bin_label) %>%
    dplyr::mutate(
      bin_start = min(bin_start, na.rm = TRUE),
      total_in_bin = sum(n)
    ) %>%
    ungroup() %>%
    arrange(bin_start, !!sym(group_by))
  
  return(complete_proportions)
}

# ===========================================================================
# Function: create_pie_chart
# ===========================================================================

# Function to create a single pie chart for a given bin
# Requires: library(dplyr); library(ggplot2); library(ggtext); library(tidyr); library(rlang)
create_pie_chart <- function(data,
                             bin_value,
                             color_palette,
                             group_by = "KPMP_celltype",
                             caption_groups = NULL,
                             digits = 1) {
  
  grp <- rlang::sym(group_by)
  
  bin_data <- data %>%
    dplyr::filter(bin_start == bin_value)
  
  if (nrow(bin_data) == 0) {
    stop(sprintf("No rows found for bin_start == %s", bin_value))
  }
  if (!group_by %in% names(bin_data)) {
    stop(sprintf("Column '%s' not found in data.", group_by))
  }
  
  # Build caption table (one row per group)
  cap_df <- bin_data %>%
    dplyr::select(!!grp, proportion) %>%
    dplyr::group_by(!!grp) %>%
    dplyr::summarise(proportion = sum(proportion, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(proportion))
  
  # Optional: enforce caption order/subset
  if (!is.null(caption_groups)) {
    cap_df <- cap_df %>%
      dplyr::filter(!!grp %in% caption_groups) %>%
      dplyr::mutate(!!grp := factor(!!grp, levels = caption_groups)) %>%
      dplyr::arrange(!!grp)
    
    # also apply to plotting data so legend/order matches
    bin_data[[group_by]] <- factor(bin_data[[group_by]], levels = caption_groups)
  }
  
  # Caption HTML
  caption_text <- paste(
    vapply(seq_len(nrow(cap_df)), function(i) {
      g   <- as.character(cap_df[[group_by]][i])
      pct <- round(cap_df$proportion[i], digits)
      col <- if (!is.null(color_palette[[g]])) color_palette[[g]] else "#000000"
      sprintf("<span style='color:%s'>%s: %s%%</span>", col, g, pct)
    }, character(1)),
    collapse = "<br>"
  )
  
  # Limit palette to groups present
  present_groups <- unique(as.character(bin_data[[group_by]]))
  pal_used <- color_palette[names(color_palette) %in% present_groups]
  
  ggplot(bin_data, aes(x = "", y = proportion, fill = !!grp)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    theme_void() +
    labs(fill = NULL, caption = caption_text) +
    theme(plot.caption = ggtext::element_markdown(hjust = 0.5, face = "bold")) +
    scale_fill_manual(values = pal_used)
}