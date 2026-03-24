

# INSTABILITY UNDER SMALL PERTURBATIONS ########################################
################################################################################
################################################################################

# RESOLUTION SPECIFIC TAU ######################################################

get_tau_mha <- function(tau_percent) {
  cell_area_nominal <- c(
    "0.2deg" = 0.05,  
    "0.4deg" = 0.10,  
    "1deg"   = 0.30   
  )
  
  tau_frac <- tau_percent / 100
  tau_frac * cell_area_nominal
}

# CLASSIFICATION FOR ONE ENSEMBLE + ONE TAU ####################################

classify_cells <- function(dt, dataset_subset, tau_percent) {
  
  rr <- unique(dt$resolution)
  if (length(rr) != 1) stop("dt must contain only one resolution")
  
  tau_mha <- unname(get_tau_mha(tau_percent)[rr])
  
  A_mat <- as.matrix(dt[, ..dataset_subset])
  storage.mode(A_mat) <- "numeric"
  
  pres_mat <- A_mat > tau_mha
  n_present <- rowSums(pres_mat, na.rm = TRUE)
  
  out <- dt[, .(lon, lat, resolution)]
  out[, `:=`(tau_percent = tau_percent,
             ensemble_size = length(dataset_subset),
             n_present = n_present,
             classified_present = n_present > 0)]
  
  out
}

# PERTURBED CLASSIFICATIONS ####################################################

# On a full ensemble, leave-one-out ensemble and different taus ----------------

build_perturbation_runs <- function(dt,
                                    dataset_names,
                                    tau_vals = c(1, 3, 5),
                                    include_full = TRUE,
                                    include_loo = TRUE) {
  
  rr <- unique(dt$resolution)
  if (length(rr) != 1) stop("dt must contain only one resolution")
  
  runs <- list()
  idx <- 1L
  
  # full ensemble ----------------------------------------------
  
  if (include_full) {
    for (tt in tau_vals) {
      tmp <- classify_cells(
        dt = dt,
        dataset_subset = dataset_names,
        tau_percent = tt
      )
      tmp[, perturbation:= paste0("full_tau", tt)]
      runs[[idx]] <- tmp
      idx <- idx + 1L
    }
  }
  
  # leave-one-out ensembles -------------------------------------
  
  if (include_loo) {
    for (drop_name in dataset_names) {
      subset_j <- setdiff(dataset_names, drop_name)
      for (tt in tau_vals) {
        tmp <- classify_cells(
          dt = dt,
          dataset_subset = subset_j,
          tau_percent = tt
        )
        tmp[, perturbation:= paste0("drop_", drop_name, "_tau", tt)]
        runs[[idx]] <- tmp
        idx <- idx + 1L
      }
    }
  }
  
  rbindlist(runs, fill = TRUE)
}

# COMPUTE CELL LEVEL INSTABILITY ###############################################

compute_cell_instability <- function(class_dt) {
  
  class_dt[, .(n_runs = .N,
               n_present_runs = sum(classified_present, na.rm = TRUE),
               share_present = mean(classified_present, na.rm = TRUE),
               
               # number of transitions after sorting perturbations 
               # alphabetically ------------------------------------------------
               
               n_changes = {
                 x <- classified_present
                 sum(x[-1] != x[-length(x)])
               },
               
               # fraction of runs not in the modal class -----------------------
               
               flip_rate = 1 - max(mean(classified_present), mean(!classified_present)),
               
               # fully stable classes-------------------------------------------
               
               always_absent = all(!classified_present),
               always_present = all(classified_present),
               unstable = any(classified_present) & !all(classified_present)
               
  ), by = .(lon, lat, resolution)]
}


# RUN INSTABILITY ANALYSIS #####################################################

run_instability_analysis <- function(dt,
                                     dataset_names,
                                     tau_vals = c(1, 3, 5),
                                     include_full = TRUE,
                                     include_loo = TRUE) {
  
  res_levels <- unique(dt$resolution)
  
  out_class <- rbindlist(lapply(res_levels, function(rr) {
    dt_rr <- dt[resolution == rr]
    build_perturbation_runs(
      dt = dt_rr,
      dataset_names = dataset_names,
      tau_vals = tau_vals,
      include_full = include_full,
      include_loo = include_loo
    )
  }), fill = TRUE)
  
  out_instability <- compute_cell_instability(out_class)
  
  list(classifications = out_class, instability = out_instability)
}

# SUMMARY TABLES ###############################################################

summarize_instability <- function(instability_dt) {
  
  instability_dt[, .(
    n_cells = .N,
    
    frac_always_absent = mean(always_absent),
    frac_always_present = mean(always_present),
    frac_unstable = mean(unstable),
    
    mean_flip_rate = mean(flip_rate, na.rm = TRUE),
    median_flip_rate = median(flip_rate, na.rm = TRUE),
    
    q05_flip_rate = quantile(flip_rate, 0.05, na.rm = TRUE),
    q25_flip_rate = quantile(flip_rate, 0.25, na.rm = TRUE),
    q75_flip_rate = quantile(flip_rate, 0.75, na.rm = TRUE),
    q95_flip_rate = quantile(flip_rate, 0.95, na.rm = TRUE)
    
  ), by = resolution][order(resolution)]
}

# Conditional summary (only on cells ever detected as irrigated)

summarize_instability_relevant <- function(instability_dt) {
  
  instability_dt[n_present_runs > 0, .(
    
    n_cells = .N,
    frac_unstable = mean(unstable),
    mean_flip_rate = mean(flip_rate, na.rm = TRUE),
    median_flip_rate = median(flip_rate, na.rm = TRUE),
    q05_flip_rate = quantile(flip_rate, 0.05, na.rm = TRUE),
    q25_flip_rate = quantile(flip_rate, 0.25, na.rm = TRUE),
    q75_flip_rate = quantile(flip_rate, 0.75, na.rm = TRUE),
    q95_flip_rate = quantile(flip_rate, 0.95, na.rm = TRUE)
  ), by = resolution][order(resolution)]
}

# CONVERGENCE TEST #############################################################
################################################################################
################################################################################

# METRICS FOR SUBSET OF DATASETS ###############################################

compute_subset_agreement <- function(dt, dataset_subset, tau_percent) {
  
  rr <- unique(dt$resolution)
  if (length(rr) != 1) stop("dt must contain only one resolution")
  
  tau_mha_map <- get_tau_mha(tau_percent)
  tau_mha <- unname(tau_mha_map[rr])
  
  A_mat <- as.matrix(dt[, ..dataset_subset])
  storage.mode(A_mat) <- "numeric"
  
  pres_mat <- A_mat > tau_mha
  
  n_present <- rowSums(pres_mat, na.rm = TRUE)
  k <- length(dataset_subset)
  
  all_absent <- n_present == 0
  all_present <- n_present == k
  full_agree <- all_absent | all_present
  disagreement <- n_present > 0 & n_present < k
  relevant <- n_present > 0
  
  p_present <- n_present / k
  
  # Strict ambiguity zone -----------------------------------------
  
  amb_40_60 <- relevant & p_present >= 0.4 & p_present <= 0.6
  
  # Broad instability zone ---------------------------------------
  
  inst_20_80 <- relevant & p_present >= 0.2 & p_present <= 0.8
  
  data.table(resolution = rr,
             tau_percent = tau_percent,
             tau_mha = tau_mha,
             k = k,
             n_cells = nrow(dt),
             n_relevant = sum(relevant),
    
             # Unconditional -------------------
             
            frac_any_irrig = mean(relevant),
            frac_full_agree = mean(full_agree),
            frac_disagree = mean(disagreement),
    
             # Conditional on at least one positive detection
            
            frac_disagree_cond =
              if (sum(relevant) > 0) mean(disagreement[relevant]) else NA_real_,
            
            frac_amb_40_60_cond =
              if (sum(relevant) > 0) mean(amb_40_60[relevant]) else NA_real_,
            
            frac_inst_20_80_cond =
              if (sum(relevant) > 0) mean(inst_20_80[relevant]) else NA_real_)
}


# RANDOM-SUBSET TEST ###########################################################

run_agreement_vs_information <- function(dt,
                                         dataset_names,
                                         tau_percent = 1,
                                         k_values = 2:length(dataset_names),
                                         n_reps = 250,
                                         seed = 123) {
  
  set.seed(seed)
  
  out_list <- vector("list", length(k_values) * n_reps)
  idx <- 1L
  
  for (k in k_values) {
    for (r in seq_len(n_reps)) {
      
      dataset_subset <- sample(dataset_names, size = k, replace = FALSE)
      
      tmp <- compute_subset_agreement(dt = dt,
                                      dataset_subset = dataset_subset,
                                      tau_percent = tau_percent)
      
      tmp[, `:=`(rep = r,
                 subset_type = "random",
                 dataset_subset = paste(sort(dataset_subset), collapse = ";"))]
      
      out_list[[idx]] <- tmp
      idx <- idx + 1L
    }
  }
  
  rbindlist(out_list, fill = TRUE)
}


# EXACT TEST WITH ALL COMBINATIONS #############################################

run_agreement_vs_information_exact <- function(dt,
                                               dataset_names,
                                               tau_percent = 1,
                                               k_values = 2:length(dataset_names)) {
  
  out_list <- list()
  idx <- 1L
  
  for (k in k_values) {
    combs <- combn(dataset_names, k, simplify = FALSE)
    
    for (j in seq_along(combs)) {
      subset_j <- combs[[j]]
      
      tmp <- compute_subset_agreement(dt = dt,
                                      dataset_subset = subset_j,
                                      tau_percent = tau_percent)
      
      tmp[, `:=`(comb_id = j,
                 subset_type = "exact",
                 dataset_subset = paste(sort(subset_j), collapse = ";"))]
      
      out_list[[idx]] <- tmp
      idx <- idx + 1L
    }
  }
  
  rbindlist(out_list, fill = TRUE)
}


# RUN ACROSS RESOLUTIONS #######################################################

run_by_resolution <- function(dt,
                              dataset_names,
                              tau_percent = 1,
                              mode = c("random", "exact"),
                              n_reps = 250,
                              seed = 123,
                              k_values = 2:length(dataset_names)) {
  
  mode <- match.arg(mode)
  res_levels <- unique(dt$resolution)
  
  out <- rbindlist(lapply(res_levels, function(rr) {
    
    dt_rr <- dt[resolution == rr]
    
    if (mode == "random") {
      run_agreement_vs_information(
        dt = dt_rr,
        dataset_names = dataset_names,
        tau_percent = tau_percent,
        k_values = k_values,
        n_reps = n_reps,
        seed = seed
      )
    } else {
      run_agreement_vs_information_exact(
        dt = dt_rr,
        dataset_names = dataset_names,
        tau_percent = tau_percent,
        k_values = k_values
      )
    }
  }), fill = TRUE)
  
  out[]
}


# SUMMARIZE RESULTS ###########################################################

summarize_agreement_vs_information <- function(res_dt) {
  
  res_dt[, .(
    
    n_subsets = .N,
    
    frac_any_irrig_mean = mean(frac_any_irrig, na.rm = TRUE),
    frac_any_irrig_q05 = quantile(frac_any_irrig, 0.05, na.rm = TRUE),
    frac_any_irrig_q95 = quantile(frac_any_irrig, 0.95, na.rm = TRUE),
    
    frac_full_agree_mean = mean(frac_full_agree, na.rm = TRUE),
    frac_full_agree_q05 = quantile(frac_full_agree, 0.05, na.rm = TRUE),
    frac_full_agree_q95 = quantile(frac_full_agree, 0.95, na.rm = TRUE),
    
    frac_disagree_mean = mean(frac_disagree, na.rm = TRUE),
    frac_disagree_q05 = quantile(frac_disagree, 0.05, na.rm = TRUE),
    frac_disagree_q95 = quantile(frac_disagree, 0.95, na.rm = TRUE),
    
    frac_disagree_cond_mean = mean(frac_disagree_cond, na.rm = TRUE),
    frac_disagree_cond_q05 = quantile(frac_disagree_cond, 0.05, na.rm = TRUE),
    frac_disagree_cond_q95 = quantile(frac_disagree_cond, 0.95, na.rm = TRUE),
    
    frac_amb_40_60_cond_mean = mean(frac_amb_40_60_cond, na.rm = TRUE),
    frac_amb_40_60_cond_q05 = quantile(frac_amb_40_60_cond, 0.05, na.rm = TRUE),
    frac_amb_40_60_cond_q95 = quantile(frac_amb_40_60_cond, 0.95, na.rm = TRUE),
    
    frac_inst_20_80_cond_mean = mean(frac_inst_20_80_cond, na.rm = TRUE),
    frac_inst_20_80_cond_q05 = quantile(frac_inst_20_80_cond, 0.05, na.rm = TRUE),
    frac_inst_20_80_cond_q95 = quantile(frac_inst_20_80_cond, 0.95, na.rm = TRUE)
    
  ), by = .(resolution, k, tau_percent, tau_mha)][order(resolution, tau_percent, k)]
}


################################################################################
################################################################################


