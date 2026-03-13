
# COMPUTE TAU INCLUDING WEIGHTS FOR DATASETS ###################################


compute_tau_weighted_fun <- function(dt,
                                     tau_mha,
                                     dataset_names,
                                     weight_label,
                                     weight_scenarios) {

  tmp <- copy(dt)

  # Raw area matrix (rows = cells, cols = datasets) ---------------------------

  A_mat <- as.matrix(tmp[, ..dataset_names])

  # Select the weight vector based on the label -------------------------------

  w_vec <- weight_scenarios[[weight_label]]
  if (is.null(w_vec)) {
    stop("Unknown weight label in compute_tau_weighted_fun: ", weight_label)
  }

  # Weights restricted to these datasets --------------------------------------

  w <- w_vec[dataset_names]

  if (any(is.na(w)) || length(w) != length(dataset_names)) {
    stop(
      "Weights missing or length mismatch for datasets under weight = ",
      weight_label, ". Missing: ",
      paste(dataset_names[is.na(w)], collapse = ", ")
    )
  }

  # Normalise (so they sum to 1) ----------------------------------------------

  w <- w / sum(w)

  # Presence / absence at this tau --------------------------------------------

  pres_mat <- A_mat > tau_mha

  # Weighted "counts" per cell (0–1) ------------------------------------------
  # p_present = sum_d w_d * 1(A_d > tau)

  p_present <- as.numeric(pres_mat %*% w)
  p_absent  <- 1 - p_present

  # Slight tolerance for exact 0 / 1 comparisons
  tol <- 1e-10

  # Derived indicators ---------------------------------------------------------

  any_irrig <- p_present > 0
  no_irrig <- p_present <= tol
  all_irrig <- p_present >= 1 - tol
  disagreement <- p_present > 0 & p_present < 1 - tol
  share_present <- p_present

  # Presence classes -----------------------------------------------------------

  presence_class <- fifelse(share_present <= tol, "None",
                            fifelse( share_present >= 1 - tol, "All",
                                     fifelse(share_present < 0.5, "Minority presence",
                                             "Majority presence")))

  # Majority structure (weighted analogues) -----------------------------------

  n_datasets <- length(dataset_names)
  max_w <- max(w)
  min_w <- min(w)

  # Bare majority:
  # original: n_present == ceiling(n/2)
  # weighted approximation: just above 0.5 but not too large

  bare_majority <- share_present > 0.5 & share_present <= 0.5 + 0.5 / n_datasets

  # Strong majority:
  # original: n_present >= n - 1
  # weighted analogue: all but (at most) the smallest-weight dataset
  strong_majority <- share_present >= (1 - min_w) & share_present > 0

  # Disagreement structure: minor vs major (weighted) -------------------------

  minor_disagreement <- disagreement & (share_present <= max_w | p_absent <= max_w)

  # Major disagreement: both camps have more than one dataset's worth of weight

  major_disagreement <- disagreement & (share_present > max_w & p_absent > max_w)

  # Intensity-based decomposition: marginal vs existential --------------------
  # Prefer using precomputed A_min if present in dt (from your main loop)

  if ("A_min" %in% names(tmp)) {
    A_min <- tmp[["A_min"]]
  } else {
    A_min <- do.call(pmin, c(as.data.frame(A_mat), list(na.rm = TRUE)))
  }

  # marginal: all datasets see some irrigation, but not all exceed tau

  marginal <- disagreement & (A_min > 0)

  # existential: at least one dataset reports zero irrigation, others exceed tau

  existence <- disagreement & (A_min == 0)

  total_disagree <- sum(disagreement)
  frac_marginal <- mean(marginal)
  frac_existence <- mean(existence)

  share_marginal <- if (total_disagree > 0)
    mean(marginal[disagreement])
  else
    NA_real_

  share_existence <- if (total_disagree > 0)
    mean(existence[disagreement])
  else
    NA_real_

  # Summary (same structure as compute_tau_fun) -------------------------------

  data.table(frac_any_irrig = mean(any_irrig),
             frac_no_irrig = mean(no_irrig),
             frac_disagree = mean(disagreement),
             frac_none = mean(presence_class == "None"),
             frac_all = mean(presence_class == "All"),
             frac_minority = mean(presence_class == "Minority presence"),
             frac_majority = mean(presence_class == "Majority presence"),
             frac_bare_majority = mean(bare_majority),
             frac_strong_majority = mean(strong_majority),
             frac_minor_disagree = mean(minor_disagreement),
             frac_major_disagree = mean(major_disagreement),
             frac_marginal = frac_marginal,
             frac_existence = frac_existence,
             share_marginal = share_marginal,
             share_existence = share_existence)
}
