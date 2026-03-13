
# COMPUTE TAU FUNCTION #########################################################

compute_tau_fun <- function(dt, tau_mha, dataset_names) {
  tmp <- copy(dt)

  # Raw area matrix (rows = cells, cols = datasets) ----------------------------
  A_mat <- as.matrix(tmp[, ..dataset_names])

  # Presence / absence at this tau ---------------------------------------------

  pres_mat <- A_mat > tau_mha

  # Counts per cell ------------------------------------------------------------

  n_present <- rowSums(pres_mat, na.rm = TRUE)
  n_datasets <- length(dataset_names)
  n_absent <- n_datasets - n_present

  # Derived indicators ---------------------------------------------------------

  any_irrig <- n_present > 0
  no_irrig <- n_present == 0
  all_irrig <- n_present == n_datasets
  disagreement <- n_present > 0 & n_present < n_datasets
  share_present <- n_present / n_datasets

  # Presence classes -----------------------------------------------------------

  presence_class <- fifelse(n_present == 0,
                            "None",
                            fifelse(
                              n_present == n_datasets,
                              "All",
                              fifelse(share_present < 0.5, "Minority presence",
                                      "Majority presence")
                            ))

  # Majority structure ---------------------------------------------------------

  bare_majority <- n_present == ceiling(n_datasets / 2)
  strong_majority <- n_present >= (n_datasets - 1) & n_present > 0

  # Disagreement structure: minor vs major -------------------------------------

  # Minor disagreement: exactly one dataset disagrees with the rest

  minor_disagreement <- disagreement &
    (n_present == 1L | n_absent == 1L)

  # Major disagreement: at least 2 datasets on each side

  major_disagreement <- disagreement &
    (n_present >= 2L & n_absent >= 2L)

  # intensity-based decomposition: marginal vs existential ---------------------
  # minimum irrigated area across datasets for each cell

  A_min <- do.call(pmin, c(as.data.frame(A_mat), list(na.rm = TRUE)))

  # marginal: all datasets see some irrigation, but not all exceed tau

  marginal   <- disagreement & (A_min > 0)

  # existential: at least one dataset reports zero irrigation, others exceed tau

  existence  <- disagreement & (A_min == 0)

  total_disagree <- sum(disagreement)

  frac_marginal   <- mean(marginal)
  frac_existence  <- mean(existence)

  share_marginal  <- if (total_disagree > 0)
    mean(marginal[disagreement])
  else
    NA_real_
  share_existence <- if (total_disagree > 0)
    mean(existence[disagreement])
  else
    NA_real_

  # Summary --------------------------------------------------------------------

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
