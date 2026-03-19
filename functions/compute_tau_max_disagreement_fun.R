
# FUNCTION TO COMPUTE TAU MAX DISAGREEMENT #####################################

# tau max = largest tau such that:
# 1) at least one dataset reports <= threshold ha (treated as zero), and
# 2) at least one dataset reports >= tau
# Cells where all datasets are > threshold ha return NA.

compute_tau_max_existential_fun <- function(dt,
                                            dataset_names,
                                            taus_mha,
                                            taus_ha,
                                            zero_tol_ha) {

  # Convert tolerance to Mha
  zero_tol_mha <- zero_tol_ha / 1e6

  # Matrix: rows = cells, cols = datasets (values in Mha)
  A <- as.matrix(dt[, ..dataset_names])

  n_cells <- nrow(A)
  n_tau <- length(taus_mha)

  # Identify zeroish values (≤ zero_tol_ha)
  is_zero <- A <= zero_tol_mha
  has_zero <- rowSums(is_zero, na.rm = TRUE) > 0
  all_nonzero <- rowSums(is_zero, na.rm = TRUE) == 0  # no dataset ≤ threshold ha

  # Identify cells where everything is essentially zero
  has_nonzero <- rowSums(A > zero_tol_mha, na.rm = TRUE) > 0
  all_zeroish <- !has_nonzero

  # Disagreement matrix: cells × tau
  disagree_mat <- matrix(FALSE, n_cells, n_tau)

  for (i in seq_along(taus_mha)) {

    tau <- taus_mha[i]

    # Datasets ≥ tau
    is_large <- A >= tau
    has_large <- rowSums(is_large, na.rm = TRUE) > 0

    # Existential disagreement:
    # at least one ≤ threshold ha AND at least one ≥ tau
    disagree_mat[, i] <- has_zero & has_large
  }

  # Largest tau where existential disagreement exists
  idx_max <- apply(disagree_mat, 1, function(x) {
    if (!any(x)) return(NA_integer_)
    max(which(x))
  })

  # Force NA for structurally consistent cells
  idx_max[all_nonzero] <- NA_integer_  # universally irrigated
  idx_max[all_zeroish] <- NA_integer_  # universally non-irrigated

  data.table(lon = dt$lon,
             lat = dt$lat,
             tau_max_ha = taus_ha[idx_max],
             tau_max_mha = taus_mha[idx_max])
}
