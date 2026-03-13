
# FUNCTIONS TO REGRID IRRIGATED AREA MAPS ######################################

# Helpers ----------------------------------------------------------------------

assert_unique_xy <- function(dti, dset) {
  dup <- dti[, .N, by = .(lon, lat)][N > 1L]
  if (nrow(dup) > 0L) {
    stop(sprintf(
      "Dataset '%s' has duplicate lon/lat points (n=%d). Aggregate intentionally first.",
      dset, nrow(dup)
    ))
  }
  invisible(TRUE)
}

clamp01 <- function(r) clamp(r, lower = 0, upper = 1, values = TRUE)

# Explicitly aligned lon/lat templates -----------------------------------------

make_template_ll <- function(res_deg, extent_ll = ext(-180, 180, -90, 90)) {
  rast(extent_ll, crs = "EPSG:4326", resolution = c(res_deg, res_deg))
}

# MAKE THE 0.2º RASTER #########################################################

make_raster_ll_02 <- function(dset) {
  dti <- dt[dataset == dset, .(lon, lat, mha)]

  # Enforce one value per cell-center (or aggregate intentionally)
  assert_unique_xy(dti, dset)
  v <- vect(dti[, .(lon, lat)], crs = "EPSG:4326")
  v$mha <- dti$mha
  r <- rasterize(v, template_ll_02, field = "mha", fun = "sum", background = 0)
  names(r) <- dset
  r
}

# ADD COUNTRIES AND CONTINENT ##################################################

add_admin <- function(d) {
  chunk_size <- 200000L
  idx <- split(seq_len(nrow(d)), ceiling(seq_len(nrow(d)) / chunk_size))
  d[, country:= NA_character_]
  for (i in idx) {
    d[i, country:= lonlat_to_country(lon, lat, countries_sf)]
  }

  country_to_continent(d, country_col = "country",
                       standardize_country_name = TRUE)

  d <- d[!is.na(country)]
  d
}

# DROP ALL ZEROES OR NA CELLS ##################################################

drop_all_zero_or_na_cells <- function(d) {
  tmp <- dcast(d, lon + lat + continent + country + code ~ dataset, value.var = "mha")

  meta_cols <- c("lon", "lat", "country", "code", "continent")
  dataset_cols <- setdiff(names(tmp), meta_cols)

  tmp[, keep := rowSums(!is.na(.SD) & (.SD != 0)) > 0L, .SDcols = dataset_cols]
  tmp2 <- tmp[keep == TRUE][, keep := NULL]

  melt(tmp2, id.vars = meta_cols, measure.vars = dataset_cols,
       variable.name = "dataset", value.name = "mha",
       variable.factor = FALSE)
}








# ------------------------------------------------------------------------------
# QA: THRESHOLD STABILITY (flip rates + prevalence change)
# ------------------------------------------------------------------------------
# Why this QA exists:
#   Your downstream analysis thresholds irrigated area at detectability threshold τ
#   (e.g., "irrigation exists if A > τ").
#   If harmonization changes which cells cross τ (even by tiny numerical amounts),
#   it can artificially increase or decrease disagreement.
#
# Therefore, we quantify:
#   - flip_up: fraction of cells that were <= τ before, but > τ after (absent->present)
#   - flip_down: fraction of cells that were > τ before, but <= τ after (present->absent)
#   - prevalence: fraction of cells with A > τ, before vs after
#
# Inputs:
#   mha_before / mha_after: SpatRaster stacks with identical layers (datasets) and identical grid.
#   Values are irrigated area in Mha (million hectares).
#   tau_ha is in hectares and converted to Mha internally.

tau_flip_rates <- function(mha_before, mha_after, tau_ha = c(0, 100, 1000)) {

  # Defensive checks: we only compare like-for-like rasters
  stopifnot(inherits(mha_before, "SpatRaster"), inherits(mha_after, "SpatRaster"))
  stopifnot(nlyr(mha_before) == nlyr(mha_after))
  stopifnot(all(names(mha_before) == names(mha_after)))

  # Loop over datasets (layers) and thresholds τ
  rbindlist(lapply(seq_len(nlyr(mha_before)), function(i) {

    # Extract the i-th dataset layer before/after harmonization
    b <- mha_before[[i]]
    a <- mha_after[[i]]

    rbindlist(lapply(tau_ha, function(tau) {

      # Convert τ from hectares to Mha (since rasters are in Mha)
      thr <- tau / 1e6

      # Presence/absence classification:
      #   TRUE if irrigated area exceeds detectability threshold τ
      bb <- b > thr
      aa <- a > thr

      # Define valid cells:
      #   We only compare cells where both before and after are not NA
      valid <- !is.na(bb) & !is.na(aa)
      n <- global(valid, "sum", na.rm = TRUE)[1, 1]
      # Why global(valid, "sum"):
      #   counts number of valid cells; we use it as denominator for flip rates.

      # Flip-up: absent->present due to harmonization
      flip_up   <- global(valid & (!bb) & aa, "sum", na.rm = TRUE)[1, 1] / n

      # Flip-down: present->absent due to harmonization
      flip_down <- global(valid & bb & (!aa), "sum", na.rm = TRUE)[1, 1] / n

      # Return a tidy row for this dataset and threshold
      data.table(dataset = names(mha_before)[i],
                 tau_ha = tau,
                 flip_up = flip_up,
                 flip_down = flip_down)
    }))
  }))[order(dataset, tau_ha)]
}


tau_prevalence_change <- function(mha_before, mha_after, tau_ha = c(0, 100, 1000)) {

  # Same defensive checks as flip rates
  stopifnot(inherits(mha_before, "SpatRaster"), inherits(mha_after, "SpatRaster"))
  stopifnot(nlyr(mha_before) == nlyr(mha_after))
  stopifnot(all(names(mha_before) == names(mha_after)))

  rbindlist(lapply(seq_len(nlyr(mha_before)), function(i) {
    b <- mha_before[[i]]
    a <- mha_after[[i]]

    rbindlist(lapply(tau_ha, function(tau) {
      thr <- tau / 1e6  # ha -> Mha

      # Valid where both rasters have data
      valid <- !is.na(b) & !is.na(a)
      n <- global(valid, "sum", na.rm = TRUE)[1, 1]

      # Prevalence = fraction of valid cells classified as present
      p_before <- global(valid & (b > thr), "sum", na.rm = TRUE)[1, 1] / n
      p_after  <- global(valid & (a > thr), "sum", na.rm = TRUE)[1, 1] / n

      data.table(dataset = names(mha_before)[i],
                 tau_ha = tau,
                 p_before = p_before,
                 p_after = p_after,
                 delta = p_after - p_before)
      # Why delta:
      #   large positive delta means harmonization increases irrigation footprint
      #   at that threshold; large negative delta means it removes irrigation.
    }))
  }))[order(dataset, tau_ha)]
}

flip_summary <- function(flip_dt) {
  # Summarize flip-rate distributions across datasets at each τ.
  # Why:
  #   The per-dataset flips can be noisy; the median/min/max are SI-friendly summaries.
  flip_dt[, .(
    flip_up_median = median(flip_up),
    flip_up_min = min(flip_up),
    flip_up_max = max(flip_up),
    flip_down_median = median(flip_down),
    flip_down_min = min(flip_down),
    flip_down_max = max(flip_down)
  ), by = tau_ha][order(tau_ha)]
}


# ------------------------------------------------------------------------------
# QA: OPERATOR SENSITIVITY (bilinear vs near) + OPTIONAL RENORMALIZATION
# ------------------------------------------------------------------------------
# Why:
#   Reprojection/resampling requires an interpolation operator.
#   - Bilinear smooths values and can spread tiny positives into neighboring zeros,
#     which matters a lot at low τ (especially τ=0).
#   - Nearest neighbour preserves categorical structure better (no smoothing),
#     making it preferable for threshold-based analyses.
#
# We therefore run BOTH operators and quantify:
#   - global-sum conservation before/after renormalization
#   - threshold flips and prevalence changes
#   - changes in disagreement curves vs baseline


run_old_roundtrip <- function(f_ll02, template_ea, template_ll_02,
                              method = c("bilinear", "near"),
                              clamp_to_01 = TRUE) {

  method <- match.arg(method)

  # Start with fractions on 0.2° lat/lon grid
  f0 <- f_ll02

  # Clamp to [0,1] to enforce valid fractional coverage.
  # Why:
  #   Projection/interpolation can create tiny overshoots (<0 or >1) due to numerical effects.
  if (clamp_to_01) f0 <- clamp01(f0)

  # Step 1: project from geographic lat/lon to equal-area grid.
  # Why:
  #   Doing resampling in an equal-area CRS helps ensure that "area" concepts are consistent
  #   and avoids distortions associated with varying cell area in geographic coordinates.
  f_ea <- project(f0, template_ea, method = method)
  if (clamp_to_01) f_ea <- clamp01(f_ea)

  # Step 2: project back from equal-area to the target 0.2° lat/lon grid.
  # Why:
  #   Our analysis is conducted on regular lat/lon grids (0.2°, then aggregated).
  f_back <- project(f_ea, template_ll_02, method = method)
  if (clamp_to_01) f_back <- clamp01(f_back)

  # Convert fraction back to irrigated area (Mha).
  # Why:
  #   The analysis uses irrigated area per cell (Mha), not fractions.
  A_ll02_back_ha <- cellSize(f_back, unit = "ha")       # cell area in hectares
  mha_back <- (f_back * A_ll02_back_ha) / 1e6           # ha -> Mha

  # Preserve layer names to keep dataset identities aligned
  names(mha_back) <- names(f_ll02)

  list(f_ea = f_ea, f_back = f_back, mha = mha_back)
}


renormalize_to_totals <- function(mha_target, mha_reference) {
  # Why:
  #   Even with careful processing, reprojection/resampling can lose a small amount of mass.
  #   To pre-empt reviewer criticism ("your totals changed"), we enforce exact conservation.
  #
  # Mechanism:
  #   For each dataset layer i, compute scale_i = total_reference_i / total_target_i
  #   and multiply the target layer by scale_i.

  ref <- global(mha_reference, "sum", na.rm = TRUE)[, 1]
  tgt <- global(mha_target, "sum", na.rm = TRUE)[, 1]
  scale <- ref / tgt

  out <- mha_target
  for (i in seq_len(nlyr(out))) out[[i]] <- out[[i]] * scale[i]
  names(out) <- names(mha_target)
  out
}


# Convert a multi-layer raster stack into a wide data.table (lon, lat, dataset columns)
stack_to_dtwide <- function(r) {
  as.data.table(as.data.frame(r, xy = TRUE, na.rm = FALSE)) |>
    setnames(c("x", "y"), c("lon", "lat"))
}

# Compute the full τ curve for a raster stack
compute_curve_from_stack <- function(r, tau_grid) {
  dtw <- stack_to_dtwide(r)
  dataset_names <- setdiff(names(dtw), c("lon", "lat"))

  # compute_tau_fun() returns metrics (frac_disagree, etc.) at a given tau_mha
  out <- tau_grid[
    ,
    compute_tau_fun(dtw, tau_mha, dataset_names),
    by = .(tau_label, tau_ha, tau_mha)
  ]
  out[]
}

# Summarize: worst-case deviation between two curves for one metric
summarize_delta <- function(delta_dt, metric) {
  stopifnot(metric %in% names(delta_dt))
  tmp <- copy(delta_dt)
  tmp[, absd := abs(get(metric))]

  worst <- tmp[which.max(absd),
               .(tau_label, tau_ha, tau_mha,
                 delta = get(metric),
                 abs_delta = absd)]

  data.table(metric = metric,
             max_abs_delta = worst$abs_delta,
             tau_ha_at_max = worst$tau_ha,
             tau_label_at_max = worst$tau_label,
             delta_at_max = worst$delta)
}
