
# FUNCTION TO APPLY CROPLAND MASK ##############################################

filter_long_dt_by_crop_mask_aggregated <- function(dt,
                                                   crop_native,
                                                   rule = c("fraction", "any"),
                                                   threshold = 0.01) {

  rule <- match.arg(rule)
  dt   <- data.table::copy(dt)

  # binary crop mask: crop = 1, non-crop = 0
  crop_bin <- crop_native
  v <- terra::values(crop_bin, mat = FALSE)
  v <- ifelse(!is.na(v) & v > 0, 1, 0)
  terra::values(crop_bin) <- v

  res_lookup <- c("0.2deg" = 0.2,
                  "0.4deg" = 0.4,
                  "1deg"   = 1.0)

  if (!all(unique(dt$resolution) %in% names(res_lookup))) {
    stop("Unknown resolution labels: ",
         paste(setdiff(unique(dt$resolution), names(res_lookup)), collapse = ", "))
  }

  out <- vector("list", length(unique(dt$resolution)))
  k <- 1

  for (rr in unique(dt$resolution)) {

    dt_sub <- dt[resolution == rr]
    target_res <- unname(res_lookup[rr])

    fact_x <- round(target_res / terra::res(crop_bin)[1])
    fact_y <- round(target_res / terra::res(crop_bin)[2])

    # aggregate native mask to target resolution
    if (rule == "any") {
      crop_aggr <- terra::aggregate(
        crop_bin,
        fact = c(fact_x, fact_y),
        fun = max,
        na.rm = TRUE
      )
    } else {
      crop_aggr <- terra::aggregate(
        crop_bin,
        fact = c(fact_x, fact_y),
        fun = mean,
        na.rm = TRUE
      )
    }

    # find which aggregated raster cell each dt row falls into
    cell_id <- terra::cellFromXY(crop_aggr, as.matrix(dt_sub[, .(lon, lat)]))

    # extract aggregated crop values for those cells
    vals <- rep(NA_real_, length(cell_id))
    ok   <- !is.na(cell_id)
    vals[ok] <- terra::values(crop_aggr, mat = FALSE)[cell_id[ok]]

    # apply filter
    if (rule == "any") {
      keep <- !is.na(vals) & vals > 0
    } else {
      keep <- !is.na(vals) & vals >= threshold
    }

    out[[k]] <- dt_sub[keep]
    k <- k + 1
  }

  rbindlist(out, use.names = TRUE)
}
