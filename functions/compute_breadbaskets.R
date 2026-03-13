

# FUNCTION TO COMPUTE BREADBASKET LOSSES #######################################

compute_breadbaskets <- function(res_tag, crops_aligned, dt_pres,
                                 regions, grain_crops, land_mask) {

  dt_pres_res <- dt_pres[resolution == res_tag,
                         .(lon, lat, country, code, continent, n_pos)]

  crops_join <- merge(crops_aligned, dt_pres_res,
                      by = c("lon", "lat", "country", "code", "continent"),
                      all.x = TRUE)

  crops_join[is.na(n_pos), n_pos := 0L]

  if (!"grain_prod" %in% names(crops_join)) {
    crops_join[, grain_prod := rowSums(.SD, na.rm = TRUE), .SDcols = grain_crops]
  }

  coords <- cbind(crops_join$lon, crops_join$lat)
  land_vals <- terra::extract(land_mask, coords)[, 1]
  crops_join[, land:= !is.na(land_vals)]
  crops_join_land <- crops_join[land == TRUE]

  ks <- 10:5

  rbindlist(
    lapply(names(regions), function(region_name) {
      e <- regions[[region_name]]

      dt_region <- crops_join_land[
        lon >= xmin(e) & lon <= xmax(e) &
          lat >= ymin(e) & lat <= ymax(e)
      ]

      totals <- vapply(
        ks,
        function(k) dt_region[n_pos >= k, sum(grain_prod, na.rm = TRUE)],
        numeric(1)
      )
      names(totals) <- as.character(ks)

      base <- totals["5"]

      if (is.na(base) || base <= 0) {
        frac <- rep(NA_real_, length(totals))
        loss <- rep(NA_real_, length(totals))
      } else {
        frac <- totals / base
        loss <- 1 - totals / base
      }

      data.table(region = region_name,
                 resolution = res_tag,
                 k = ks,
                 total_prod = as.numeric(totals),
                 frac_vs_5 = as.numeric(frac),
                 loss_vs_5 = as.numeric(loss))
    }),
    use.names = TRUE,
    fill = TRUE
  )
}
