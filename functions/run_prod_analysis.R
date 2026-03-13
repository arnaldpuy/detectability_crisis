

################################################################################
## FUNCTION: production under agreement rule for any resolution
################################################################################

run_prod_analysis <- function(res_tag, crops_aligned, dt_pres, crop_cols) {

  dt_pres_res <- dt_pres[resolution == res_tag,
                         .(lon, lat, country, code, continent, n_pos)]

  crops_join <- merge(
    crops_aligned,
    dt_pres_res,
    by = c("lon", "lat", "country", "code", "continent"),
    all.x = TRUE
  )

  crops_join[is.na(n_pos), n_pos := 0L]

  ks <- 10:5

  ## ---- Per-crop totals ----
  prod_agreement_wide <- rbindlist(
    lapply(ks, function(k) {
      crops_join[n_pos >= k,
                 lapply(.SD, sum, na.rm = TRUE),
                 .SDcols = crop_cols][, k := k][]
    }),
    use.names = TRUE,
    fill = TRUE
  )

  prod_agreement_long <- melt(
    prod_agreement_wide,
    id.vars = "k",
    variable.name = "crop",
    value.name = "production"
  )

  baseline_dt <- prod_agreement_long[
    k == 5,
    .(crop, baseline_prod = production)
  ]

  prod_agreement_long <- merge(
    prod_agreement_long,
    baseline_dt,
    by = "crop",
    all.x = TRUE
  )

  prod_agreement_long[
    , `:=`(
      frac_vs_5 = fifelse(baseline_prod > 0,
                          production / baseline_prod,
                          NA_real_),
      loss_vs_5 = fifelse(baseline_prod > 0,
                          1 - production / baseline_prod,
                          NA_real_)
    )
  ]

  prod_agreement_long[, resolution := res_tag]

  ## ---- Global production ----
  crops_join[
    , total_crop_prod := rowSums(.SD, na.rm = TRUE),
    .SDcols = crop_cols
  ]

  prod_total <- rbindlist(
    lapply(ks, function(k) {
      crops_join[n_pos >= k,
                 .(k = k,
                   total_production = sum(total_crop_prod, na.rm = TRUE))]
    })
  )

  baseline_total <- prod_total[k == 5, total_production][1]

  prod_total[
    , `:=`(
      frac_vs_5 = total_production / baseline_total,
      loss_vs_5 = 1 - total_production / baseline_total
    )
  ]
  prod_total[, resolution := res_tag]

  list(crop = prod_agreement_long,
       total = prod_total)
}
