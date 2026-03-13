


# Funciton to reshape wide + run compute_tau_fun -------------------------------

run_tau_pipeline_fun <- function(path, tau_grid) {
  dt <- fread(path)

  dt_wide <- dcast(dt, lon + lat + country + code + continent ~ dataset,
                   value.var = "mha")

  meta_cols <- c("lon", "lat", "country", "code", "continent")
  dataset_names <- setdiff(names(dt_wide), meta_cols)

  tau_grid[
    ,
    compute_tau_fun(dt_wide, tau_mha, dataset_names),
    by = .(tau_label, tau_ha, tau_mha)
  ]
}
