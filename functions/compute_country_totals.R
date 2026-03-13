
# FUNCTION TO COMPUTE COUNTRY TOTALS ###########################################

compute_country_totals <- function(cell_dt, dataset_names) {
  rbindlist(lapply(dataset_names, function(d) {
    cell_dt[, .(irrigated_area_mha = sum(get(d), na.rm = TRUE)), country] %>%
      .[, `:=`(
        dataset = d,
        rank = frank(-irrigated_area_mha, ties.method = "average")
      )]
  }))
}
