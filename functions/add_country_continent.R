
# FUNCTION TO ADD COUNTRIES AND CONTINENT BASED ON COORDINATES #################

add_country_continent <- function(dt, countries_sf, chunk_size = 200000L) {
  dt[, country := NA_character_]
  idx <- split(seq_len(nrow(dt)),
               ceiling(seq_len(nrow(dt)) / chunk_size))
  for (i in idx) {
    dt[i, country:= lonlat_to_country(lon, lat, countries_sf)]
  }
  country_to_continent(dt, country_col = "country", standardize_country_name = TRUE)
  dt
}
