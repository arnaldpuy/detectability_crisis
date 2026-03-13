
# FUNCTION TO TRANSLATE LON LAT TO COUNTRY NAMES ###############################

lonlat_to_country <- function(lon, lat, countries_sf) {
  pts <- st_as_sf(data.frame(lon = lon, lat = lat),
                  coords = c("lon", "lat"),
                  crs = 4326,
                  remove = FALSE)

  joined <- st_join(pts, countries_sf, join = st_within, left = TRUE)

  joined$name
}
