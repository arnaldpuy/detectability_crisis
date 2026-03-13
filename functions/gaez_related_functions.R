
# GAEZ-RELATED FUNCTIONS (TO REGRID CROP DATA FROM GAEZ) #######################


# Helper to regrid GAEZ to a regular lon/lat grid (sum) ------------------------

# Regrid native 5' GAEZ production to a regular lon/lat grid at resolution `res`.
# Production is extensive, so we use sum.

regrid_gaez_dt <- function(dt, res) {
  dt <- copy(dt)

  # Bin indices (assuming grid is anchored at -180, -90, as in your other work)
  dt[, lon_bin := floor((lon + 180) / res)]
  dt[, lat_bin := floor((lat + 90)  / res)]

  # Coarse grid cell centres
  dt[, lon_c := -180 + (lon_bin + 0.5) * res]
  dt[, lat_c := -90  + (lat_bin + 0.5) * res]

  # Identify crop columns
  helper_cols <- c("lon", "lat", "lon_bin", "lat_bin", "lon_c", "lat_c")
  crop_cols   <- setdiff(names(dt), helper_cols)

  # Sum production within each coarse cell
  dt_out <- dt[
    ,
    lapply(.SD, sum, na.rm = TRUE),
    by = .(lon = lon_c, lat = lat_c),
    .SDcols = crop_cols
  ]

  # Optional: round coords slightly to avoid floating-point mismatch
  dt_out[, `:=`(
    lon = round(lon, 6),
    lat = round(lat, 6)
  )]

  setorder(dt_out, lat, lon)
  dt_out[]
}




# Helper to get template grid cells from a raster-------------------------------

# Extract (lon, lat) grid cell centres from one layer of a template raster
get_template_coords <- function(template_path) {
  r <- rast(template_path)

  # use first layer; we just need coordinates
  coords <- as.data.table(
    as.data.frame(r[[1]], xy = TRUE, cells = FALSE, na.rm = FALSE)
  )

  setnames(coords, c("x", "y"), c("lon", "lat"))

  # Reduce to unique lon/lat grid
  coords <- unique(coords[, .(lon, lat)])

  # Round to match regridded coords rounding
  coords[, `:=`(
    lon = round(lon, 6),
    lat = round(lat, 6)
  )]

  setorder(coords, lat, lon)
  coords[]
}





# 3. Full function: regrid + align to a given template -------------------------

# Regrid GAEZ to resolution `res`, then align to `template_path` grid cells.
# Returns a data.table with:
#   lon, lat + one column per crop
# for exactly the same (lon, lat) as the template raster.
regrid_gaez_to_template <- function(dt_native, res, template_path) {
  # 1) Aggregate 5' -> coarse grid by sum
  dt_aggr <- regrid_gaez_dt(dt_native, res = res)

  # 2) Get template coordinates
  coords_template <- get_template_coords(template_path)

  # 3) Merge template coords with aggregated GAEZ
  dt_merged <- merge(
    coords_template,
    dt_aggr,
    by = c("lon", "lat"),
    all.x = TRUE,
    sort  = FALSE
  )

  # 4) Replace NA crops with 0 (no production in that coarse cell)
  crop_cols <- setdiff(names(dt_merged), c("lon", "lat"))
  for (j in crop_cols) {
    set(dt_merged, which(is.na(dt_merged[[j]])), j, 0)
  }

  setorder(dt_merged, lat, lon)
  dt_merged[]
}
