
# FUNCTIONS FOR THE EVAPOTRANSPIRATION WORKFLOW ################################

# HELPERS ----------------------------------------------------------------------

# Read CESM timmean (extensionless files)---------------------------------------

read_cesm_timmean <- function(fpath, out_name = NULL, to_mmday = TRUE) {
  r <- terra::rast(fpath)
  r <- terra::rotate(r)
  if (to_mmday) r <- r * 86400 # Unit conversion from mm/s to mm/d)
  if (!is.null(out_name)) names(r) <- out_name
  r
}

# Resample a raster to each template--------------------------------------------

resample_to_templates <- function(r, templates, method = "bilinear") {
  lapply(templates, function(tmpl) terra::resample(r, tmpl, method = method))
}

# Convert raster to DT using TEMPLATE cell ids (no float joins)
# - gives lon/lat from the template grid (xyFromCell)
# - cell ids are stable within each resolution ---------------------------------

rast_to_dt_with_cell <- function(rr, tmpl, value_name = names(rr)[1]) {
  stopifnot(terra::ncell(rr) == terra::ncell(tmpl))

  vals <- terra::values(rr, mat = FALSE)
  cell <- seq_len(terra::ncell(rr))
  xy   <- terra::xyFromCell(tmpl, cell)

  dt <- data.table(cell = cell, lon  = xy[, 1], lat  = xy[, 2])
  dt[, (value_name):= vals]
  dt[]
}

# Convert dt_pres (lon/lat) to template cell ids using cellFromXY (nearest)
# This avoids float equality issues entirely.-----------------------------------

pres_to_cell <- function(dt_pres_res, tmpl) {
  dt <- copy(dt_pres_res)

  # Map each row to a template cell.
  dt[, cell := terra::cellFromXY(tmpl, cbind(lon, lat))]

  # Drop rows that failed mapping (should be 0)
  dt <- dt[!is.na(cell)]

  # If dt_pres has repeated lon/lat/cell due to extra keys (country, etc.),
  # collapse to one row per cell for the join (keep n_pos and dataset fractions).
  # Choose max n_pos per cell (conservative: "cell has irrigation if any record says so").
  dataset_cols <- setdiff(names(dt), c("country","code","continent","resolution","lon","lat"))

  num_cols <- dataset_cols[sapply(dt[, ..dataset_cols], is.numeric)]

  dt_cell <- dt[, lapply(.SD, max, na.rm = TRUE), by = .(cell), .SDcols = num_cols]
  dt_cell[]
}


# Pairwise ET (IRRk - NOIk) on each template grid -----------------------------

delta_pair_resampled <- function(irr_f, noi_f, templates) {
  tag_irr <- sub(".*_(IRR0[1-3]).*", "\\1", basename(irr_f))
  tag_noi <- sub("IRR", "NOI", tag_irr)

  irr <- read_cesm_timmean(irr_f, out_name = paste0("ET_", tag_irr, "_mmday"))
  noi <- read_cesm_timmean(noi_f, out_name = paste0("ET_", tag_noi, "_mmday"))

  irr_rs <- resample_to_templates(irr, templates, method = "bilinear")
  noi_rs <- resample_to_templates(noi, templates, method = "bilinear")

  mapply(function(a, b) {
    d <- a - b
    names(d) <- paste0("dET_", tag_irr, "_mmday")
    d
  }, irr_rs, noi_rs, SIMPLIFY = FALSE)
}

# area per cell from a template (m^2) ------------------------------------------

cell_area_dt <- function(tmpl) {
  a <- terra::cellSize(tmpl, unit = "m")   # m^2 for lon/lat grids
  data.table(
    cell = 1:terra::ncell(tmpl),
    area_m2 = as.vector(a)
  )
}
