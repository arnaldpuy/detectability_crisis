
# PRELIMINARY FUNCTIONS ########################################################
################################################################################

sensobol::load_packages(c("data.table", "terra", "here", "countrycode",
                          "rworldmap", "sp", "sf", "rnaturalearth", "readxl",
                          "tidyverse", "scales", "ncdf4"))

# Source all helper functions --------------------------------------------------

r_functions <- list.files(path = here("functions"), pattern = "\\.R$",
                          full.names = TRUE)

invisible(lapply(r_functions, source))

# Reproducibility --------------------------------------------------------------

set.seed(123L)

# Countries layer (WGS84, only name + geometry) --------------------------------

countries_sf <- rnaturalearth::ne_countries(scale = "medium",
                                            returnclass = "sf") %>%
  .[, c("name", "geometry")]

# Output folder for all original data tables -----------------------------------

out_dir <- here("original_datasets")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Columns to keep for gridded products -----------------------------------------

cols_to_retrieve <- c("lon", "lat", "country", "continent", "mha", "dataset")

# =============================================================================
# MEIER ET AL.
#   - 1 km grid; keep class 1 pixels and convert to Mha per cell
# =============================================================================

f_meier <- "./irrigated_area_datasets/meier_et_al/global_irrigated_areas.tif"
r_meier <- rast(f_meier)

meier_dt <- as.data.table(as.data.frame(r_meier, xy = TRUE, na.rm = TRUE))
setnames(meier_dt, old = c("x", "y", "global_irrigated_areas"),
         new = c("lon", "lat", "class"))

# Keep irrigated areas ---------------------------------------------------------
# LEGEND:
# 0 = no irrigated area",
# 1 = downscaled Siebert et al. 2013",
# 2 = low suitability + high NDVI + NDVI vegetation course",
# 3 = potential multiple cropping < actual multiple cropping",
# 4 = cropland (ESA-CCI-LC/GlobCover) + low suitability",

meier_dt <- meier_dt[class %in% 1:4]

# Pixel area at 30 arc-sec (1/120 deg) ----------------------------------------

dlat <- 1 / 120
dlon <- 1 / 120
R <- 6371008.8  # m

meier_dt[, cell_area_m2:= (R^2) * (dlon * pi / 180) * (sin((lat + dlat / 2) * pi / 180) -
              sin((lat - dlat / 2) * pi / 180))]
meier_dt[, cell_area_m2:= abs(cell_area_m2)]
meier_dt[, cell_area_ha:= cell_area_m2 / 1e4]
meier_dt[, cell_area_mha:= cell_area_ha / 1e6]

# Countries + continents -------------------------------------------------------

meier_dt <- add_country_continent(meier_dt, countries_sf, chunk_size = 500000L)

# Dataset + mha column ---------------------------------------------------------

meier_dt[, dataset:= "meier"]
meier_dt[, mha:= cell_area_mha]

meier_dt <- meier_dt[, ..cols_to_retrieve]
meier_dt <- na.omit(meier_dt)

# Check totals------------------------------------------------------------------

# Original paper reports 367 Mha (p.1124)
meier_dt[, sum(mha)]

# Export------------------------------------------------------------------------

fwrite(meier_dt, file.path(out_dir, "meier_dt.csv"))

# =============================================================================
# GIAM (Thenkabail et al. 2009)
#   - 30 arc-sec grid; class-based IAFs → TAAI + AIA, we use TAAI Mha
# =============================================================================

f_giam <- "./irrigated_area_datasets/giam/giam_28_classes_global.tif"
r_giam <- rast(f_giam)

giam_dt <- as.data.table(as.data.frame(r_giam, xy = TRUE, na.rm = TRUE))

# Drop non-irrigated class 0 and ensure integer class codes --------------------

giam_dt <- giam_dt[giam_28_classes_global != 0]
giam_dt[, giam_28_classes_global:= as.integer(giam_28_classes_global)]
setnames(giam_dt, c("x", "y"), c("lon", "lat"))

# Pixel area (30 arc-sec) ------------------------------------------------------

dlat <- 1 / 120
dlon <- 1 / 120
R <- 6371008.8

giam_dt[, cell_area_m2:= abs((R^2) * (dlon * pi / 180) *
                                (sin((lat + dlat / 2) * pi / 180) -
                                   sin((lat - dlat / 2) * pi / 180)))]
giam_dt[, cell_area_ha:= cell_area_m2 / 1e4]
giam_dt[, cell_area_mha:= cell_area_m2 / 1e10]

# IAF table (class → seasonal fractions) ---------------------------------------

iaf_xlsx <- "./irrigated_area_datasets/giam/GIAM_IAF_Table.xlsx"
iaf <- as.data.table(readxl::read_excel(iaf_xlsx, sheet = "GIAM_IAF_Table"))
setnames(iaf, tolower(names(iaf)))

iaf[, s1_mean:= rowMeans(.SD, na.rm = TRUE), .SDcols = c("s1_hri", "s1_spdt")]
iaf[, s2_mean:= rowMeans(.SD, na.rm = TRUE), .SDcols = c("s2_hri", "s2_spdt")]
iaf[, cont_mean:= rowMeans(.SD, na.rm = TRUE), .SDcols = c("cont_hri", "cont_spdt")]

for (cc in c("s1_mean", "s2_mean", "cont_mean")) {
  iaf[is.nan(get(cc)), (cc) := NA_real_]
}

setkey(iaf, class)
setkey(giam_dt, giam_28_classes_global)

giam_dt[iaf, `:=`(iaf_gee = i.iaf_gee, iaf_s1 = i.s1_mean,
                  iaf_s2 = i.s2_mean, iaf_cont = i.cont_mean)]

# Net irrigated area (TAAI) + annualized area (AIA) ----------------------------

giam_dt[, irrig_taai_mha:= cell_area_mha * iaf_gee]

giam_dt[, irrig_aia_mha:= cell_area_mha * (fcoalesce(iaf_s1, 0) + fcoalesce(iaf_s2, 0) +
                                             fcoalesce(iaf_cont, 0))]
giam_dt[, irrig_taai_ha:= irrig_taai_mha * 1e6]
giam_dt[, irrig_aia_ha:= irrig_aia_mha  * 1e6]

# Countries + continents -------------------------------------------------------

giam_dt <- add_country_continent(giam_dt, countries_sf, chunk_size = 500000L)

# Dataset + mha (use TAAI) -----------------------------------------------------

giam_dt <- giam_dt[!is.na(country) & !is.na(continent)]
giam_dt[, dataset:= "giam"]
giam_dt[, mha:= irrig_taai_mha]

giam_dt <- giam_dt[, ..cols_to_retrieve]

# Check totals------------------------------------------------------------------

# Original paper reports 398 Mha (Table 3)
giam_dt[, sum(mha)]

# Export------------------------------------------------------------------------

fwrite(giam_dt, file.path(out_dir, "giam_dt.csv"))

# =============================================================================
# GMIA v5.0 (Siebert et al.)
#   - 5 arc-min grid; AEI × AAI% → actually irrigated Mha per cell
# =============================================================================

r_gmia_aei <- rast("./irrigated_area_datasets/gmia/gmia_v5_aei_ha.asc")
gmia_dt <- as.data.table(as.data.frame(r_gmia_aei, xy = TRUE, na.rm = FALSE))

r_gmia_pct <- rast("./irrigated_area_datasets/gmia/gmia_v5_aai_pct_aei.asc")
gmia_dt_pct <- as.data.table(as.data.frame(r_gmia_pct, xy = TRUE, na.rm = FALSE))

setkey(gmia_dt, x, y)
setkey(gmia_dt_pct, x, y)
gmia_dt <- gmia_dt[gmia_dt_pct]

# Remove non-irrigated cells and convert percentage to fraction ----------------

gmia_dt[, gmia_v5_aai_pct_aei:= gmia_v5_aai_pct_aei / 100]
gmia_dt[, gmia_v5_aei_mha:= gmia_v5_aei_ha / 1e6]

setnames(gmia_dt, c("x", "y"), c("lon", "lat"))

# Countries + continents -------------------------------------------------------

gmia_dt <- add_country_continent(gmia_dt, countries_sf, chunk_size = 200000L)
gmia_dt <- na.omit(gmia_dt)

# Dataset + mha (AEI × AAI fraction) ------------------------------------------

gmia_dt[, dataset:= "gmia"]
gmia_dt[, mha:= gmia_v5_aei_mha * gmia_v5_aai_pct_aei]

gmia_dt <- gmia_dt[, ..cols_to_retrieve]

# Check totals------------------------------------------------------------------

# Original paper reports 255 Mha (p. 14)
gmia_dt[, sum(mha)]

# Export------------------------------------------------------------------------

fwrite(gmia_dt, file.path(out_dir, "gmia_dt.csv"))

# =============================================================================
# NAGARAJ ET AL. (2021)
#   - 5 arc-min; classes → representative fractions × fixed 8604 ha per pixel
# =============================================================================

r_nag <- rast("./irrigated_area_datasets/nagaraj_et_al/v3b_combined_2015.tif")
nagaraj_dt <- as.data.table(as.data.frame(r_nag, xy = TRUE, na.rm = FALSE))

# Keep irrigation classes only -------------------------------------------------
nagaraj_dt <- nagaraj_dt[classification %in% c(1, 2)]

pixel_area_ha <- 8604  # fixed by paper

# Irrigation classes represent the fraction of a 5-arc-minute pixel that is
# irrigated. Pixels classified as low–medium irrigation correspond to approximately
# 1–20% of the pixel being irrigated.Pixels classified as high-intensity irrigation
# correspond to more than 20% of the pixel being irrigated.
# To estimate irrigated area, we assign representative irrigated fractions
# of 10% to the low–medium irrigation class and 35% to the high-intensity class.
# Irrigated area is then computed by multiplying these fractions by the
# nominal 5-arc-minute pixel area (8604 ha).
# This procedure is consistent with the class definitions and training
# thresholds described in Nagaraj et al. (2021).

nagaraj_dt[, irrig_frac:= fifelse(
  classification == 1, 0.10, # low–medium: 1–20%
  fifelse(classification == 2, 0.35, NA_real_)  # high: >20%
)]

nagaraj_dt[, irrig_area_ha:= irrig_frac * pixel_area_ha]
nagaraj_dt[, irrig_area_mha:= irrig_area_ha / 1e6]

setnames(nagaraj_dt, c("x", "y"), c("lon", "lat"))

# Countries + continents -------------------------------------------------------

nagaraj_dt <- add_country_continent(nagaraj_dt, countries_sf, chunk_size = 200000L)

nagaraj_dt[, dataset:= "nagaraj"]
nagaraj_dt[, mha:= irrig_area_mha]
nagaraj_dt <- na.omit(nagaraj_dt)

nagaraj_dt <- nagaraj_dt[, ..cols_to_retrieve]

# Check totals------------------------------------------------------------------

# Original paper does not report a total
nagaraj_dt[, sum(mha)]

# Export------------------------------------------------------------------------

fwrite(nagaraj_dt, file.path(out_dir, "nagaraj_dt.csv"))

# =============================================================================
# MIRCA 2000
#   - Maximum area equipped for irrigation per cell (ha → Mha)
# =============================================================================

r_mirca <- rast("./irrigated_area_datasets/MIRCA_2000/MAX_CROPPED_AREA_IRC_HA.ASC")
mirca_dt <- as.data.table(as.data.frame(r_mirca, xy = TRUE, na.rm = FALSE))

mirca_dt <- mirca_dt[MAX_CROPPED_AREA_IRC_HA != 0]
mirca_dt[, mha := MAX_CROPPED_AREA_IRC_HA / 1e6]

setnames(mirca_dt, c("x", "y"), c("lon", "lat"))

# Countries + continents -------------------------------------------------------

mirca_dt <- add_country_continent(mirca_dt, countries_sf, chunk_size = 500000L)

mirca_dt[, dataset:= "mirca2000"]
mirca_dt <- mirca_dt[, ..cols_to_retrieve]
mirca_dt <- mirca_dt[!is.na(country)]

# Check totals------------------------------------------------------------------

# Original paper does not report a total
mirca_dt[, sum(mha)]

# Export------------------------------------------------------------------------

fwrite(mirca_dt, file.path(out_dir, "mirca_dt.csv"))

# =============================================================================
# GAEZ+ 2015
#   - Summed irrigated harvested area per cell, capped by physical cell area
# =============================================================================

gaez_files <- list.files("./irrigated_area_datasets/GAEZ+_2015",
                         pattern = "\\.tif$", full.names = TRUE)
r_gaez  <- rast(gaez_files)
gaez_dt <- as.data.table(as.data.frame(r_gaez, xy = TRUE, na.rm = FALSE))

crop_cols <- grep("^GAEZAct2015_HarvArea_.*_Irrigated$",
                  names(gaez_dt), value = TRUE)

# Total irrigated harvested area per cell (1000 ha) ----------------------------

gaez_dt[, irrig_iharv_1000ha:= rowSums(.SD, na.rm = TRUE), .SDcols = crop_cols]

gaez_dt[, `:=`(irrig_iharv_ha = irrig_iharv_1000ha * 1000,
               irrig_iharv_mha = irrig_iharv_1000ha / 1000)]

# Physical cell area (approx, 0.0833° grid) -----------------------------------

dlat <- 1 / 12
dlon <- 1 / 12
R <- 6371008.8

gaez_dt[, cell_area_ha:= abs((R^2) * (dlon * pi / 180) *
                               (sin((y + dlat / 2) * pi / 180) -
                                  sin((y - dlat / 2) * pi / 180))) / 1e4]

# Cap irrigated area by physical extent ---------------------------------------

gaez_dt[, irrig_extent_ha:= pmin(cell_area_ha, irrig_iharv_ha)]
gaez_dt[, irrig_extent_mha:= irrig_extent_ha / 1e6]

gaez_dt <- gaez_dt[!is.na(irrig_iharv_1000ha)]
setnames(gaez_dt, c("x", "y"), c("lon", "lat"))

# Countries + continents -------------------------------------------------------

gaez_dt <- add_country_continent(gaez_dt, countries_sf, chunk_size = 200000L)
gaez_dt <- gaez_dt[!is.na(country)]

gaez_dt[, dataset:= "gaez_v4"]
gaez_dt[, mha:= irrig_extent_mha]

gaez_dt <- gaez_dt[, ..cols_to_retrieve]

# Check totals------------------------------------------------------------------

# Original paper does not report a total
gaez_dt[, sum(mha)]

# Export------------------------------------------------------------------------

fwrite(gaez_dt, file.path(out_dir, "gaez_dt.csv"))

# =============================================================================
# SPAM 2010
#   - Sum irrigated crop areas per cell, capped by cell area (Mha)
# =============================================================================

spam_files <- list.files("./irrigated_area_datasets/SPAM_2010",pattern = "\\.tif$",
                         full.names = TRUE)
r_spam  <- rast(spam_files)
spam_dt <- as.data.table(as.data.frame(r_spam, xy = TRUE, na.rm = FALSE))

irrig_cols <- grep("^spam2010V2r0_global_A_.*_I$", names(spam_dt), value = TRUE)
spam_dt[, irrig_crop_area_ha:= rowSums(.SD, na.rm = TRUE), .SDcols = irrig_cols]

setnames(spam_dt, c("x", "y"), c("lon", "lat"))

# Cell area on SPAM grid (ha) --------------------------------------------------

template <- r_spam[[1]]
cell_area_ha_r <- cellSize(template, unit = "ha")

spam_dt[, cell:= cellFromXY(template, cbind(lon, lat))]
spam_dt[, cell_area_ha:= values(cell_area_ha_r)[cell]]

spam_dt[, irrig_extent_ha:= pmin(irrig_crop_area_ha, cell_area_ha)]
spam_dt[, irrig_extent_mha:= irrig_extent_ha / 1e6]

# Countries + continents -------------------------------------------------------

spam_dt <- add_country_continent(spam_dt, countries_sf, chunk_size = 200000L)
spam_dt <- spam_dt[!is.na(country)]

spam_dt[, dataset:= "spam"]
spam_dt[, mha:= irrig_extent_mha]

spam_dt <- spam_dt[, ..cols_to_retrieve]

# Check totals------------------------------------------------------------------

# Original paper does not report a total
spam_dt[, sum(mha)]

# Export------------------------------------------------------------------------

fwrite(spam_dt, file.path(out_dir, "spam_dt.csv"))

# =============================================================================
# GRIPC
#   - Irrigated area per cell given directly in ha
# =============================================================================

asc_gripc <- "./irrigated_area_datasets/GRIPC/GRIPC_irrigated_area.asc"
r_gripc <- rast(asc_gripc)

gripc_dt <- as.data.table(as.data.frame(r_gripc, xy = TRUE, na.rm = FALSE))
setnames(gripc_dt, c("x", "y"), c("lon", "lat"))

# Rename raster value to irrig_area_ha -----------------------------------------

setnames(gripc_dt, names(gripc_dt)[3], "irrig_area_ha")

# Cell area (ha) for fraction (optional diagnostic) ---------------------------

cell_ha_r <- cellSize(r_gripc, unit = "ha")
gripc_dt[, cell_area_ha:= values(cell_ha_r)]

gripc_dt[, `:=`(irrig_area_mha = irrig_area_ha / 1e6,
                irrig_frac = irrig_area_ha / cell_area_ha)]

# Keep only irrigated cells ----------------------------------------------------

gripc_dt <- gripc_dt[!is.na(irrig_area_ha)]

# Countries + continents -------------------------------------------------------

gripc_dt <- add_country_continent(gripc_dt, countries_sf, chunk_size = 200000L)
gripc_dt <- gripc_dt[!is.na(country) & !is.na(continent)]

gripc_dt[, dataset:= "gripc"]
gripc_dt[, mha:= irrig_area_mha]

gripc_dt <- gripc_dt[, ..cols_to_retrieve]

# Check totals------------------------------------------------------------------

# Original paper does not report a total
gripc_dt[, sum(mha)]

# Export------------------------------------------------------------------------

fwrite(gripc_dt, file.path(out_dir, "gripc_dt.csv"))

# =============================================================================
# MIRCA-OS (2005)
#   - Sum irrigated crop areas per cell (Mha)
# =============================================================================

mirca_os_path <- "./irrigated_area_datasets/MIRCA_OS/mirca_os_mmgag/Maximum Monthly Growing Area Grids/2005/5-arcminute"
tif_files <- list.files(path = mirca_os_path, pattern = "\\.tif$", full.names = TRUE)

ir_files <- tif_files[grepl("_ir\\.tif$", tif_files)]
ir_stack <- rast(ir_files)

mirca_os_irr_ha <- sum(ir_stack, na.rm = TRUE)
names(mirca_os_irr_ha) <- "irr_ha"

mirca_os_dt <- as.data.table(as.data.frame(mirca_os_irr_ha, xy = TRUE, na.rm = TRUE))
setnames(mirca_os_dt, c("x", "y", "irr_ha"), c("lon", "lat", "irr_ha"))

mirca_os_dt[, mha := irr_ha / 1e6]

# Countries + continents -------------------------------------------------------

mirca_os_dt <- add_country_continent(mirca_os_dt, countries_sf, chunk_size = 200000L)
mirca_os_dt <- mirca_os_dt[!is.na(country) & !is.na(continent)]

mirca_os_dt[, dataset:= "mirca_os"]
mirca_os_dt <- mirca_os_dt[, ..cols_to_retrieve]

# Check totals------------------------------------------------------------------

# Original paper does not report a total
mirca_os_dt[, sum(mha)]

# Export------------------------------------------------------------------------

fwrite(mirca_os_dt, file.path(out_dir, "mirca_os_dt.csv"))

# =============================================================================
# LUH2 (2015 snapshot)
#   - Fraction irrigated crops × cell area to irrigated area per cell (Mha)
# =============================================================================

states_file <- "./irrigated_area_datasets/LUH2/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp119-2-1-f_gn_2015-2100.nc"
mgmt_file <- "./irrigated_area_datasets/LUH2/multiple-management_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp119-2-1-f_gn_2015-2100.nc"

nc <- nc_open(states_file)
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
time <- ncvar_get(nc, "time")
years_all <- 2015 + time
nc_close(nc)

year_to_layer <- function(year) which(years_all == year)

# Land-use fractions -----------------------------------------------------------

c3ann <- rast(states_file, subds = "c3ann")
c4ann <- rast(states_file, subds = "c4ann")
c3per <- rast(states_file, subds = "c3per")
c4per <- rast(states_file, subds = "c4per")
c3nfx <- rast(states_file, subds = "c3nfx")

# Irrigation fractions ---------------------------------------------------------

ir_c3ann <- rast(mgmt_file, subds = "irrig_c3ann")
ir_c4ann <- rast(mgmt_file, subds = "irrig_c4ann")
ir_c3per <- rast(mgmt_file, subds = "irrig_c3per")
ir_c4per <- rast(mgmt_file, subds = "irrig_c4per")
ir_c3nfx <- rast(mgmt_file, subds = "irrig_c3nfx")

k <- year_to_layer(2015)

c3ann_y <- c3ann[[k]]; c4ann_y <- c4ann[[k]]
c3per_y <- c3per[[k]]; c4per_y <- c4per[[k]]
c3nfx_y <- c3nfx[[k]]

ir_c3ann_y <- ir_c3ann[[k]]; ir_c4ann_y <- ir_c4ann[[k]]
ir_c3per_y <- ir_c3per[[k]]; ir_c4per_y <- ir_c4per[[k]]
ir_c3nfx_y <- ir_c3nfx[[k]]

fill <- 1.00000002004088e+20
clamp_fill <- function(r) { r[r >= fill / 10] <- NA; r }

c3ann_y <- clamp_fill(c3ann_y)
c4ann_y <- clamp_fill(c4ann_y)
c3per_y <- clamp_fill(c3per_y)
c4per_y <- clamp_fill(c4per_y)
c3nfx_y <- clamp_fill(c3nfx_y)
ir_c3ann_y <- clamp_fill(ir_c3ann_y)
ir_c4ann_y <- clamp_fill(ir_c4ann_y)
ir_c3per_y <- clamp_fill(ir_c3per_y)
ir_c4per_y <- clamp_fill(ir_c4per_y)
ir_c3nfx_y <- clamp_fill(ir_c3nfx_y)

irr_frac <- c3ann_y * ir_c3ann_y +
  c4ann_y * ir_c4ann_y +
  c3per_y * ir_c3per_y +
  c4per_y * ir_c4per_y +
  c3nfx_y * ir_c3nfx_y

cell_area_ha_luh <- cellSize(irr_frac, unit = "ha")
irr_area_ha_r <- irr_frac * cell_area_ha_luh
names(irr_area_ha_r) <- "irrig_area_ha"

luh2_dt <- as.data.table(as.data.frame(irr_area_ha_r, xy = TRUE, na.rm = FALSE))
setnames(luh2_dt, c("x", "y", "irrig_area_ha"), c("lon", "lat", "irrig_area_ha"))

# Drop NA ----------------------------------------------------------------------

luh2_dt <- luh2_dt[!is.na(irrig_area_ha)]

# Countries + continents (chunked) ---------------------------------------------

luh2_dt <- add_country_continent(luh2_dt, countries_sf, chunk_size = 200000L)
luh2_dt <- luh2_dt[!is.na(country) & !is.na(continent)]

luh2_dt[, mha:= irrig_area_ha / 1e6]
luh2_dt[, dataset:= "luh2"]

luh2_dt <- luh2_dt[, ..cols_to_retrieve]

# Check totals------------------------------------------------------------------

# Original paper does not report a total
luh2_dt[, sum(mha)]

# Export------------------------------------------------------------------------

fwrite(luh2_dt, file.path(out_dir, "luh2_dt.csv"))

meier_dt <- fread("/Users/arnaldpuy/Documents/papers/detectability_irrigated_areas/code_detectability_irrigated_areas/original_datasets/meier_dt.csv")
mirca_dt <- fread("/Users/arnaldpuy/Documents/papers/detectability_irrigated_areas/code_detectability_irrigated_areas/original_datasets/mirca_dt.csv")

# MERGE ALL DATASETS AND EXPORT ################################################

# Rbind ------------------------------------------------------------------------

dada <- rbind(gaez_dt,
              giam_dt,
              gmia_dt,
              gripc_dt,
              luh2_dt,
              meier_dt,
              mirca_dt,
              mirca_os_dt,
              nagaraj_dt,
              spam_dt)

names_datasets <- unique(dada$dataset)
dada <- dada[continent %in% c("Africa", "Americas", "Oceania", "Europe", "Asia")]
fwrite(dada, "./datasets/irrigated_areas/irrigated_areas.csv")

################################################################################
################################################################################

# EXCLUDED DATASETS AFTER PRELIMINARY ASSESSMENTS ##############################

# =============================================================================
# GirrEO
#   - Keep only irrigation within cropland mask
#   - Fraction irrigated × cell area (ha) → irrigated Mha per cell
# =============================================================================

# GirrEO raster ----------------------------------------------------------------

girr_path <- "./irrigated_area_datasets/GirrEO/Assimila_IrrigatedAreas_candidate_EVI_idx_Global_2023_irrigation_mask_0.1666667Deg.tif"
girr_rast <- rast(girr_path)

# Cropland mask ( 1 km, 0.0089286°) --------------------------------------------

# Adjust the path & layer selection to your actual file
crop_path  <- "./irrigated_area_datasets/GirrEO/asap_mask_crop_v02.tif"
crop_native <- rast(crop_path)

# Align extent with GirrEO (y extent: -90 to 90 vs -56 to 75 in crop mask)
# This crops the crop mask to the GirrEO extent (no change in this case for x,
# but keeps things explicit and avoids extent mismatches).
crop_native_aligned <- crop(crop_native, ext(girr_rast))

# Resample cropland mask to GirrEO grid (nearest neighbour for categorical mask)
crop_mask_girr <- resample(crop_native_aligned, girr_rast, method = "near")

# Mask GirrEO irrigation by cropland
# Assumption: 1 = cropland, 0 (or NA) = non-cropland.
# We turn non-cropland (0) into NA so it’s dropped later.
girr_rast_crop <- mask(girr_rast, crop_mask_girr, maskvalues = 0, updatevalue = NA)

# -----------------------------------------------------------------------------
# Convert to data.table
# -----------------------------------------------------------------------------
girr_dt <- as.data.table(as.data.frame(girr_rast_crop, xy = TRUE, na.rm = FALSE))
setnames(girr_dt,
         old = c("x", "y", names(girr_rast_crop)),
         new = c("lon", "lat", "irrigated"))

# Cell area (ha) ---------------------------------------------------------------

cell_area_m2_girr <- cellSize(girr_rast, unit = "m")
cell_area_ha_girr <- cell_area_m2_girr / 1e4

area_dt <- as.data.table(as.data.frame(cell_area_ha_girr, xy = TRUE, na.rm = FALSE))
setnames(area_dt, c("x", "y", names(cell_area_ha_girr)),
         c("lon", "lat", "cell_area_ha"))

setkey(girr_dt, lon, lat)
setkey(area_dt, lon, lat)
girr_dt <- girr_dt[area_dt]

# Treat remaining NAs in 'irrigated' as zero (ocean, non-irrigated cropland, outside mask coverage)
girr_dt[is.na(irrigated), irrigated := 0]

# Irrigated area per cell ------------------------------------------------------

girr_dt[, irr_ha:= irrigated * cell_area_ha]
girr_dt[, irr_mha:= irr_ha / 1e6]

# Countries + continents -------------------------------------------------------

girr_dt <- add_country_continent(girr_dt, countries_sf, chunk_size = 200000L)
girr_dt <- girr_dt[!is.na(country) & !is.na(continent)]

girr_dt[, dataset := "girreo"]
girr_dt[, mha:= irr_mha]

girr_dt <- girr_dt[, ..cols_to_retrieve]
girr_dt[, sum(mha)]

# =============================================================================
# ZOHAIB ET AL.
#   - Binary/class-based mask; irrigated cells take full cell area (Mha)
# =============================================================================

zohaib_path <- "./irrigated_area_datasets/zohaib_et_al/global_actual_irrigated_area.tif"
r_class <- rast(zohaib_path)

# Classes:
# 0 = non-irrigated / masked
# 1 = SM-based, 2 = LST-based, 3 = AL-based, 4 = combined (higher confidence)
irri_classes <- 1:4

r_irri_bin <- app(r_class, fun = function(x) ifelse(x %in% irri_classes, 1, 0))

cell_area_ha_zohaib <- cellSize(r_irri_bin, unit = "ha")

irri_area_ha  <- r_irri_bin * cell_area_ha_zohaib
irri_area_mha <- irri_area_ha / 1e6

s_all <- c(r_class, r_irri_bin, irri_area_ha, irri_area_mha)
names(s_all) <- c("class", "irrigated_flag", "irrigated_ha", "irrigated_mha")

zohaib_dt <- as.data.table(as.data.frame(s_all, xy = TRUE, na.rm = FALSE))
setnames(zohaib_dt, c("x", "y"), c("lon", "lat"))

# Keep irrigated class 4
# (Combined irrigated area (detected by two or more indicators) ----------------

zohaib_dt <- zohaib_dt[class %in% 4]

# Countries + continents -------------------------------------------------------

zohaib_dt <- add_country_continent(zohaib_dt, countries_sf, chunk_size = 200000L)
zohaib_dt <- zohaib_dt[!is.na(country) & !is.na(continent)]

zohaib_dt[, dataset:= "zohaib"]
zohaib_dt[, mha:= irrigated_mha]

zohaib_dt <- zohaib_dt[, ..cols_to_retrieve]
zohaib_dt <- na.omit(zohaib_dt)

zohaib_dt[, sum(mha)]

fwrite(zohaib_dt, file.path(out_dir, "zohaib_dt.csv"))

# END RUNS #####################################################################
################################################################################
################################################################################


