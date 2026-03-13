## ----setup, include=FALSE-----------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, dev = "pdf", cache = TRUE)


## ----crop_regridding, message=FALSE, results = "hide"-------------------------------------------

# ARRANGE CROP DATASETS ########################################################
################################################################################

sensobol::load_packages(c("data.table", "terra", "magrittr", "here", "rnaturalearth",
                          "sf", "countrycode"))

# Source all helper functions --------------------------------------------------

r_functions <- list.files(path = here("functions"), pattern = "\\.R$",
                          full.names = TRUE)

invisible(lapply(r_functions, source))

# Countries layer (WGS84, only name + geometry) --------------------------------

countries_sf <- rnaturalearth::ne_countries(scale = "medium",
                                            returnclass = "sf") %>%
  .[, c("name", "geometry")]

# Columns to keep for gridded products -----------------------------------------

cols_to_retrieve <- c("lon", "lat", "country", "continent", "mha", "dataset")

# GAEZ 2015 ####################################################################
################################################################################

# LIST ALL IRRIGATED CROP PRODUCTION FILES -------------------------------------

# Base directory ---------------------------------------------------------------
gaez_dir <- "./crop_datasets/GAEZ_2015/GAEZ._2015_crop_production"

# List the irrigated rasters ---------------------------------------------------

crop_files <- list.files(gaez_dir,
                         pattern = "^GAEZAct2015_Production_.*_Irrigated\\.tif$",
                         full.names = TRUE)

# CREATE CLEAN CROP NAMES ------------------------------------------------------

crop_names_raw <- basename(crop_files)

crop_names <- sub("^GAEZAct2015_Production_", "", crop_names_raw)
crop_names <- sub("_Irrigated\\.tif$", "", crop_names)
crop_names <- tolower(crop_names)
crop_names <- gsub("\\.+", "_", crop_names)
crop_names <- make.names(crop_names, unique = TRUE)

crop_names
# e.g. "banana", "barley", "cassava", ..., "yamsandotherroots"

#  READ AS SPATRASTER STACK ----------------------------------------------------

r_crops <- rast(crop_files)
names(r_crops) <- crop_names

r_crops

#  Convert to wide DT ----------------------------------------------------------

gaez_dt_native <- as.data.table(as.data.frame(
  r_crops, xy = TRUE, cells = FALSE, na.rm = FALSE))

setnames(gaez_dt_native, c("x", "y"), c("lon", "lat"))

# Treat no-data as zero production ---------------------------------------------

crop_cols <- setdiff(names(gaez_dt_native), c("lon", "lat"))

for (j in crop_cols) {
  set(gaez_dt_native, which(is.na(gaez_dt_native[[j]])), j, 0)
}

# Paths to our irrigated area stacks -------------------------------------------

template_02_path <- "./mha_stack_ll_02_aligned_NEAR.tif"
template_04_path <- "./mha_stack_ll_04_aligned_NEAR.tif"
template_1_path  <- "./mha_stack_ll_10_aligned_NEAR.tif"

# CREATE FILES -----------------------------------------------------------------

# 0.2°--------------------------------------------------------------------------

gaez_crops_02_aligned <- regrid_gaez_to_template(dt_native = gaez_dt_native,
                                                 res = 0.2,
                                                 template_path = template_02_path)

gaez_crops_02_aligned <- add_country_continent(gaez_crops_02_aligned, countries_sf,
                                               chunk_size = 200000L) %>%
  .[!is.na(country)]

# 0.4° -------------------------------------------------------------------------

gaez_crops_04_aligned <- regrid_gaez_to_template(dt_native = gaez_dt_native,
                                                 res = 0.4,
                                                 template_path = template_04_path)

gaez_crops_04_aligned <- add_country_continent(gaez_crops_04_aligned, countries_sf,
                                               chunk_size = 200000L) %>%
  .[!is.na(country)]

# 1° ---------------------------------------------------------------------------

gaez_crops_1_aligned <- regrid_gaez_to_template(dt_native = gaez_dt_native,
                                                res = 1.0,
                                                template_path = template_1_path)

gaez_crops_1_aligned <- add_country_continent(gaez_crops_1_aligned, countries_sf,
                                              chunk_size = 200000L) %>%
  .[!is.na(country)]

# EXPORT FILES -----------------------------------------------------------------

fwrite(gaez_crops_02_aligned, "./datasets/crops/gaez2015_crops_irrigated_0p2deg_aligned.csv")
fwrite(gaez_crops_04_aligned, "./datasets/crops/gaez2015_crops_irrigated_0p4deg_aligned.csv")
fwrite(gaez_crops_1_aligned,  "./datasets/crops/gaez2015_crops_irrigated_1deg_aligned.csv")

################################################################################
#################### END RUNS ##################################################
################################################################################


