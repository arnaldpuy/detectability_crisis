
# CLEAN PRODUCTION FUN #########################################################

clean_production_raster <- function(r, max_valid_value = 1e9) {
  r[r < 0] <- NA
  r[r > max_valid_value] <- NA
  r
}

# HARMONIZE PRODUCTION #########################################################

harmonize_one_production_to_file <- function(file, crop, template_r, out_file,
                                             method = "bilinear",
                                             max_valid_value = 1e9,
                                             overwrite = TRUE) {
  message("Harmonizing ", crop, " -> ", basename(out_file))
  
  r <- rast(file)
  names(r) <- crop
  
  r <- clean_production_raster(r, max_valid_value = max_valid_value)
  
  orig_total <- as.numeric(global(r, "sum", na.rm = TRUE)[1, 1])
  
  prod_target <- resample(r, template_r, method = method)
  
  prod_target <- clamp(prod_target, lower = 0, values = TRUE)
  prod_target <- clean_production_raster(prod_target, max_valid_value = max_valid_value)
  
  new_total <- as.numeric(global(prod_target, "sum", na.rm = TRUE)[1, 1])
  fac <- if (is.finite(orig_total / new_total) && new_total != 0) orig_total / new_total else 1
  prod_target <- prod_target * fac
  names(prod_target) <- crop
  
  writeRaster(prod_target,
              filename = out_file,
              overwrite = overwrite,
              gdal = c("COMPRESS=LZW"))
  
  rm(r, prod_target)
  gc()
  
  data.table(crop = crop,
             file_out = out_file,
             production_tons_raw_clean = orig_total)
}

# LIST AVAILABLE CROP NAMES IN A FOLDER ########################################

get_available_crop_names <- function(folder, res_tag) {
  files <- list.files(folder,
                      pattern = paste0("_Production_ll_", res_tag, "\\.tif$"),
                      full.names = FALSE)
  
  sub(paste0("_Production_ll_", res_tag, "\\.tif$"), "", files)
}

# BUILD ONE TABLE WITH TOTAL CEREAL PRODUCTION #################################

# This reads only the selected cereal rasters, merges them by lon/lat,
# and computes a single column grain_prod = sum of selected cereals.

build_cereal_sum_table <- function(folder, res_tag,
                                   grain_targets = c("wheat", "rice", "maize", "barley", "millet", "sorghum"),
                                   grain_aliases = c("othercereals", "cerealnes")) {
  
  available <- get_available_crop_names(folder, res_tag)
  
  other_cereal_name <- intersect(grain_aliases, available)
  
  if (length(other_cereal_name) > 0) {
    
    selected_crops <- c(grain_targets, other_cereal_name[1])
    
  } else {
    
    selected_crops <- grain_targets
  }
  
  selected_crops <- intersect(selected_crops, available)
  
  if (length(selected_crops) == 0) {
    stop("No requested cereal crops found in ", folder, " for res_tag = ", res_tag)
  }
  
  message("Using crops for ", res_tag, ": ", paste(selected_crops, collapse = ", "))
  
  dt_list <- lapply(selected_crops, function(crop) {
    f <- file.path(folder, paste0(crop, "_Production_ll_", res_tag, ".tif"))
    r <- rast(f)
    
    dt <- as.data.table(as.data.frame(r, xy = TRUE, na.rm = FALSE))
    setnames(dt, c("x", "y", names(r)), c("lon", "lat", crop))
    dt
  })
  
  dt_wide <- Reduce(function(x, y) merge(x, y, by = c("lon", "lat"), all = TRUE), dt_list)
  
  dt_wide[, grain_prod := rowSums(.SD, na.rm = TRUE), .SDcols = selected_crops]
  dt_wide <- dt_wide[!is.na(grain_prod) & grain_prod > 0, .(lon, lat, grain_prod)]
  
  attr(dt_wide, "selected_crops") <- selected_crops
  dt_wide
}


# COMPUTE BREADBASKET LOSSES ###################################################

compute_breadbaskets_earthstat <- function(res_tag, grain_prod_dt, dt_pres, 
                                           regions_dt, template_r) {
  
  # Keep only the target resolution-------------------------------------------
  
  dt_pres_res <- unique( dt_pres[resolution == res_tag, .(lon, lat, n_pos)])
  
  # Map irrigation table to template cell IDs -------------------------------
  dt_pres_res[, cell := cellFromXY(template_r, as.matrix(.SD)), 
              .SDcols = c("lon", "lat")]
  dt_pres_res <- dt_pres_res[!is.na(cell)]
  
  # If multiple rows map to same cell, keep the maximum n_pos ---------------
  
  dt_pres_cell <- dt_pres_res[, .(n_pos = max(n_pos, na.rm = TRUE)), by = cell]
  
  # Map grain production table to template cell IDs--------------------------
  
  grain_dt <- copy(grain_prod_dt)
  grain_dt[, cell := cellFromXY(template_r, as.matrix(.SD)), 
           .SDcols = c("lon", "lat")]
  grain_dt <- grain_dt[!is.na(cell)]
  
  # Aggregate grain production by cell in case of duplicates----------------
  
  grain_cell <- grain_dt[, .(grain_prod = sum(grain_prod, na.rm = TRUE)), by = cell]
  
  # Recover exact template cell centres for regional clipping----------------
  
  xy <- as.data.table(xyFromCell(template_r, grain_cell$cell))
  setnames(xy, c("x", "y"), c("lon", "lat"))
  grain_cell <- cbind(grain_cell, xy)
  
  # Merge by cell index, not by lon/lat
  
  prod_join <- merge(grain_cell, dt_pres_cell, by = "cell", all.x = TRUE)
  prod_join[is.na(n_pos), n_pos:= 0L]
  
  ks <- 10:1
  
  out <- rbindlist(
    lapply(seq_len(nrow(regions_dt)), function(i) {
      
      rr <- regions_dt[i]
      
      dt_region <- prod_join[
        lon >= rr$xmin & lon <= rr$xmax &
          lat >= rr$ymin & lat <= rr$ymax
      ]
      
      totals <- vapply(
        ks,
        function(k) dt_region[n_pos >= k, sum(grain_prod, na.rm = TRUE)],
        numeric(1)
      )
      names(totals) <- as.character(ks)
      
      base_1 <- totals["1"]
      base_5 <- totals["5"]
      
      if (is.na(base_1) || base_1 <= 0) {
        frac_vs_1 <- rep(NA_real_, length(totals))
        loss_vs_1 <- rep(NA_real_, length(totals))
      } else {
        frac_vs_1 <- totals / base_1
        loss_vs_1 <- 1 - totals / base_1
      }
      
      if (is.na(base_5) || base_5 <= 0) {
        frac_vs_5 <- rep(NA_real_, length(totals))
        loss_vs_5 <- rep(NA_real_, length(totals))
      } else {
        frac_vs_5 <- totals / base_5
        loss_vs_5 <- 1 - totals / base_5
      }
      
      data.table(
        region = rr$region,
        resolution = res_tag,
        k = ks,
        total_prod = as.numeric(totals),
        frac_vs_1 = as.numeric(frac_vs_1),
        loss_vs_1 = as.numeric(loss_vs_1),
        frac_vs_5 = as.numeric(frac_vs_5),
        loss_vs_5 = as.numeric(loss_vs_5)
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )
  
  # loss_vs_5 only makes sense for k >= 5
  out[k < 5, `:=`(
    frac_vs_5 = NA_real_,
    loss_vs_5 = NA_real_
  )]
  
  out
}
