## ----setup, include=FALSE-----------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, dev = "pdf", cache = TRUE)


## ----warning=FALSE, message=FALSE, results = "hide"---------------------------------------------

# Load libraries ---------------------------------------------------------------

sensobol::load_packages(c("data.table", "scales", "cowplot", "tidyverse", "sensobol",
                          "here", "terra", "rnaturalearth", "sf", "countrycode",
                          "wesanderson", "benchmarkme", "parallel"))

# Create custom theme ----------------------------------------------------------

theme_AP <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA),
          legend.key.width = unit(0.4, "cm"),
          legend.key.height = unit(0.5, "lines"),
          legend.key.spacing.y = unit(0, "lines"),
          legend.box.spacing = unit(0, "pt"),
          legend.spacing.y  = unit(0.1, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size = 7.3),
          axis.title.y = element_text(size = 7.3),
          axis.title = element_text(size = 10),
          plot.title = element_text(size = 8),
          strip.text.x = element_text(size = 7.4),
          strip.text.y = element_text(size = 7.4),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(margin = margin(t = 1.5, b = 1.5)))
}


# Source all .R files in the "functions" folder -------------------------------

r_functions <- list.files(path = here("functions"),
                          pattern = "\\.R$", full.names = TRUE)

lapply(r_functions, source)

# Set seed ---------------------------------------------------------------------

seed <- 123

# Load countries as sf (medium scale is a good balance) -----------------------

countries_sf <- ne_countries(scale = "medium", returnclass = "sf")

# Keep only what is needed -----------------------------------------------------

countries_sf <- countries_sf[, c("name", "geometry")]

# Define colors and palettes and plotting datasets -----------------------------

tau_palette <- c("Agreement"= "#2ECC71",
                 "Low disagreement" = "#F7DC6F",
                 "Moderate disagreement" = "#F39C12",
                 "High disagreement" = "#E74C3C",
                 "Extreme disagreement" = "red")

res_palette <- c("0.2deg" = "lightblue",
                 "0.4deg" = "#117A65",
                 "1deg"   = "darkblue")


## ----uasa---------------------------------------------------------------------------------------

# LOAD ALL DATASETS ############################################################

# List all resolution files ----------------------------------------------------

files <- list.files("./datasets/irrigated_areas_regridded",
                    pattern = "irrigated_areas_regridded_.*\\.csv",
                    full.names = TRUE)

# Load them, bind them and label them ------------------------------------------

dt <- rbindlist(lapply(files, function(f) {
  x <- fread(f)
  res <- sub(".*_(\\d+)\\.csv", "\\1", f)
  x[, resolution:= paste0(as.numeric(res) / 10, "deg")]
  x
}))

# Load cropland mask -----------------------------------------------------------

crop_path  <- "./irrigated_area_datasets/GirrEO/asap_mask_crop_v02.tif"
crop_native <- rast(crop_path)

# Retrieve only cells within cropland ------------------------------------------

dt <- filter_long_dt_by_crop_mask_aggregated(dt = dt, crop_native = crop_native,
                                             rule = "any")

# DEFINE TAU GRID ##############################################################

# Vector with the datasets used ------------------------------------------------

dataset_names <- unique(dt$dataset)
new_dataset_names <- c("GAEZ+2015", "GIAM", "GMIA", "GRIPC", "LUH2", "Meier",
                       "MIRCA 2000", "MIRCA OS", "Nagaraj", "SPAM2010")

# Define a dense grid of tau values from 1 ha (10^0) to 1e5 ha (10^5) ----------

log_tau_min <- 0
log_tau_max <- 5
log_step <- 0.1

taus_ha_nonzero <- 10^(seq(log_tau_min, log_tau_max, by = log_step))
taus_ha <- c(0, taus_ha_nonzero)
taus_mha <- taus_ha * 1e-6

# Labels: "any_positive", then "1ha", "1.26ha", ..., "1e+05ha"---

tau_labels <- ifelse(taus_ha == 0, "any_positive", sprintf("%gha", signif(taus_ha, 3)))

# Define tau grid ------------------------------------

tau_grid <- data.table(tau_label = tau_labels, tau_ha = taus_ha, tau_mha = taus_mha)

# TRANSFORM THE SAMPLE MATRIX ##################################################

# Define the settings ---------------------------------------------------------

N <- 2^13
params <- c("resolution", "tau", "exclusion")
matrices <- c("A", "B", "AB")
R <- 10^3
boot <- TRUE

# Resolution vector ------------------------------------------------------------

resolution_vec <- c("0.2deg", "0.4deg", "1deg")
n_resolution_vec <- length(resolution_vec)

# Tau vector -------------------------------------------------------------------

tau_vec <- tau_grid$tau_ha
n_tau_vec <- length(tau_vec)

# Dataset exclusion vector -----------------------------------------------------

dataset_vec <- c("all", dataset_names)
n_dataset_vec <- length(dataset_vec)

# Create the sample matrix -----------------------------------------------------

mat <- data.table(sobol_matrices(matrices = matrices, N = N, params = params))

# Transform resolution column --------------------------------------------------

mat[, resolution:= pmin(n_resolution_vec, 1L + floor(resolution * n_resolution_vec))]
mat[, resolution:= resolution_vec[resolution]]

# Transform tau column ---------------------------------------------------------

mat[, tau:= pmin(n_tau_vec, 1L + floor(tau * n_tau_vec))]
mat[, tau:= tau_vec[tau]]

# Transform exclusion column ---------------------------------------------------

mat[, exclusion:= pmin(n_dataset_vec, 1L + floor(exclusion * n_dataset_vec))]
mat[, exclusion:= dataset_vec[exclusion]]

# PRE-COMPUTE ALL COMBINATIONS #################################################

# Compute all possibilities by looping over resolutions and tau values ---------

# Set parallel backend --------------------------

resolutions <- sort(unique(dt$resolution))
ncores <- max(1, parallel::detectCores() - 1)

# loop over resolutions and scenarios in parallel

tau_results <- rbindlist(
  parallel::mclapply(resolutions, function(res) {

    dt_res_long <- dt[resolution == res]

    # datasets actually present at this resolution
    ds_res <- sort(unique(dt_res_long$dataset))

    rbindlist(
      lapply(dataset_vec, function(scn) {

        # which datasets to use -----------------------------------------------

        if (scn == "all") {

          use_ds <- ds_res

        } else {

          use_ds <- setdiff(ds_res, scn)
        }

        # if only one dataset left, skip ---------------------------------------
        if (length(use_ds) < 2L) {

          return(NULL)
        }

        # reshape to wide ------------------------------------------------------

        dt_res_wide <- dcast(dt_res_long[dataset %in% use_ds],
          lon + lat + country + code + continent ~ dataset,
          value.var = "mha", fill = 0
        )

        dt_res_wide[, `:=`(A_min = do.call(pmin, c(.SD, list(na.rm = TRUE))),
                  A_max = do.call(pmax, c(.SD, list(na.rm = TRUE)))),
           .SDcols = use_ds]

        # compute tau ----------------------------------------------------------

        out <- tau_grid[, compute_tau_fun(dt_res_wide, tau_mha, use_ds),
                        .(tau_label, tau_ha, tau_mha)]


        # tag with resolution + scenario ---------------------------------------

        out[, `:=`(resolution = res, scenario = scn)]

        out
      }),
      use.names = TRUE, fill = TRUE
    )
  }, mc.cores = ncores),
  use.names = TRUE, fill = TRUE)

# ARRANGE OUTPUT ###############################################################

# Rename tau_ha to tau to match mat --------------------------------------------

setnames(tau_results, "tau_ha", "tau")

# Join all three disagreement metrics onto mat ---------------------------------

mat[tau_results, on = .(resolution, tau, exclusion = scenario),
    `:=`(frac_disagree = i.frac_disagree,
         frac_major_disagree = i.frac_major_disagree,
         frac_minor_disagree = i.frac_minor_disagree,
         share_existence = i.share_existence,
         share_marginal = i.share_marginal)]

# Reshape just these three for plotting ----------------------------------------

mat_long <- melt(mat, id.vars = c("resolution", "tau", "exclusion"),
                 measure.vars = c("frac_disagree", "frac_major_disagree",
                                  "frac_minor_disagree", "share_existence",
                                  "share_marginal"),
                 variable.name = "agreement", value.name = "frac_value")


## ----plot_uasa, dependson="uasa", fig.height=3, fig.width=3-------------------------------------

# PREPARE DATA TO PLOT UNCERTAINTY #############################################

# Calculate mean and quantiles by tau × agreement ------------------------------

mat_sum <- mat_long[, .(frac_mean = mean(frac_value, na.rm = TRUE),
                        frac_median = median(frac_value, na.rm = TRUE),
                        frac_lwr = quantile(frac_value, 0.025, na.rm = TRUE),
                        frac_upr = quantile(frac_value, 0.975, na.rm = TRUE)),
                    .(tau, agreement)]

# Get values at tau = 0, 1, 10, 500, 1000, 5000, 10000 -------------------------

target_tau <- c(0, 1, 10, 500, 1000, 1600, 2000, 5000, 10000, 13000)

mat_sum[, {dt_sub <- .SD
idx <- sapply(target_tau, function(x) {
  which.min(abs(log10(tau + 1e-12) - log10(x + 1e-12)))
      }
    )
    dt_sub[idx]
  },
  agreement
] %>%
  .[, .(agreement, tau, frac_mean, frac_lwr, frac_upr)]

# Plot uncertainty -------------------------------------------------------------

agreement_cols <- c(frac_disagree = "black",
                    frac_major_disagree = "#B2182B",
                    frac_minor_disagree = "#2166AC"  )

plot_uncertainty1 <- mat_sum[!agreement %in% c("share_existence", "share_marginal")] %>%
  ggplot(.,
  aes(x = tau, y = frac_median, colour = agreement, fill = agreement)) +
  geom_ribbon(aes(ymin = frac_lwr, ymax = frac_upr), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_x_log10(labels = label_log(digits = 2)) +
  scale_colour_manual(values = agreement_cols,
                      breaks = c("frac_disagree",
                                 "frac_major_disagree",
                                 "frac_minor_disagree"),
                      labels = c("Total",
                                 "Major",
                                 "Minor"),
                      name = "Disagreement") +
  scale_fill_manual(values = agreement_cols,
                    breaks = c("frac_disagree",
                               "frac_major_disagree",
                               "frac_minor_disagree"),
                    labels = c("Total",
                               "Major",
                               "Minor"),
                    name = "Disagreement") +
  labs(x = expression(tau~"(ha)"), y = "Fraction of cells") +
  theme_AP() +
  theme(legend.position = c(0.65, 0.77))

plot_uncertainty1

# Set colors for existential and marginal disagreement -------------------------

exist_cols <- c("share_existence" = "#1f78b4", "share_marginal" = "#e66101")

plot_agreement1 <- mat_sum[agreement %in% c("share_existence", "share_marginal")] %>%
  ggplot(., aes(x = tau, y = frac_median, colour = agreement, fill = agreement)) +
  geom_ribbon(aes(ymin = frac_lwr, ymax = frac_upr), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = exist_cols,
                      labels = c(share_existence = "Existential",
                                 share_marginal = "Marginal"),
                      name = "Disagreement") +
  scale_fill_manual(values = exist_cols,
                    labels = c(share_existence = "Existential",
                               share_marginal = "Marginal"),
                    name = "Disagreement") +
  scale_x_log10(labels = label_log(digits = 2)) +
  labs(x = expression(tau~"(ha)"), y = "Fraction of cells") +
  theme_AP() +
  theme(legend.position = c(0.3, 0.7))


plot_agreement1

# Around tau > 10^4 ha, the existential and marginal components intersect. That
# intersection has a clear conceptual meaning:
# -It is the tolerance level at which disagreement about whether irrigation exists
# at all becomes less dominant than disagreement about how much irrigation exists.
# -Still small relative to grid-cell area. 1000 ha is only 2–3% of a 0.2º grid cell.
# That makes it a conservative tolerance as we are not allowing massive within-cell
# divergence.

# SENSITIVITY ANALYSIS #########################################################

# Sobol' indices ---------------------------------------------------------------

ind <- sobol_indices(matrices = matrices, N = N, params = params, Y
                     = mat[, frac_disagree], R = R, boot = TRUE)

print(ind)

plot_indices1 <- plot(ind) +
  labs(x = "", y = "Fraction of variance") +
  theme_AP() +
  scale_x_discrete(labels = c(exclusion = "exclusion",
                              resolution = "resolution",
                              tau = expression(tau)),
                   guide = guide_axis(n.dodge = 2)) +
  theme(legend.position = c(0.3, 0.8))

plot_indices1

# SCATTERPLOTS OF OUTPUT AGAINS INTPUT #########################################

label_map <- setNames(new_dataset_names, dataset_names)

a1 <- ggplot(mat[!exclusion == "all", c("exclusion", "frac_disagree")],
            aes(exclusion, frac_disagree)) +
  geom_boxplot() +
  stat_summary_bin(fun = "mean", geom = "point",
                   colour = "red", size = 1) +
  scale_x_discrete(labels = label_map) +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  theme_AP() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "exclusion", y = "Fraction disagree")

a1

b1 <- ggplot(mat[1:N, c("resolution", "frac_disagree")],
       aes(resolution, frac_disagree)) +
  geom_boxplot() +
  scale_x_discrete(labels = c(`0.2deg` = "0.2º",
                              `0.4deg` = "0.4º",
                              `1deg` = "1.0º")) +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  stat_summary_bin(fun = "mean", geom = "point",
                   colour = "red", size = 1) +
  theme_AP() +
  labs(x = "resolution", y = "")

b1

c1 <- ggplot(mat[1:N, c("tau", "frac_disagree")],
       aes(tau, frac_disagree)) +
  geom_point(size = 0.1, alpha = 0.1) +
  stat_summary_bin(fun = "mean", geom = "point",
                   colour = "red", size = 1) +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  theme_AP() +
  scale_x_log10(labels = label_log(digits = 2)) +
  labs(x = expression(tau~"(ha)"), y = "")

c1


## ----all_maps, dependson=c("uasa", "plot_uasa"), fig.height=1.5, fig.width=5.8------------------

# PLOT DATASETS ################################################################

tmp <- dt[resolution == "0.2deg"] %>%
  .[, .(mha = sum(mha)), dataset] %>%
  .[, continent:= "GLOBAL"]

dt[resolution == "0.2deg"] %>%
  .[, .(mha = sum(mha)), .(continent, dataset)] %>%
  rbind(., tmp) %>%
  .[, continent:= factor(continent, levels = c("Africa", "Americas",
                                               "Asia", "Europe", "Oceania",
                                               "GLOBAL"))] %>%
  ggplot(., aes(dataset, mha, fill = dataset)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired", label = label_map, name = "") +
  facet_wrap(~continent, ncol = 6, scales = "free_y") +
  theme_AP() +
  labs(x = "", y = "Mha") +
  theme(legend.position = "top") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())



## ----merge_uasa, dependson="plot_uasa", fig.height=3.5, fig.width=5.5---------------------------

# MERGE UA/SA PLOTS ############################################################

top <- plot_grid(plot_uncertainty1, plot_agreement1, plot_indices1,
                 ncol = 3, labels = "auto", rel_widths = c(0.4, 0.3, 0.3))
bottom <- plot_grid(a1, b1, c1, ncol = 3, rel_widths = c(0.5, 0.25, 0.25), labels = "d")
plot_grid(top, bottom, ncol = 1)


## ----effects_coarsening, dependson="uasa", fig.height=1.8, fig.width=1.8------------------------

# DOES ZOOMING OUT REDUCE DISAGREEMENT? ########################################

effect_resolution_dt <- tau_results[scenario == "all", .(tau, frac_disagree, resolution)]

# Order resolution -------------------------------------------------------------

effect_resolution_dt[, resolution:= factor(resolution,
                                           levels = c("0.2deg", "0.4deg", "1deg"))]

# To wide table ----------------------------------------------------------------

wide_dt <- dcast(effect_resolution_dt, tau ~ resolution, value.var = "frac_disagree")

# Rename -----------------------------------------------------------------------

setnames(wide_dt, old = c("0.2deg", "0.4deg", "1deg"), new = c("res_0.2", "res_0.4", "res_1"))

# Compute absolute change (fine - coarse) --------------------------------------

wide_dt[, abs_change:= res_0.2 - res_1]

# Compute elasticity (log-log) -------------------------------------------------

# Area ratio from 0.2 to 1
area_ratio <- 25
wide_dt[, elasticity:= ifelse(res_0.2 > 0 & res_1 > 0,
                              log(res_1 / res_0.2) / log(area_ratio), NA_real_)]

# Plot disagreement curves -----------------------------------------------------

p1 <- ggplot(effect_resolution_dt, aes(tau, frac_disagree, colour = resolution)) +
  geom_line(linewidth = 0.8) +
  scale_x_log10(labels = label_log(digits = 2)) +
  scale_color_manual(values = res_palette) +
  labs(x = expression(tau~"(ha)"),
       y = "Fraction of cells",
       colour = "Resolution") +
  theme_AP() +
  theme(legend.position = c(0.3, 0.35))

p1

# Plot absolute change ---------------------------------------------------------

p2 <- ggplot(wide_dt, aes(tau, abs_change)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  scale_x_log10(labels = label_log(digits = 2)) +
  labs(x = expression(tau~"(ha)"),
       y = "Change disagreement") +
  theme_AP()

p2

# Plot elasticity --------------------------------------------------------------

p3 <- ggplot(wide_dt[!is.na(elasticity)],
       aes(x = tau, y = elasticity)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  scale_x_log10(labels = label_log(digits = 2)) +
  labs(x = expression(tau~"(ha)"),
       y = "Elasticity disagreement") +
  theme_AP()

p3



## ----plot_effects_resolution, dependson="effects_coarsening", fig.height=1.8, fig.width=5.3-----

# Merge ------------------------------------------------------------------------

plot_grid(p1, p2, p3, ncol = 3, labels = "auto")

# Increasing spatial aggregation from 0.2 to 1 (25x coarsening) reduces
# existential disagreement by c20 percentage points at near-zero
# thresholds. this effect vanishes around tau = 300–500 ha, beyond which
# aggregation increases disagreement. At higher thresholds disagreement is
# greater at 1 than at 0.2. Coarse-scale aggregation can
# amplify rather than resolve structural spatial inconsistencies between datasets.

# These results refut the criticism of "if you zoom out, the maps converge".

# ALSO: Increasing cell area by a factor of 25 (0.2 to 1) yields an elasticity
# of approximately -0.07 at tau =  0, implying that a 1% increase in spatial
# aggregation reduces disagreement by only 0.07%. Elasticity approaches zero
# at moderate thresholds and becomes positive at higher thresholds.



## ----uasa_on_mapping_paradigm, dependson = "uasa", echo=FALSE, fig.height=3, fig.width=3--------

# MAPPING PARADIGM EXCLUSION UA/SA #############################################

ds_groups <- list(remote_sensing  = c("giam", "gripc"),
                  statistical_disagg  = c("gmia", "meier", "mirca2000", "mirca_os"),
                  landuse_allocation = c("spam", "gaez_v4"),
                  machine_learning = c("nagaraj"),
                  IAM_based = c("luh2"))

mapping_paradigm_vec <- names(ds_groups)
n_mapping_paradigm_vec <- length(mapping_paradigm_vec)

# Create the sample matrix -----------------------------------------------------

mat <- data.table(sobol_matrices(matrices = matrices, N = N, params = params))

# Transform resolution column --------------------------------------------------

mat[, resolution:= pmin(n_resolution_vec, 1L + floor(resolution * n_resolution_vec))]
mat[, resolution:= resolution_vec[resolution]]

# Transform tau column ---------------------------------------------------------

mat[, tau:= pmin(n_tau_vec, 1L + floor(tau * n_tau_vec))]
mat[, tau:= tau_vec[tau]]

# Transform exclusion column ---------------------------------------------------

mat[, exclusion:= pmin(n_mapping_paradigm_vec, 1L + floor(exclusion * n_mapping_paradigm_vec))]
mat[, exclusion:= mapping_paradigm_vec[exclusion]]

# PRE-COMPUTE ALL COMBINATIONS #################################################

# Compute all possibilities by looping over resolutions and tau values ---------

# Set parallel backend --------------------------

resolutions <- sort(unique(dt$resolution))
ncores <- max(1, parallel::detectCores() - 1)

# ds_groups: named list of dataset keys to exclude by mapping category
# mat: data.table with columns resolution, tau (in ha), exclusion (category name)

idcols <- c("lon","lat","country","code","continent")

tau_results <- rbindlist(
  parallel::mclapply(resolutions, function(res) {

    dt_res_long <- dt[resolution == res]
    if (nrow(dt_res_long) == 0L) return(NULL)

    # scenarios to run for this resolution
    mat_res <- mat[resolution == res]
    if (nrow(mat_res) == 0L) return(NULL)

    # datasets present at this resolution
    ds_res <- sort(unique(dt_res_long$dataset))

    # ---- CAST ONCE per resolution  ----

    dt_res_wide_all <- dcast(dt_res_long,
                             lon + lat + country + code + continent ~ dataset,
                             value.var = "mha", fill = 0)

    # loop exclusion groups (few) instead of mat rows (many)
    excl_levels <- unique(mat_res$exclusion)

    rbindlist(lapply(excl_levels, function(excl_grp) {

      # datasets to exclude for this group
      excl_ds <- ds_groups[[excl_grp]]
      if (is.null(excl_ds)) excl_ds <- character(0)

      use_ds <- setdiff(ds_res, excl_ds)
      if (length(use_ds) < 2L) return(NULL)

      # taus used for THIS (resolution, exclusion-group)
      taus_ha <- unique(mat_res[exclusion == excl_grp, tau])
      taus_mha <- taus_ha * 1e-6

      # ---- compute A_min/A_max ONCE per exclusion group ----
      dt_w <- dt_res_wide_all[, c(idcols, use_ds), with = FALSE]

      dt_w[, `:=`(A_min = do.call(pmin, c(.SD, list(na.rm = TRUE))),
                  A_max = do.call(pmax, c(.SD, list(na.rm = TRUE)))),
           .SDcols = use_ds]

      # ---- now only loop over taus
      out_grp <- rbindlist(lapply(seq_along(taus_mha), function(k) {

        tau_ha  <- taus_ha[k]
        tau_mha <- taus_mha[k]

        out <- compute_tau_fun(dt_w, tau_mha, use_ds)

        out[, `:=`(
          tau_label = paste0(format(tau_ha, scientific = FALSE, trim = TRUE), "ha"),
          tau_ha = tau_ha,
          tau_mha = tau_mha,
          resolution = res,
          exclusion  = excl_grp
        )]

        out
      }), use.names = TRUE, fill = TRUE)

      out_grp
    }), use.names = TRUE, fill = TRUE)

  }, mc.cores = ncores),
  use.names = TRUE, fill = TRUE
)

# ARRANGE OUTPUT ##############################################################

# Rename tau_ha to tau to match mat --------------------------------------------

setnames(tau_results, "tau_ha", "tau")

# Join all three disagreement metrics onto mat ---------------------------------

mat[tau_results, on = .(resolution, tau, exclusion = exclusion),
    `:=`(frac_disagree = i.frac_disagree,
         frac_major_disagree = i.frac_major_disagree,
         frac_minor_disagree = i.frac_minor_disagree,
         share_existence = i.share_existence,
         share_marginal = i.share_marginal)]

# Reshape just these three for plotting ----------------------------------------

mat_long <- melt(mat, id.vars = c("resolution", "tau", "exclusion"),
                 measure.vars = c("frac_disagree", "frac_major_disagree",
                                  "frac_minor_disagree", "share_existence",
                                  "share_marginal"),
                 variable.name = "agreement", value.name = "frac_value")

# PREPARE DATA TO PLOT UNCERTAINTY #############################################

# Calculate mean and quantiles by tau × agreement ------------------------------

mat_sum <- mat_long[, .(frac_mean = mean(frac_value, na.rm = TRUE),
                        frac_median = median(frac_value, na.rm = TRUE),
                        frac_lwr = quantile(frac_value, 0.025, na.rm = TRUE),
                        frac_upr = quantile(frac_value, 0.975, na.rm = TRUE)),
                    .(tau, agreement)]

mat_sum[, {dt_sub <- .SD
idx <- sapply(target_tau, function(x) {
  which.min(abs(log10(tau + 1e-12) - log10(x + 1e-12)))
      }
    )
    dt_sub[idx]
  },
  agreement
]

# Plot uncertainty -------------------------------------------------------------

agreement_cols <- c(frac_disagree = "black",
                    frac_major_disagree = "#B2182B",
                    frac_minor_disagree = "#2166AC"  )

plot_uncertainty <- mat_sum[!agreement %in% c("share_existence", "share_marginal")] %>%
  ggplot(.,
  aes(x = tau, y = frac_median, colour = agreement, fill = agreement)) +
  geom_ribbon(aes(ymin = frac_lwr, ymax = frac_upr), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_x_log10(labels = label_log(digits = 2)) +
  scale_colour_manual(values = agreement_cols,
                      breaks = c("frac_disagree",
                                 "frac_major_disagree",
                                 "frac_minor_disagree"),
                      labels = c("Total",
                                 "Major (\u2265 2 vs \u2265 2)",
                                 "Minor (one outlier)"),
                      name = "Disagreement") +
  scale_fill_manual(values = agreement_cols,
                    breaks = c("frac_disagree",
                               "frac_major_disagree",
                               "frac_minor_disagree"),
                    labels = c("Total",
                               "Major (\u2265 2 vs \u2265 2)",
                               "Minor (one outlier)"),
                    name = "Disagreement") +
  labs(x = expression(tau~"(ha)"), y = "Fraction of cells") +
  theme_AP() +
  theme(legend.position = c(0.65, 0.77))

plot_uncertainty

plot_agreement <- mat_sum[agreement %in% c("share_existence", "share_marginal")] %>%
  ggplot(., aes(x = tau, y = frac_median, colour = agreement, fill = agreement)) +
  geom_ribbon(aes(ymin = frac_lwr, ymax = frac_upr), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = exist_cols,
                      labels = c(share_existence = "Existential",
                                 share_marginal = "Marginal"),
                      name = "Disagreement") +
  scale_fill_manual(values = exist_cols,
                    labels = c(share_existence = "Existential",
                               share_marginal = "Marginal"),
                    name = "Disagreement") +
  scale_x_log10(labels = label_log(digits = 2)) +
  labs(x = expression(tau~"(ha)"), y = "Fraction of cells") +
  theme_AP() +
  theme(legend.position = c(0.6, 0.8))


plot_agreement

# SENSITIVITY ANALYSIS #########################################################

# Sobol' indices ---------------------------------------------------------------

ind <- sobol_indices(matrices = matrices, N = N, params = params, Y
                     = mat[, frac_disagree], R = R, boot = TRUE)

plot_indices2 <- plot(ind) +
  theme_AP() +
  labs(x = "", y = "Fraction of variance") +
  scale_x_discrete(labels = c(exclusion = "exclusion",
                              resolution = "resolution",
                              tau = expression(tau)),
                   guide = guide_axis(n.dodge = 2)) +
  theme(legend.position = c(0.3, 0.8))

plot_indices2

# SCATTERPLOTS OF OUTPUT AGAINS INTPUT #########################################

label_map <- setNames(new_dataset_names, dataset_names)

a2<- ggplot(mat[!exclusion == "all", c("exclusion", "frac_disagree")],
            aes(exclusion, frac_disagree)) +
  geom_boxplot() +
  stat_summary_bin(fun = "mean", geom = "point",
                   colour = "red", size = 1) +
  scale_x_discrete(labels = label_map) +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  theme_AP() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "exclusion", y = "Fraction disagree")

a2

b2 <- ggplot(mat[1:N, c("resolution", "frac_disagree")],
            aes(resolution, frac_disagree)) +
  geom_boxplot() +
  scale_x_discrete(labels = c(`0.2deg` = "0.2º",
                              `0.4deg` = "0.4º",
                              `1deg` = "1.0º")) +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  stat_summary_bin(fun = "mean", geom = "point",
                   colour = "red", size = 1) +
  theme_AP() +
  labs(x = "resolution", y = "")

b2

c2 <- ggplot(mat[1:N, c("tau", "frac_disagree")],
            aes(tau, frac_disagree)) +
  geom_point(size = 0.1, alpha = 0.1) +
  stat_summary_bin(fun = "mean", geom = "point",
                   colour = "red", size = 1) +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  theme_AP() +
  scale_x_log10(labels = label_log(digits = 2)) +
  labs(x = expression(tau~"(ha)"), y = "")

c2


## ----merge_dataset_classes, fig.height=3.8, fig.width=5.5---------------------------------------

# MERGE UA/SA PLOTS ############################################################

top <- plot_grid(plot_uncertainty, plot_agreement, plot_indices2,
                 ncol = 3, labels = "auto", rel_widths = c(0.4, 0.3, 0.3))
bottom <- plot_grid(a2, b2, c2, ncol = 3, rel_widths = c(0.4, 0.3, 0.3), labels = "d")
plot_grid(top, bottom, ncol = 1)


## ----uasa_on_weights, dependson = "uasa", echo=FALSE, fig.height=3, fig.width=3-----------------

# UA/SA ON THE WEIGHTING SCHEME ################################################

year_vector <- c(2015, 2009, 2013, 2015, 2020, 2018,
                 2010, 2025, 2021, 2020)

ds_meta <- data.table(dataset = dataset_names, year =  year_vector)

# Building the weighting schemes -----------------------------------------------

# Linear weights with year (rescaled to avoid zeros)
ds_meta[, w_linear_raw:= year - min(year) + 1]
ds_meta[, w_linear:= w_linear_raw / sum(w_linear_raw)]

# Exponential bias toward newer datasets
lambda <- 0.1  # tuning parameter; larger to stronger bias toward recent
ds_meta[, w_exp_raw:= exp(lambda * (year - max(year)))]
ds_meta[, w_exp:= w_exp_raw / sum(w_exp_raw)]

# Put into list ----------------------------------------------------------------

weight_scenarios <- list(linear = ds_meta[, setNames(w_linear, dataset)],
                         exp = ds_meta[, setNames(w_exp, dataset)])

n_weight <- length(weight_scenarios)

# Create the sample matrix -----------------------------------------------------

params <- c("resolution", "tau", "exclusion", "weight")
mat <- data.table(sobol_matrices(matrices = matrices, N = N, params = params))

# Transform resolution column --------------------------------------------------

mat[, resolution:= pmin(n_resolution_vec, 1L + floor(resolution * n_resolution_vec))]
mat[, resolution:= resolution_vec[resolution]]

# Transform tau column ---------------------------------------------------------

mat[, tau:= pmin(n_tau_vec, 1L + floor(tau * n_tau_vec))]
mat[, tau:= tau_vec[tau]]

# Transform exclusion column ---------------------------------------------------

mat[, exclusion:= pmin(n_dataset_vec, 1L + floor(exclusion * n_dataset_vec))]
mat[, exclusion:= dataset_vec[exclusion]]

#Transform weights column ------------------------------------------------------

mat[, weight:= pmin(n_weight, 1L + floor(weight * n_weight))]
mat[, weight:= names(weight_scenarios)[weight]]

# RUN LOOP TO CHECK INFLUENCE OF WEIGHTS #######################################

idcols <- c("lon", "lat", "country", "code", "continent")
tau_results <- rbindlist(
  parallel::mclapply(resolutions, function(res) {

    # data for this resolution -------------------------------------------------
    dt_res_long <- dt[resolution == res]
    if (nrow(dt_res_long) == 0L) return(NULL)

    # rows of the design matrix for this resolution ---------------------------
    mat_res <- mat[resolution == res]
    if (nrow(mat_res) == 0L) return(NULL)

    # datasets actually present at this resolution -----------------------------
    ds_res <- sort(unique(dt_res_long$dataset))

    rbindlist(
      lapply(dataset_vec, function(scn) {

        # which datasets to use -----------------------------------------------

        if (scn == "all") {
          use_ds <- ds_res
        } else {
          # leave-one-out: exclude the scenario dataset
          use_ds <- setdiff(ds_res, scn)
        }

        # if only one dataset left, skip --------------------------------------
        if (length(use_ds) < 2L) return(NULL)

        # design rows for THIS (resolution, scenario) -------------------------
        mat_scn <- unique(mat_res[exclusion == scn, .(tau, weight)])
        if (nrow(mat_scn) == 0L) return(NULL)

        # reshape to wide ONCE per (res, scenario) ----------------------------
        dt_res_wide <- dcast(
          dt_res_long[dataset %in% use_ds],
          lon + lat + country + code + continent ~ dataset,
          value.var = "mha",
          fill = 0
        )

        # A_min / A_max ONCE per (res, scenario) ------------------------------
        dt_res_wide[, `:=`(
          A_min = do.call(pmin, c(.SD, list(na.rm = TRUE))),
          A_max = do.call(pmax, c(.SD, list(na.rm = TRUE)))
        ), .SDcols = use_ds]

        # now loop over tau–weight combos -------------------------------------
        out_scn <- rbindlist(
          lapply(seq_len(nrow(mat_scn)), function(k) {

            tau_ha <- mat_scn$tau[k]
            tau_mha <- tau_ha * 1e-6
            weight_lab <- mat_scn$weight[k]

            out <- compute_tau_weighted_fun(
              dt = dt_res_wide,
              tau_mha = tau_mha,
              dataset_names = use_ds,
              weight_label = weight_lab,
              weight_scenarios = weight_scenarios )

            out[, `:=`(
              tau_label = paste0(format(tau_ha, scientific = FALSE, trim = TRUE), "ha"),
              tau_ha = tau_ha,
              tau_mha = tau_mha,
              resolution = res,
              scenario = scn,
              weight = weight_lab
            )]

            out
          }),
          use.names = TRUE, fill = TRUE
        )

        out_scn
      }),
      use.names = TRUE, fill = TRUE
    )
  }, mc.cores = 1),
  use.names = TRUE, fill = TRUE
)

# ARRANGE OUTPUT ###############################################################

# Rename tau_ha to tau to match mat --------------------------------------------

setnames(tau_results, "tau_ha", "tau")

# Join all three disagreement metrics onto mat ---------------------------------

mat[tau_results, on = .(resolution, tau, exclusion = scenario),
    `:=`(frac_disagree = i.frac_disagree,
         frac_major_disagree = i.frac_major_disagree,
         frac_minor_disagree = i.frac_minor_disagree,
         share_existence = i.share_existence,
         share_marginal = i.share_marginal)]

# Reshape just these three for plotting ----------------------------------------

mat_long <- melt(mat, id.vars = c("resolution", "tau", "exclusion"),
                 measure.vars = c("frac_disagree", "frac_major_disagree",
                                  "frac_minor_disagree", "share_existence",
                                  "share_marginal"),
                 variable.name = "agreement", value.name = "frac_value")


## ----plot_uasa, dependson="uasa"------------------------------------------------------------------

# PREPARE DATA TO PLOT UNCERTAINTY #############################################

# Calculate mean and quantiles by tau × agreement ------------------------------

mat_sum <- mat_long[, .(frac_mean = mean(frac_value, na.rm = TRUE),
                        frac_median = median(frac_value, na.rm = TRUE),
                        frac_lwr = quantile(frac_value, 0.025, na.rm = TRUE),
                        frac_upr = quantile(frac_value, 0.975, na.rm = TRUE)),
                    .(tau, agreement)]

# Get values at tau = 0, 1, 10, 100, 1000 --------------------------------------

mat_sum[tau %in% c(0, 1, 10, 100, 1000)]

# Plot uncertainty -------------------------------------------------------------

agreement_cols <- c(frac_disagree = "black",
                    frac_major_disagree = "#B2182B",
                    frac_minor_disagree = "#2166AC"  )

plot_uncertainty <- mat_sum[!agreement %in% c("share_existence", "share_marginal")] %>%
  ggplot(.,
         aes(x = tau, y = frac_median, colour = agreement, fill = agreement)) +
  geom_ribbon(aes(ymin = frac_lwr, ymax = frac_upr), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_x_log10(labels = label_log(digits = 2)) +
  scale_colour_manual(values = agreement_cols,
                      breaks = c("frac_disagree",
                                 "frac_major_disagree",
                                 "frac_minor_disagree"),
                      labels = c("Total",
                                 "Major (\u2265 2 vs \u2265 2)",
                                 "Minor (one outlier)"),
                      name = "Disagreement") +
  scale_fill_manual(values = agreement_cols,
                    breaks = c("frac_disagree",
                               "frac_major_disagree",
                               "frac_minor_disagree"),
                    labels = c("Total",
                               "Major (\u2265 2 vs \u2265 2)",
                               "Minor (one outlier)"),
                    name = "Disagreement") +
  labs(x = expression(tau~"(ha)"), y = "Fraction of cells") +
  theme_AP() +
  theme(legend.position = c(0.65, 0.77))

plot_uncertainty

# Set colors for existential and marginal disagreement -------------------------

exist_cols <- c("share_existence" = "#1f78b4", "share_marginal" = "#e66101")

plot_agreement <- mat_sum[agreement %in% c("share_existence", "share_marginal")] %>%
  ggplot(., aes(x = tau, y = frac_median, colour = agreement, fill = agreement)) +
  geom_ribbon(aes(ymin = frac_lwr, ymax = frac_upr), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = exist_cols,
                      labels = c(share_existence = "Existential",
                                 share_marginal = "Marginal"),
                      name = "Disagreement") +
  scale_fill_manual(values = exist_cols,
                    labels = c(share_existence = "Existential",
                               share_marginal = "Marginal"),
                    name = "Disagreement") +
  scale_x_log10(labels = label_log(digits = 2)) +
  labs(x = expression(tau~"(ha)"), y = "Fraction of cells") +
  theme_AP() +
  theme(legend.position = c(0.6, 0.8))


plot_agreement

# SENSITIVITY ANALYSIS #########################################################

# Sobol' indices ---------------------------------------------------------------

ind <- sobol_indices(matrices = matrices, N = N, params = params, Y
                     = mat[, frac_disagree], R = R, boot = TRUE)

plot_indices3 <- plot(ind) +
  theme_AP() +
  scale_x_discrete(labels = c(exclusion = "exclusion",
                              resolution = "resolution",
                              tau = expression(tau)),
                   guide = guide_axis(n.dodge = 2)) +
  theme(legend.position = c(0.3, 0.8))

plot_indices3

# SCATTERPLOTS OF OUTPUT AGAINS INTPUT #########################################

label_map <- setNames(new_dataset_names, dataset_names)

a3 <- ggplot(mat[!exclusion == "all", c("exclusion", "frac_disagree")],
            aes(exclusion, frac_disagree)) +
  geom_boxplot() +
  stat_summary_bin(fun = "mean", geom = "point",
                   colour = "red", size = 1) +
  scale_x_discrete(labels = label_map) +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  theme_AP() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "exclusion", y = "Fraction disagree")

a3

b3 <- ggplot(mat[1:N, c("resolution", "frac_disagree")],
            aes(resolution, frac_disagree)) +
  geom_boxplot() +
  scale_x_discrete(labels = c(`0.2deg` = "0.2º",
                              `0.4deg` = "0.4º",
                              `1deg` = "1.0º")) +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  stat_summary_bin(fun = "mean", geom = "point",
                   colour = "red", size = 1) +
  theme_AP() +
  labs(x = "resolution", y = "")

b3

c3 <- ggplot(mat[1:N, c("tau", "frac_disagree")],
            aes(tau, frac_disagree)) +
  geom_point(size = 0.1, alpha = 0.1) +
  stat_summary_bin(fun = "mean", geom = "point",
                   colour = "red", size = 1) +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  theme_AP() +
  scale_x_log10(labels = label_log(digits = 2)) +
  labs(x = expression(tau~"(ha)"), y = "")

c3

d3 <- ggplot(mat[1:N, c("weight", "frac_disagree")],
            aes(weight, frac_disagree)) +
  geom_point(size = 0.1, alpha = 0.1) +
  stat_summary_bin(fun = "mean", geom = "point",
                   colour = "red", size = 1) +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  theme_AP() +
  labs(x = "weight", y = "")

d3


## ----plot_weights, dependson="uasa_on_weights", fig.height=3.8, fig.width=5.5-------------------

# MERGE UA/SA PLOTS ############################################################

top <- plot_grid(plot_uncertainty, plot_agreement, plot_indices3,
                 ncol = 3, labels = "auto", rel_widths = c(0.4, 0.3, 0.3))
bottom <- plot_grid(a3, b3, c3, d3, ncol = 4, rel_widths = c(0.35, 0.25, 0.2, 0.2), labels = "d")
plot_grid(top, bottom, ncol = 1)


## ----tau_max, dependson="uasa"------------------------------------------------------------------

# COMPUTE TAU MAX ##############################################################
################################################################################

# Tau max is the largest irrigation threshold at which at least one dataset
# classifies the cell as essentially unirrigated  while another classifies
# it as irrigated. Tau max thus quantifies how large the ontological disagreement
# on irrigation is. We also set NA to two cases:

# - Universally non-irrigated cells: all datasets report areas <= threshold.
# - Universally irrigated cells: all datasets report areas >threshold.

# NA cells are cells for which there is agreement irrigation exists or it is absent.

## Compute tau_max_results -----------------------------------------------------

tau_max_results <- rbindlist(
  parallel::mclapply(resolutions, function(res) {

    dt_res_long <- dt[resolution == res]

    # datasets present at this resolution
    ds_res <- sort(unique(dt_res_long$dataset))

    # resolution-specific detectability thresholds (in ha)----------------------

    zero_tol_ha <- switch(res,
      "0.2deg" = 2450, # 1% of cell at 0.2°
      "0.4deg" = 9850, # 1% of cell at 0.4°
      "1deg" = 61600, # 1% of cell at 1°
      stop("Unknown resolution: ", res))

    rbindlist(
      lapply(dataset_vec, function(scn) {

        # which datasets to use -----------------------------------------------

        if (scn == "all") {
          use_ds <- ds_res
        } else {
          use_ds <- setdiff(ds_res, scn)
        }

        # if only one dataset left, skip --------------------------------------

        if (length(use_ds) < 2L) return(NULL)

        # reshape to wide ----------------------------------------------------


        dt_res_wide <- dcast(
          dt_res_long[dataset %in% use_ds],
          lon + lat + country + code + continent ~ dataset,
          value.var = "mha", fill = 0
        )

        tau_max_dt <- tau_grid[,
          compute_tau_max_existential_fun(dt_res_wide, use_ds, tau_mha, tau_ha,
                                          zero_tol_ha = zero_tol_ha)
        ]

        tau_max_dt[, `:=`(resolution = res, scenario = scn)]

        tau_max_dt
      }),
      use.names = TRUE, fill = TRUE
    )
  }, mc.cores = ncores),
  use.names = TRUE, fill = TRUE
)


# Add the countries ------------------------------------------------------------

chunk_size <- 500000L
idx <- split(seq_len(nrow(tau_max_results)),
             ceiling(seq_len(nrow(tau_max_results)) / chunk_size))

tau_max_results[, country:= NA_character_]
for (i in idx) {
  tau_max_results[i, country:= lonlat_to_country(lon, lat, countries_sf)]
}

# Add the continents -----------------------------------------------------------

country_to_continent(tau_max_results)

# Retrieve the round with all datasets and exclude the rest --------------------

tau_max_all <- tau_max_results[scenario == "all"]


## ----define_plots_taumax, dependson="tau_max"---------------------------------------------------

# DEFINE PLOTS TAU MAX #########################################################

# Bin the tau max results ------------------------------------------------------

tau_max_results[, det_thresh:= fcase(
  resolution == "0.2deg", 2450,
  resolution == "0.4deg", 9850,
  resolution == "1deg", 61600)]

tau_max_results[, tau_bin:= fcase(
  is.na(tau_max_ha), "Agreement",
  tau_max_ha <= 2 * det_thresh, "Low disagreement",
  tau_max_ha <= 5 * det_thresh, "Moderate disagreement",
  tau_max_ha <= 10 * det_thresh, "High disagreement",
  default = "Extreme disagreement"
)]

tau_max_results[, tau_bin:= factor(tau_bin,levels = c("Agreement",
                                                       "Low disagreement",
                                                       "Moderate disagreement",
                                                       "High disagreement",
                                                       "Extreme disagreement"),
                                    ordered = TRUE)]

# Proportion of cells disagreeing at 1% threshold-------------------------------

tau_max_results[, .N, .(tau_bin, resolution)] %>%
  .[, fraction := N / sum(N), by = resolution] %>%
  .[!tau_bin == "Agreement", sum(fraction), resolution]

thresh_df <- unique(tau_max_results[!is.na(det_thresh), .(resolution, det_thresh)])

# ANALYSIS OF DATA -------------------------------------------------------------

# Fraction of high and extreme ontological disagreement
# across resolutions in each continent -----------------------------------------

dt_all <- tau_max_results[scenario == "all"]
dt_all[, N_total := .N, by = .(resolution, continent)]
frac_by_rcb <- dt_all[, .(frac = .N / unique(N_total)),
                      .(resolution, continent, tau_bin)]

frac_by_rcb[, .(min = min(frac), max = max(frac)), .(continent, tau_bin)] %>%
  .[tau_bin %in% c("Extreme disagreement", "High disagreement")] %>%
  .[order(-min)] %>%
  .[, .(min = min(min), max = max(max)), .(continent)]

# distribution of detectability limits

tau_max_results[, .N, .(tau_max_ha, resolution, scenario)]

tmp <- copy(tau_max_results)
tmp[is.na(tau_max_ha), tau_max_ha:= 0.1]

plot_ecdf <- ggplot(tmp, aes(x = tau_max_ha, color = resolution)) +
  stat_ecdf(size = 1) +
  scale_x_log10(labels = label_log(digits = 2)) +
  scale_color_manual(values = res_palette) +
  labs(x = expression(tau[max]~"(ha)"), y = "Fraction grid cells",
       color = "Resolution") +
  theme_AP() +
  theme(legend.position = c(0.65, 0.25))

plot_ecdf

# How much of the world is perfectably detectable vs fundamentally ambiguous?

plot_stacked <- tau_max_results[, .N, by = .(resolution, tau_bin)] %>%
  .[, frac := N / sum(N), by = resolution] |>
  ggplot(aes(resolution, frac, fill = tau_bin)) +
  geom_col() +
  scale_y_continuous() +
  scale_fill_manual(values = tau_palette, name = expression(tau[max])) +
  labs(y = "Fraction grid cells", x = NULL, fill = expression(tau[max])) +
  theme_AP() +
  theme(legend.position = "none")

plot_stacked

# Boxplots of tau max per resolution -------------------------------------------

plot_box <- ggplot(tau_max_results[!is.na(tau_max_ha) & tau_max_ha > 0],
                   aes(resolution, tau_max_ha, fill = resolution)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_hline(data = thresh_df, aes(yintercept = det_thresh, color = resolution),
             linetype = "dashed", linewidth = 0.7, show.legend = FALSE) +
  scale_y_log10(labels = scales::label_log(digits = 2),
                name   = expression(tau[max]~"(ha)")) +
  scale_fill_manual(values = res_palette, guide = "none") +
  scale_color_manual(values = res_palette, guide = "none") +
  labs(
    x = NULL,
  ) +
  theme_AP() +
  theme()

plot_box


## ----merge_tau_max, dependson=c("tau_max", "define_plots_taumax"), fig.height=3, fig.width=5----

# MERGE TAU MAX PLOTS ##########################################################

top <- plot_grid(plot_ecdf, plot_stacked, plot_box, ncol = 1,
                 rel_heights = c(0.35, 0.35, 0.3),
          labels = "auto")

top


## ----plot_map, dependson=c("tau_max", "merge_tau_max", "define_plots_taumax")-------------------

# The map ######################################################################

plot_raster <- ggplot(tau_max_results, aes(lon, lat, fill = tau_bin)) +
  geom_raster() +
  coord_equal(expand = FALSE) +
  facet_wrap(~ resolution, ncol = 1) +
  scale_fill_manual(values = tau_palette, name = expression(tau[max])) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_AP() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 11),
        legend.position = "top") +
  labs(x = "longitude", y = "latitude") +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

plot_raster


## ----merge_plot_tau_max_map, dependson=c("plot_map", "merge_tau_max"), fig.width=5.5------------

# MERGE PLOTS ##################################################################

plot_grid(top, plot_raster, rel_widths = c(0.3, 0.7), ncol = 2, labels = c("", "d"))


## ----datasets_tau_max, dependson="tau_max", fig.height=1.5, fig.width=5-------------------------

# WHICH DATASETS DRIVE DEEP DISAGREEMENT? --------------------------------------

# Check leave one-out scenarios ------------------------------------------------

loo_stats <- tau_max_results[resolution == "0.2deg",
                             .(median_tau = median(tau_max_ha, na.rm = TRUE)),
                             scenario] %>%
  .[order(-median_tau)]

loo_stats

# Most datasets do not drive median existential disagreement, but three datasets
# reduce the median diagreement when removed: GIAM, LUH2, Nagaraj (c. 20% reduction)
# . They are the ones contributing most disagreement to existential contradictions
# at the global median level (they most often say "there is irrigation" while
# others say "there is no irrigation" or the other way around). So disagreement
# is systemic but three datasets aplify it further.

# Full stats -------------------------------------------------------------------

loo_stats_full <- tau_max_results[, .(
    median_tau = median(tau_max_ha, na.rm = TRUE),
    p75_tau = quantile(tau_max_ha, 0.75, na.rm = TRUE),
    p90_tau = quantile(tau_max_ha, 0.90, na.rm = TRUE),
    share_gt_1000 = mean(tau_max_ha >= 1000, na.rm = TRUE)),
    scenario] %>%
  .[order(-p90_tau)]

loo_stats_full

# In half of all grid cells worldwide, at least one dataset reports >=1000 ha
# of irrigation while another reports none. Removing most datasets does not
# change anything. Removing GIAM, Nagaraj and LUH2 drops the median disagreement
# and the extreme values, as they particularly increase existential contradictions
# in the upper half of the distribution.

# In which regions is each dataset structurally responsible for
# deep existential contradictions? -----------------------------------------------

# Compute spatial delta -----------

tau_all <- tau_max_results[scenario == "all",.(lon, lat, resolution,
                                               tau_all = tau_max_ha)]

# Keep only datasets of interest--

tau_loo <- tau_max_results[scenario %in% c("giam", "luh2", "nagaraj"),
                           .(lon, lat, country, continent, resolution, scenario,
                             tau_loo = tau_max_ha)]

# Merge ---------------------------

tau_delta <- merge(tau_loo, tau_all, by = c("lon", "lat", "resolution"),
                   all.x = TRUE)

# Compute difference --------------

tau_delta[, delta:= tau_all - tau_loo]

# Let us define "strong amplification" of differences --------------------------
# Differences larger than 1000 ha.

tau_delta[, strong:= delta >= 1000]

# Compute overlap --------------------------------------------------------------

influence_matrix <- dcast(tau_delta[resolution == "0.2deg"],
                          lon + lat ~ scenario, value.var = "strong",
                          fill = FALSE)

influence_matrix[, n_datasets:= rowSums(.SD), .SDcols = c("giam", "luh2", "nagaraj")]

table(influence_matrix$n_datasets)

# 0: deep disagreement in the cell not driven by these three.
# 1: dataset-specific amplification.
# 2: partial overlap.
# 3: same hotspot amplified by all three.

# Since there are zero cells where two datasets simultaneously aplify deep
# disagreement (> 1000 ha) or all three do, it means that their influence is
# spatially distinct. 22,615 cells show strong amplification but always due
# to exactly one dataset.

# Where does amplification concentrate? ----------------------------------------

regional_amp <- tau_delta[resolution == "0.2deg",
                          .(median_delta = median(delta, na.rm = TRUE),
                            share_strong = mean(delta >= 1000, na.rm = TRUE)),
                          .(continent, scenario)] %>%
  .[order(-median_delta)]


regional_amp
# Table shows which continents each dataset most strongly affects:
# - Europe: GIAM dominates existential disagreement in Europe, LUH2 barely affects it.
# - Asia: all three datasets contribute, LUH2 contributes less existential disagreement
# - Ocenia: GIAM contributes the largest existential disagreement.
# - Americas: Nagaraj.
# - Africa: Overall amplification much lower, all contribute equally.
# So overall there is no global problem dataset! And since they all show
# median delta of zero everywhere, these datasets do not globally modify the
# disagreement but create region-specific deep contradictions.

ggplot(regional_amp, aes(scenario, share_strong)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~continent, ncol = 5) +
  scale_y_continuous(breaks = breaks_pretty(n = 2)) +
  theme_AP() +
  labs(x = "", y = "Fraction cells")


# Which dataset drives the 22,615 cells? ---------------------------------------

influence_matrix[, dominant := fifelse(giam, "GIAM",
                                       fifelse(luh2, "LUH2",
                                               fifelse(nagaraj, "Nagaraj",
                                                       NA_character_)))]

table(influence_matrix$dominant)


## ----tau_weighted, dependson=c("uasa", "tau_max"), fig.height=2, fig.width=2.3------------------

# DETECTABILITY OF IRRIGATED AREA ACROSS MAPS ##################################
################################################################################

# UNWEIGHTED: ALL CELLS COUNT EQUALLY ##########################################
# QUESTION ANSWERED: HOW MANY GRID CELLS DISAGREE? -----------------------------

# Make dataset wide and retrieve only the 0.2deg resolution one ----------------

dt_wide <- dt[resolution == "0.2deg"] %>%
  dcast(., lon + lat + country + code + continent ~ dataset, value.var = "mha")

# Retrieve 0.2deg resolution from tau max dataset and remove unwanted cols -----

tau_max_results_clean <- tau_max_results[resolution == "0.2deg" & scenario == "all"]
tau_max_results_clean[, resolution:= NULL]
tau_max_results_clean[, scenario:= NULL]

# Merge ------------------------------------------------------------------------

dt2 <- merge(dt_wide, tau_max_results_clean[, .(lon, lat, tau_max_mha)],
             by = c("lon", "lat"), all.x = TRUE)

setnames(dt2, dataset_names, new_dataset_names)

# Convert to ha for easier interpretation --------------------------------------

dt2[, tau_max_ha:= tau_max_mha * 1e6]

# Which percentage of the world is in disagreement on whether there is
# or there is not irrigation? --------------------------------------------------

frac_cells_nonid <- dt2[, mean(!is.na(tau_max_mha))]
frac_cells_nonid

summary(dt2$tau_max_ha)

# Minimum = 398 ha, the detectability threshold for 0.2. No disagreement below it.
#
# Median disagreement (50%) is eight times the detectability threshold (400 ha)
# 75% of disagreeing cells exceed 2.5x detectability.
# Half exceed 8x detectability

# Even after imposing a conservative 1% rule some cells show contradictions at
# tens of thousands of hectares; these are large irrigation landscapes, not
# fragmented patches.

# Existential contradictions occur at scales that are agronomically and
# hydrologically meaningful and they are not confined to marginal lands. Include
# cells where irrigation is either extensive or completely absent depending on
# the dataset.Ddisagreement is not about "how much" irrigation exists but about
# whether large irrigation systems exist at all in certain locations.



## ----weighted, dependson="tau_weighted"---------------------------------------------------------

# WEIGHTED: CELLS WITH LARGER IRRIGATION COUNT MORE ############################
# QUESTION: HOW MUCH IRRIGATED LAND LIES IN DISAGREEMENT? ----------------------

# Keep only cells with existential disagreement -------------------------------
# (tau_max_ha non-NA)

dt2_disagree <- dt2[!is.na(tau_max_ha)]

# Sanity check: no NA tau_max_ha left
summary(dt2_disagree$tau_max_ha)
# should show NA's = 0

# Compute total and non-identifiable area per dataset -------------------------

area_nonid <- rbindlist(lapply(new_dataset_names, function(d) {

  # cells where this dataset has irrigation AND there is disagreement

  dt_d <- dt2[get(d) > 0, .(tau_max_ha, area = get(d))]

  # total irrigated area (for this dataset) in cells with disagreement

  total_area <- dt_d[, sum(area, na.rm = TRUE)]

  # irrigated area in cells where tau_max_ha is not NA; that is, in cells
  # that disagree

  nonid_area <- dt_d[!is.na(tau_max_ha), sum(area, na.rm = TRUE)]

  data.table(dataset = d, total_area_mha = total_area,
             nonid_area_mha = nonid_area, share_nonid = nonid_area / total_area)
}))

area_nonid

# Depending on dataset, between 50-70% of all global irrigated areas
# lie in areas where we do not even agree on whether there is irrigation at all.

# How much tolerance do we need to make datasets agree?-------------------------

tau_diag <- rbindlist(lapply(new_dataset_names, function(d) {

  dt_d <- dt2_disagree[get(d) > 0, .(tau_max_ha, area = get(d))]

  if (nrow(dt_d) == 0) {
    return(data.table(dataset = d,
                      tau_50 = NA_real_,
                      tau_80 = NA_real_,
                      tau_90 = NA_real_,
                      area_share_1000 = NA_real_  ))
  }

  # area-weighted empirical CDF of tau_max_ha (for this dataset in disagreeing cells)

  tmp <- dt_d[, .(area = sum(area)), by = tau_max_ha][order(tau_max_ha)]
  tmp[, cum_area  := cumsum(area)]
  tmp[, cum_share := cum_area / sum(area)]

  # tau_50: minimum tau required to make 50% of irrigated area identifiable
  # tau_80: tau required for 80%
  # tau_90: tau required for 90%

  tau_50 <- tmp[cum_share >= 0.5][1, tau_max_ha]
  tau_80 <- tmp[cum_share >= 0.8][1, tau_max_ha]
  tau_90 <- tmp[cum_share >= 0.9][1, tau_max_ha]

# share of area that is identifiable under tau_threshold_ha (e.g. 1000 ha)

  area_share_1000 <- tmp[tau_max_ha <= 1000, fifelse(.N == 0, 0,
                                                     max(cum_share, na.rm = TRUE))]

  data.table(dataset = d,
             tau_50 = tau_50,
             tau_80 = tau_80,
             tau_90 = tau_90,
             area_share_1000 = area_share_1000)
}))

tau_diag

# tau_50 ranges from 8,000 to 16,000 ha.
# tau_80 is 25,000–40,000 ha.
# tau90 is often 31,000–40,000 ha.
#
# At 0.2° resolution (40,000 ha per cell):
# tau_50 equals 12,000–16,000 ha; that is 30–40% of the grid cell.
# tau_90 equals 30,000–40,000 ha; most of the grid cell area.
#
# That means that Tto reconcile 90% of irrigated area in disagreeing cells one
# must tolerate differences approaching the size of an entire grid cell.

# area_share_1000: Under a tolerance of 1000 ha, only 0.6%–2.7% of irrigated area
# in disagreeing cells becomes identifiable.

# PARAGRAPH?
# Area-weighted quantiles of tau_maxreveal that resolving even half of irrigated
# area located in disagreeing cells would require tolerating discrepancies of
# 8,000–16,000 ha. resolving 90% would require accepting differences
# exceeding 30,000 ha in most datasets. Under a tolerance of 1,000 ha,
# less than 3% of irrigated area in disagreeing cells becomes identifiable.
# These values indicate that disagreement persists at scales corresponding to
# substantial irrigation systems rather than marginal patches.


## ----country_level, dependson="tau_weighted"----------------------------------------------------

# DOES THE DIFFERENT DATASETS IDENTIFY THE SAME TOP IRRIGATION HOTSPOTS? #######

# Identify top 10% hotspots per dataset ----------------------------------------

hotspots <- lapply(new_dataset_names, function(d) {

  tmp <- copy(dt2)
  tmp[, value:= get(d)]

  # rank cells---------------------------

  tmp <- tmp[order(-value)]
  n_hot <- ceiling(0.10 * nrow(tmp))

  tmp[, hotspot:= FALSE]
  tmp[1:n_hot, hotspot:= TRUE]

  tmp[, .(lon, lat, hotspot, dataset = d)]
})

# name the list elements by dataset --------------------------------------------

names(hotspots) <- sapply(hotspots, function(x) unique(x$dataset))

# extract hotspot cell IDs per dataset -----------------------------------------

hotspot_sets <- lapply(hotspots, function(dt) {
  dt[hotspot == TRUE, .(lon, lat)]})

# Country rankings per dataset -------------------------------------------------

country_totals <- rbindlist(lapply(new_dataset_names, function(d) {
  dt2[, .(irrigated_area = sum(get(d), na.rm = TRUE)), country] %>%
    .[, dataset:= d]}))

jaccard_cells <- function(dt1, dt2) {
  set1 <- paste(dt1$lon, dt1$lat)
  set2 <- paste(dt2$lon, dt2$lat)

  length(intersect(set1, set2)) / length(union(set1, set2))
}

# Pairwise jacqard similarity --------------------------------------------------

datasets <- names(hotspot_sets)

jaccard_results <- rbindlist(combn(datasets, 2, simplify = FALSE,
                                   FUN = function(pair) {
    d1 <- pair[1]
    d2 <- pair[2]

    data.table(dataset_1 = d1, dataset_2 = d2,
               jaccard = jaccard_cells(hotspot_sets[[d1]], hotspot_sets[[d2]]))
  })
)

jaccard_results[order(jaccard)]
jaccard_results[, .(min_jaccard = min(jaccard),
                    median_jaccard = median(jaccard),
                    mean_jaccard = mean(jaccard),
                    max_jaccard = max(jaccard))]

# results:
# 1: both datasets identify the same hotspots
# 0: datasets disagree completely
# 0.5: half overlap, half disagreement

# No pair of global irrigation datasets agrees strongly on where
# irrigation hotspots are.

# Weak agreement (0.35 - 0.45):
# -----------------------------
# one dataset says "this is a core irrigated region",
# the other dataset says "no, it is not".

# Moderate agreement (0.5-0.65)
# -------------------------------
# one third to one half of hotspots do not overlap

# High agreement (0.7 - 0.8)
# -------------------------------
# High agreement only between specific subsets because of shared lineage
# (GMIA and Meier; MIRCA2000 and Meier, etc)


## ----national_rankings, dependson=c("tau_weighted", "country_level")----------------------------

# NATIONAL RANKINGS USING ALL CELLS ############################################

# Country totals & ranks .------------------------------------------------------

country_totals_all <- compute_country_totals(dt2, new_dataset_names)

# Pairwise rank correlations ---------------------------------------------------

rank_corr_all <- compute_rank_correlations(country_totals_all)
rank_corr_all[order(-spearman)]

# Pairwise top-N overlap (e.g. top 20 irrigators) ------------------------------

topN_overlap_all <- compute_topN_overlap(country_totals_all,
                                         dataset_names = new_dataset_names,
                                         N = 20)
topN_overlap_all[order(-jaccard_topN)]

# NATIONAL RANKINGS USING ONLY "IDENTIFIABLE" CELLS ############################

# Only areas that agree --------------------------------------------------------

dt_ident <- dt2[is.na(tau_max_ha)]

# Country totals & ranks on identifiable cells only ----------------------------

country_totals_ident <- compute_country_totals(dt_ident, new_dataset_names)

# Pairwise rank correlations (identifiable-only) -------------------------------

rank_corr_ident <- compute_rank_correlations(country_totals_ident)
rank_corr_ident[order(-spearman)]

# Pairwise top-N overlap (identifiable-only) -----------------------------------

topN_overlap_ident <- compute_topN_overlap(country_totals_ident,
                                           dataset_names = new_dataset_names,
                                           N = 20)
topN_overlap_ident[order(-jaccard_topN)]


## ----rankings_change, dependson=c("national_rankings", "country_level"), fig.width=4------------

# HOW DO RANKINGS CHANGE WHEN NON-IDENTIFIABLE CELLS ARE MASKED? ###############

# How much does each country's rank shift for each dataset? --------------------

rank_change <- merge(
  country_totals_all[, .(country, dataset, rank_all = rank)],
  country_totals_ident[, .(country, dataset, rank_ident = rank)],
  by = c("country", "dataset"),
  all.x = TRUE) %>%
  .[, rank_shift := rank_ident - rank_all]

rank_change

# PLOT RESULTS #################################################################

# Arrange data -----------------------------------------------------------------

rank_corr_all[, pair:= paste(dataset_1, dataset_2, sep = "–")]
rank_corr_ident[, pair:= paste(dataset_1, dataset_2, sep = "–")]

rank_corr_long <- rbind(rank_corr_all[, .(pair, spearman, mode = "All cells")],
                        rank_corr_ident[, .(pair, spearman, mode = "Identifiable cells")])

# order pairs by Spearman (all cells) so the y-axis is readable ----------------

pair_order <- rank_corr_all[order(spearman), pair]
rank_corr_long[, pair := factor(pair, levels = pair_order)]

# Plot rank correlation --------------------------------------------------------

p_rankcorr <- ggplot(rank_corr_long, aes(spearman, pair)) +
  geom_segment(data = dcast(rank_corr_long, pair ~ mode, value.var = "spearman"),
               aes(x = `All cells`, xend = `Identifiable cells`,
                   y = pair, yend = pair),
               inherit.aes = FALSE,
               colour = "grey70", linewidth = 0.3) +
  geom_point(aes(colour = mode), size = 2, alpha = 0.5) +
  scale_x_continuous(limits = c(0.4, 1.0)) +
  scale_color_manual(values = c("Identifiable cells" = "#d73027",
                                "All cells" = "black")) +
  labs(x = "Rank correlation of \n national irrigation rankings",
       y = "",
       colour = NULL) +
  theme_AP() +
  theme(legend.position = c(0.4, 0.85))

p_rankcorr

# Only for top N countries -----------------------------------------------------

topN_overlap_all[, pair:= paste(dataset_1, dataset_2, sep = "–")]
topN_overlap_ident[, pair:= paste(dataset_1, dataset_2, sep = "–")]

topN_long <- rbind(topN_overlap_all[, .(pair, jaccard_topN, mode = "All cells")],
                   topN_overlap_ident[, .(pair, jaccard_topN, mode = "Identifiable cells")])

# order by All-cells Jaccard for consistency -----------------------------------

pair_order_topN <- topN_overlap_all[order(jaccard_topN), pair]
topN_long[, pair:= factor(pair, levels = pair_order_topN)]

# Plot -------------------------------------------------------------------------

p_topN <- ggplot(topN_long, aes(x = jaccard_topN, y = pair)) +
  geom_segment(data = dcast(topN_long, pair ~ mode, value.var = "jaccard_topN"),
               aes(x = `All cells`, xend = `Identifiable cells`,
                   y = pair, yend = pair), inherit.aes = FALSE, colour = "grey70",
               linewidth = 0.3) +
  geom_point(aes(colour = mode), size = 2, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Identifiable cells" = "#d73027",
                                "All cells" = "black")) +
  labs(x = "Jaccard similarity of \n top 20 countries",
       y = "", colour = NULL) +
  theme_AP() +
  theme(legend.position = "none")

p_topN


## ----merge, dependson="rankings_change", fig.width=5.5------------------------------------------

# MERGE ########################################################################

# (a) Spearman rank correlation of national total irrigated area rankings for all
# pairwise combinations of datasets. Black points show correlations computed across
# all grid cells, while red points show correlations restricted to identifiable
# cells (i.e., excluding cells where datasets fundamentally disagree on the
# presence of irrigation). Correlations are generally high (0.7–1.0),
# indicating broad agreement in country-level rankings. In some cases removing
# disputed cells increases the correlation (red dot to the right of black dot),
# while in others it does not (red dot to the left of black dot)
#
# (b) Jaccard similarity of the top 20 irrigated countries for each dataset pair.
# Black points again denote results using all cells and red points identifiable cells
# only. Most of the time focusing only on identifiable cells increases the
# similarity; meaning that the top 20 become more similar across datasets. In other
# words, once disputed cells are removed there is a more coherent "core" of
# globally dominant irrigation countries. Possibly this is because the strongest
# disagreement disproportionately affect marginal or mid-rank countries that
# occassionally enter or leave the top 20 depending on dataset.

# SO: removing disagreeing cells stabilizes who belongs in the top twenty, but not
# necessarily the full ranking hierarchy across al coutnries.

plot_grid(p_rankcorr, p_topN, ncol = 2, labels = "auto")


## ----shift_country_ranks, dependson=c("rankings_change", "national_rankings"), fig.height=5.4, fig.width=5----

# SHIFT IN COUNTRY RANKS #######################################################

# Total rank shift -------------------------------------------------------------

rank_change[, abs_shift:= abs(rank_shift)]
rank_change_clean <- rank_change[!is.na(abs_shift)]

# Histogram of rank shift ------------------------------------------------------

p_rankshift_hist <- ggplot(rank_change_clean, aes(abs_shift)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left") +
  labs(x = "|Change in national irrigation rank|",
       y = "Number of country–dataset combinations") +
  theme_AP()

p_rankshift_hist

## Cumulative share ------------------------------------------------------------

rank_shift_stats <- rank_change_clean[, .(share_ge_1  = mean(abs_shift >= 1),
                                          share_ge_5  = mean(abs_shift >= 5),
                                          share_ge_10 = mean(abs_shift >= 10),
                                          share_ge_20 = mean(abs_shift >= 20),
                                          share_ge_50 = mean(abs_shift >= 50))]

rank_shift_stats

# 98.4% of country–dataset combinations change rank by at least 1 position.
# No country keeps the same rank once we enforce identifiability.
#
# 90% of countries move by at least 5 rank positions
#
# 71% of country–dataset combinations move by 10 ranks or more

# FOCUS ON THE SO-CALLED "TOP IRRIGATORS" #####################################

number_countries <- 20

# Countries consistently classified as top irrigators across datasets ----------

top10_all <- country_totals_all[order(rank), .SD[1:number_countries], dataset] %>%
  .[, .N, by = country] %>%
  .[!country == ""] %>%
  .[N > 5]

top10_all

# Also with the median ---------------------------------------------------------

top10_median <- country_totals_all[, .(median_rank = median(rank)), country] %>%
  .[order(median_rank)] %>%
  .[1:number_countries]

top10_median

# See change of ranks in countries when only confirmed irrigation is
# accounted for ----------------------------------------------------------------

rank_focus <- merge(country_totals_all[, .(country, dataset, rank_all = rank)],
                    country_totals_ident[, .(country, dataset, rank_ident = rank)],
                    by = c("country", "dataset"), all = TRUE)

rank_focus[, rank_shift:= rank_ident - rank_all]
rank_focus

# Summarize --------------------------------------------------------------------

rank_focus_summary <- rank_focus[!is.na(rank_ident), .(median_shift = median(rank_shift, na.rm = TRUE),
                                     max_drop = max(rank_shift, na.rm = TRUE),
                                     max_rise = min(rank_shift, na.rm = TRUE),
                                     share_drop5 = mean(rank_shift >= 5, na.rm = TRUE),
                                     share_drop10 = mean(rank_shift >= 10, na.rm = TRUE)),
                                 country] %>%
  .[order(-median_shift)]

# Focus on top 20 irrigation hotspots ------------------------------------------

focus_countries <- top10_median$country
rank_focus_summary[country %in% focus_countries]

# Positive rank shift: worse rank if only confirmed cells are accounted for.
# Negative rank shift: better rank if only confirmed cells are accounted for.

# - median_shift: median of rank_shift across datasets.
# - max_drop: maximum worsening (largest positive rank_shift).
# - max_rise: maximum improvement (most negative rank_shift).
# - share_drop5: share of dataset cases with rank_shift >= 5.
# - share_drop10: share with rank_shift >= 10.

# VIETNAM, BANGLADESH, THAILAND, INDONESIA: in 90% of cases, Vietnam's rank gets
# at least 10 places worse. Median drop is c. 40 places. Bangladesh, Indonesia,
# Thailand show similar pattern; median shifts of 26-27 ranks and 90% of cases
# with drops larger than 90 ranks.

# CHINA, INDIA, PAKISTAN, USA: Global rank is invariant to whether we use all cells
# or only agreement cells. Irrigated area dominance is robust.

# ITALY, EGYPT, UZBEKISTAN, IRAQ, SPAIN, KAZAKHSTAN: they
# improve their rank when we restrict to confirmed cells; because when we remove
# ambiguous cells in other countries, these countries become more prominent.

# Focus on the countries whose ranking drops the most --------------------------

rank_focus_summary[order(-median_shift)][1:number_countries]

# Plot the results -------------------------------------------------------------

plot_shift_ranks <- rank_focus[country %in% focus_countries] %>%
  .[, .(country, dataset,
        `All irrigated cells` = rank_all,
        `Confirmed irrigation` = rank_ident)] %>%
  melt(., measure.vars = c("Confirmed irrigation", "All irrigated cells")) %>%
  ggplot(., aes(dataset, value, color = variable,
                group = interaction(country, variable))) +
  geom_line(linewidth = 1) +
  scale_y_reverse() +
  scale_color_manual(values = c("Confirmed irrigation" = "#d73027",
                                "All irrigated cells" = "black"),
                     name = "") +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  facet_wrap(~country) +
  labs(y = "National irrigation rank", x = NULL) +
  coord_flip() +
  theme_AP() +
  theme(legend.position = "top")

plot_shift_ranks

# Dashed (red) line much lower than solid (black)
# - country drops in rank when ambiguity is removed
# - it looks large only because of non-identifiable cells
#
# Red equals black
# - rank is robust to identifiability
# - irrigation signal is spatially coherent
#
# Red above black
# - country rises in rank when ambiguity is removed
# - it was previously masked by others’ ambiguous irrigation

# If Red and black lines cross across datasets:
# The ordering of countries changes
# Rankings are dataset-dependent
# There is no stable hierarchy
# Crossings are visual proof of rank instability.



## ----plot_tileplot, dependson="tau_weighted", fig.height=2.7, fig.width=3.5---------------------

# COMPUTE AREA RETAINED AND LOST OF TOP 20 COUNTRIES ###########################

country_loss <- rbindlist(lapply(new_dataset_names, function(d) {

  dt_country <- dt2[country %in% top10_all$country]

  dt_country[, .(dataset = d,
                 total_area = sum(get(d), na.rm = TRUE),
                 retained_area = sum(get(d)[is.na(tau_max_ha)], na.rm = TRUE),
                 lost_area = sum(get(d)[!is.na(tau_max_ha)], na.rm = TRUE)),
             country]
}))

country_loss[, share_lost:= lost_area / total_area]

# Compute mean loss per country across datasets to reorder ---------------------

country_order <- country_loss[, .(mean_share_lost = mean(share_lost, na.rm = TRUE)),
  by = country][order(-mean_share_lost)]  # decreasing order

# Set factor levels in that order-----------------------------------------------

country_loss[, country_f:= factor(country, levels = rev(country_order$country))]

# Plot tileplot ----------------------------------------------------------------

plot_tileplot <- ggplot(country_loss, aes(dataset, country_f, fill = share_lost)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(colours = c("white", "gold", "orange", "red"),
                       values  = c(0, 0.5, 0.8, 1),
                       limits  = c(0, 1),
                       name = "Share of irrigated\narea lost") +
  scale_y_discrete(name = NULL) +
  theme_AP() +
  labs(x = "", y = "") +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "right")

plot_tileplot


## ----fraction_lost_k, dependson="uasa", fig.height=2, fig.width=3-------------------------------

# COLLAPSE ASSESSMENT ACROS A CONTINUUM OF AGREEMENT ###########################

dt_share_lost <- dcast(dt, lon + lat + country + code + continent + resolution ~ dataset,
      value.var = "mha")

# Define resolution specific thresholds ----------------------------------------

zero_tol_mha_fun <- function(res) {
  switch(res,
    "0.2deg" = 2450  / 1e6, # 0.0004 Mha
    "0.4deg" = 9850 / 1e6, # 0.0016 Mha
    "1deg" = 61600 / 1e6,  # 0.01 Mha
    stop("Unknown resolution: ", res)
  )
}

# Make copy just in case -------------------------------------------------------

dt_pres <- copy(dt_share_lost)

# Apply resolution-specific thresholds -----------------------------------------

dt_pres[, (dataset_names) := {
    tol_mha <- zero_tol_mha_fun(resolution[1])
    lapply(.SD, function(x) as.integer(x >= tol_mha))
    },
    resolution,
    .SDcols = dataset_names]

# How many datasets say "irrigated" per cell -----------------------------------

dt_pres[, n_pos:= rowSums(.SD, na.rm = TRUE), .SDcols = dataset_names]

n_datasets <- length(dataset_names)  # should be 10

# Melt for k-of-n agreement ----------------------------------------------------

grid_long <- melt(dt_pres, id.vars = c("lon","lat","country","code",
                    "continent","resolution","n_pos"),
                  measure.vars = dataset_names,
                  variable.name = "dataset",
                  value.name = "present")

# Run calculations for multiple k (10 down to 5) -------------------------------

ks <- n_datasets:5
loss_by_k <- rbindlist(lapply(ks, compute_loss_for_k), use.names = TRUE)

# Summarize results ------------------------------------------------------------

country_loss_by_k <- loss_by_k[, .(min  = min(share_lost,  na.rm = TRUE),
                                   max  = max(share_lost,  na.rm = TRUE),
                                   mean = mean(share_lost, na.rm = TRUE)),
                               .(country, k)] %>%
  .[, allowed_disagree:= 10 - k]

# PLOT RESULTS: THE GLOBAL COLLAPSE CURVE ######################################

global_k <- country_loss_by_k[, .(mean_loss   = mean(mean, na.rm = TRUE),
                                  median_loss = median(mean, na.rm = TRUE)), k] %>%
  .[order(-k)] %>%
  .[, allowed_disagree:= 10 - k]

print(global_k)

plot_global_k <- global_k %>%
  melt(., measure.vars = c("median_loss", "mean_loss")) %>%
  ggplot(., aes(k, value, color = variable)) +
  geom_line() +
  geom_point()  +
  scale_color_discrete(labels = c(median_loss = "median", mean_loss = "mean"),
                       name = "") +
  labs(x = "Datasets agreeing", y = "Fraction irrigation lost") +
  theme_AP() +
  scale_x_reverse() +
  theme(legend.position = c(0.7, 0.8))

plot_global_k

# COUNTRY-LEVEL COLLAPSE PROFILES ###############################################

# Optional: keep only countries with mean loss < 1 at k=5

collapse_profiles <- country_loss_by_k[country %in% top10_all$country] %>%
  .[, allowed_disagree:= 10 - k]

print(collapse_profiles)

plot_collapse_profiles <- ggplot(collapse_profiles, aes(k, mean, group = country,
                              color = country)) +
  scale_color_discrete(name = "") +
  geom_line() +
  scale_x_reverse() +
  labs(x = "Datasets agreeing", y = "Fraction irrigation lost") +
  theme_AP()

plot_collapse_profiles

# Global distribution by k -----------------------------------------------------

plot_curves_country <- ggplot(country_loss_by_k, aes(mean, color = factor(k), group = k)) +
  geom_density() +
  labs(x = "Fraction irrigation lost", y = "Density", color = "Datasets \nagreeing") +
  theme_AP() +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  scale_color_viridis_d(option = "mako", direction = -1)

plot_curves_country


## ----merge_k, dependson="fraction_lost_k", fig.height=2.8, fig.width=5.5------------------------

# MERGE ########################################################################

left <- plot_grid(plot_global_k, plot_curves_country, ncol = 1, labels = "auto",
                  rel_heights = c(0.6, 0.4))
plot_grid(left, plot_collapse_profiles, ncol = 2, labels = c("", "c"),
          rel_widths = c(0.35, 0.65))


## ----plot_agreement, dependson="fraction_lost_k", fig.height=1.5, fig.width=5.5-----------------

# PLOTS SHOWING K-OF-10 AGREEMENT ##############################################

# Prepare data -----------------------------------------------------------------

dt_plot <- dt_pres[, .N, .(resolution, n_pos)]
dt_plot[, frac:= N / sum(N), by = resolution]

# Fraction of cells in the consensus core (>=9/10) -----------------------------

dt_plot[order(n_pos)] %>%
  .[n_pos >= 9] %>%
  .[, sum(frac), resolution]

# Fraction of cells irrigated by less than half of the datasets .---------------

dt_plot[order(n_pos)] %>%
  .[n_pos < 5 & n_pos != 0] %>%
  .[, sum(frac), resolution]

# Compute statistics -----------------------------------------------------------

dt_summary <- dt_plot[, .(total = sum(N),
                          absence = sum(N[n_pos == 0]),
                          contested = sum(N[n_pos >=1 & n_pos <=8]),
                          fewer_than_half = sum(N[n_pos < 5 & n_pos > 0]),
                          presence = sum(N[n_pos >=9])),
                      resolution]

dt_summary[, `:=`(frac_absence = absence / total,
                  frac_contested = contested / total,
                  frac_fewer_than_half = fewer_than_half / total,
                  frac_presence = presence / total)]

dt_summary

# Probability that a dataset detects irrigation in a given cell ----------------

dt_k <- dt_pres[, .N, .(resolution, n_pos)]

dt_prob <- dt_k[, .(mean_k = weighted.mean(n_pos, N),
                    detection_prob = weighted.mean(n_pos, N) / 10),
                resolution]

dt_prob

# In how many cells do fewer than half of the datasets detect irrigation? ------

dt_weak <- dt_k[, .(total_cells = sum(N),
  no_irrigation = sum(N[n_pos == 0]),
  any_irrigation = sum(N[n_pos >= 1]),
  weak_detection = sum(N[n_pos >= 1 & n_pos <= 5])), resolution]

dt_weak[, `:=`(frac_weak_all = weak_detection / total_cells,
               frac_weak_given_any = weak_detection / any_irrigation)]

dt_weak

# Plot distribution ------------------------------------------------------------

p1 <- plot_agreement <- ggplot(dt_plot, aes(n_pos, frac)) +
  geom_col() +
  facet_wrap(~resolution) +
  labs(x = "Datasets agreeing",
       y = "Fraction of cells") +
  theme_AP()

p1

# Prepare data for cumulative plot ---------------------------------------------

dt_ecdf <- dt_pres[, .N, .(resolution, n_pos)]
dt_ecdf[, frac:= N / sum(N), resolution]
setorder(dt_ecdf, resolution, n_pos)
dt_ecdf[, cdf:= cumsum(frac), resolution]

dt_ecdf[, resolution:= factor(resolution, levels = c("0.2deg","0.4deg","1deg"))]

# Plot -------------------------------------------------------------------------

p2 <- ggplot(dt_ecdf, aes(n_pos, cdf, colour = resolution)) +
  geom_step(linewidth = 0.8) +
  scale_x_continuous(breaks = 0:10) +
  scale_color_manual(values = res_palette, name = "") +
  labs(x = "Datasets agreeing", y = "Cum. grid cells",
       colour = "Resolution") +
  theme_AP() +
  theme(legend.position = c(0.7, 0.3))

p2

plot_grid(p1, p2, ncol = 2, rel_widths = c(0.7, 0.3), labels = "auto")



## ----plot_agreement2, dependson="plot_agreement", fig.height=2, fig.width=2---------------------

# Plot: consensus regime plot --------------------------------------------------

dt_regime <- dt_k[, .(absence = sum(N[n_pos == 0]),
                      contested = sum(N[n_pos >= 1 & n_pos <= 8]),
                      presence = sum(N[n_pos >= 9])), resolution]

dt_regime <- melt(dt_regime, id.vars = "resolution", variable.name = "regime",
                  value.name = "N")

dt_regime[, frac:= N / sum(N), by = resolution]

plot_regime <- ggplot(dt_regime, aes(resolution, frac, fill = regime)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(absence = "#4D4D4D",
                               contested = "#E69F00",
                               presence = "#0072B2"),
                    guide = guide_legend(nrow = 3),
                    name = "Agreement \nregime") +
  labs(x = NULL, y = "Fraction of cells") +
  theme_AP() +
  theme(legend.position = "top")


plot_regime


## ----crop_consequences, dependson="fraction_lost_k", fig.height=3, fig.width=4------------------

# GLOBAL CROP PRODUCTION UDNER DIFFERENT AGREEMENT RULES #######################
################################################################################

# Read crop datasets -----------------------------------------------------------

gaez_crops_02_aligned <- fread("./datasets/crops/gaez2015_crops_irrigated_0p2deg_aligned.csv")
gaez_crops_04_aligned <- fread("./datasets/crops/gaez2015_crops_irrigated_0p4deg_aligned.csv")
gaez_crops_1_aligned <- fread("./datasets/crops/gaez2015_crops_irrigated_1deg_aligned.csv")

# Extract crop cols-------------------------------------------------------------

crop_cols <- setdiff(names(gaez_crops_02_aligned),
                     c("lon", "lat", "country", "code", "continent"))

# Run for all grids ------------------------------------------------------------

res_list <- list("0.2deg" = gaez_crops_02_aligned,
                 "0.4deg" = gaez_crops_04_aligned,
                 "1deg"   = gaez_crops_1_aligned)

results <- lapply(names(res_list), function(r) {
  run_prod_analysis(r, res_list[[r]], dt_pres, crop_cols)
})

crop_all  <- rbindlist(lapply(results, `[[`, "crop"))
total_all <- rbindlist(lapply(results, `[[`, "total"))

total_all

# PLOT #########################################################################

# Heatmap ----------------------------------------------------------------------

plot_dt <- crop_all[!is.na(frac_vs_5)]

# Fix merged crop names only ---------------------------------------------------

plot_dt[, crop:= gsub("potatoandsweetpotato", "potato & sweet potato", crop)]
plot_dt[, crop:= gsub("yamsandotherroots", "yams & other roots", crop)]
plot_dt[, crop:= gsub("oilpalmfruit", "oil palm fruit", crop)]
plot_dt[, crop:= gsub("sugarbeet", "sugar beet", crop)]
plot_dt[, crop:= gsub("sugarcane", "sugar cane", crop)]
plot_dt[, crop:= gsub("othercereals", "other cereals", crop)]
plot_dt[, crop:= gsub("foddercrops", "fodder crops", crop)]
plot_dt[, crop:= gsub("cropsnes", "crops n.e.s.", crop)]

plot_dt[, k:= factor(k, levels = 10:5)]

plot_heatmap <- ggplot(plot_dt, aes(x = k, y = crop, fill = loss_vs_5)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = c("green", "yellow", "orange", "red"),
                       values = c(0, 0.5, 1),
                       breaks = c(0, 0.5, 1),
                       name    = "Fraction production \nlost") +
  facet_wrap(~resolution, ncol = 3) +
  labs(x = "Datasets agreeing",
       y = NULL) +
  theme_AP() +
  theme(axis.text.y = element_text(size = 7),
        legend.position = "top")

plot_heatmap


## ----global_production, dependson="crop_consequences", fig.height=1.8, fig.width=2--------------

# Global production ------------------------------------------------------------

plot_crop_total <- ggplot(total_all,
       aes(x = k, y = loss_vs_5, group = resolution, color = resolution)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = res_palette) +
  scale_x_reverse(breaks = 10:5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Datasets agreeing", y = "Fraction production lost",
       color = "Resolution") +
  theme_AP() +
  theme(legend.position = c(0.7, 0.7))

plot_crop_total


## ----global_breadbaskets, dependson="crop_consequences", fig.height=3, fig.width=4--------------

# COMPUTE BREADBASKETS #########################################################
################################################################################

# LAND SEA MASK -----------------------------------------------------------------

grain_crops <- intersect(c("wheat", "rice", "maize", "barley", "millet",
                           "sorghum", "othercereals"), names(gaez_crops_02_aligned))

lon_min <- min(gaez_crops_02_aligned$lon)
lon_max <- max(gaez_crops_02_aligned$lon)
lat_min <- min(gaez_crops_02_aligned$lat)
lat_max <- max(gaez_crops_02_aligned$lat)

template <- rast(xmin = lon_min, xmax = lon_max,
                 ymin = lat_min, ymax = lat_max,
                 resolution = 0.2,
                 crs = "EPSG:4326")

world <- vect(spData::world)
land_mask <- rasterize(world, template, field = 1, background = NA)

# Focus regions ----------------------------------------------------------------

# PREVIOUS LAND MASK EXCLUDES WATER BODIES IN THE AREAS ########################

# India: Indus Valley (Pakistan) and the Ganges Delta (Bangladesh).
# China: Hebei-Henan-Shandong axis
# USA: Ogallala Box; the Nebraska-to-Texas irrigation corridor.
# Egypt: Egyptian Nile
# Mediterranean: Spain, Italy, Greece, Turkey and North Africa

regions <- list("Indo-Gangetic Plain" = ext(68, 91, 21.5, 33),
                "North China Plain" = ext(113, 121, 34, 41),
                "US High Plains" = ext(-104, -96, 32, 44),
                "Nile Delta/Basin" = ext(29.5, 33, 24, 31.5),
                "Mekong Delta" = ext(104.5, 106.5, 8.5, 11.5),
                "Mediterranean Basin" = ext(-6, 36, 30, 46))

# Run for 0.2, 0.4, 1-----------------------------------------------------------

res_list <- list("0.2deg" = gaez_crops_02_aligned,
                 "0.4deg" = gaez_crops_04_aligned,
                 "1deg"   = gaez_crops_1_aligned)

results_regions_all <- rbindlist(
  lapply(names(res_list), function(r) {
    compute_breadbaskets(
      res_tag = r,
      crops_aligned = res_list[[r]],
      dt_pres = dt_pres,
      regions = regions,
      grain_crops = grain_crops,
      land_mask = land_mask
    )
  }),
  use.names = TRUE,
  fill = TRUE
)

# PLOT: FRACTION OF GRAIN PRODUCTION LOST ######################################

plot_fraction_lost <- ggplot(results_regions_all, aes(k, loss_vs_5,
                                                      color = resolution,
                                                      group = resolution)) +
  geom_line() +
  geom_point() +
  scale_x_reverse(breaks = 10:5) +
  scale_color_manual(values = res_palette) +
  labs(x = "Datasets agreeing",
       y = "Fraction production lost",
       color = "Resolution") +
  facet_wrap(~region, ncol = 3) +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  theme_AP() +
  theme(legend.position = "top")

plot_fraction_lost


## ----merge_crop_breadbasket, dependson=c("global_breadbaskets", "crop_consequences"), fig.height=5.3, fig.width=3.7----

# MERGE ########################################################################

plot_grid(plot_heatmap, plot_fraction_lost, ncol = 1, labels = "auto", rel_heights = c(0.55, 0.45))


## ----merge_crop_breadbasket2, dependson=c("global_breadbaskets", "crop_consequences"), fig.height=3, fig.width=5.7----

# MERGE ########################################################################

plot_fraction_lost <- plot_fraction_lost +
  facet_wrap(~region, ncol = 2)

top <- plot_grid(plot_heatmap, plot_fraction_lost, ncol = 2, labels = "auto",
                 rel_widths =  c(0.55, 0.45))

top


## ----irrigation_ET, dependson="fraction_lost_k", fig.height=2, fig.width=2----------------------

# LOAD THE HARMONIZED GRIDS ####################################################

tmpl_02 <- rast("./mha_stack_ll_02_aligned_NEAR.tif")[[1]]
tmpl_04 <- rast("./mha_stack_ll_04_aligned_NEAR.tif")[[1]]
tmpl_1 <- rast("./mha_stack_ll_10_aligned_NEAR.tif")[[1]]

templates <- list("0.2deg" = tmpl_02, "0.4deg" = tmpl_04, "1deg" = tmpl_1)

# Enforce lon/lat WGS84 on templates if needed----------------------------------

templates <- lapply(templates, function(t) {
  if (is.na(crs(t)) || crs(t) == "") crs(t) <- "EPSG:4326"
  t
})

# LOAD THE FILES ###############################################################
# Files come from the Yao et al 2025 paper (Nature Water)

# Year window ------------------------------------------------------------------

win <- "1985_2014_timmean"

# Base directory ---------------------------------------------------------------

base_dir <- "./yao_et_al_2025_datasets/water_fluxes/CESM2"

# Paste ------------------------------------------------------------------------

irr_files <- paste0(base_dir, "/CESM2_IRR0", 1:3, "_1901_2014_QFLX_EVAP_TOT_yearmean_", win)
noi_files <- paste0(base_dir, "/CESM2_NOI0", 1:3, "_1901_2014_QFLX_EVAP_TOT_yearmean_", win)

# BUILD ET STACKS + ENSEMBLE MEAN PER RESOLUTION ################################

d_list <- Map(function(irr_f, noi_f) delta_pair_resampled(irr_f, noi_f, templates),
              irr_files, noi_files)

# d_list length 3; each element is list("0.2deg"=rast, "0.4deg"=rast, "1deg"=rast)
dET_mean_rasters <- lapply(names(templates), function(res_name) {

  d_rasters <- lapply(d_list, function(one_pair) one_pair[[res_name]])
  s <- terra::rast(d_rasters) # 3-layer stack
  names(s) <- paste0("dET_IRR0", 1:3, "_mmday")

  # ensemble mean: cell-by-cell mean across 3 members
  m <- terra::app(s, fun = mean, na.rm = TRUE)
  names(m) <- paste0("CESM2_dET_mean_", win, "_", res_name, "_mmday")

  list(stack = s, mean = m)
})
names(dET_mean_rasters) <- names(templates)

# EXPORT TO DATA TABLE #########################################################

dET_dt <- lapply(names(templates), function(res_name) {
  tmpl <- templates[[res_name]]
  out <- dET_mean_rasters[[res_name]]

  dt_mean <- rast_to_dt_with_cell(out$mean, tmpl, value_name = "dET_mean_mmday")
  dt_mean[, resolution:= res_name]

  dt_members <- rbindlist(lapply(1:terra::nlyr(out$stack), function(i) {
    rr <- out$stack[[i]]
    dt <- rast_to_dt_with_cell(rr, tmpl, value_name = "dET_mmday")
    dt[, member:= names(rr)] # dET_IRR01_mmday / dET_IRR02_mmday / dET_IRR03_mmday
    dt[, resolution:= res_name]
    dt[]
  }))

  list(mean = dt_mean, members = dt_members)
})
names(dET_dt) <- names(templates)

# COMPUTE K-OF-10 AGREEMENT WITH CELL JOINS ####################################

resolutions <- c("0.2deg", "0.4deg", "1deg")
ks <- 5:10

# Compute ----------------------------------------------------------------------

out_k_all <- rbindlist(lapply(resolutions, function(res) {

  tmpl <- templates[[res]]

  # irrigation presence-------------

  dt_irrig_raw <- copy(dt_pres)[resolution == res]
  dt_irrig <- pres_to_cell(dt_irrig_raw, tmpl)

  # add cell area------------------

  dt_area <- cell_area_dt(tmpl)
  setkey(dt_irrig, cell)
  setkey(dt_area,  cell)
  dt_irrig <- dt_irrig[dt_area]  # add area_m2

  # increase ET mean table ----------

  dt_et <- copy(dET_dt[[res]]$mean)
  valcol <- setdiff(names(dt_et), c("cell","lon","lat","resolution"))[1]

  # Join: we keep all ET cells and bring irrigation info when available:
  # to summarise ONLY cells that exist in the irrigation table:
  #   dt_join <- dt_et[dt_irrig, on="cell"]
  setkey(dt_et, cell)
  dt_join <- dt_irrig[dt_et]

  rbindlist(lapply(ks, function(k) {

    dt_sub <- dt_join[n_pos >= k]

    # restrict weights to cells with non-missing ET and non-missing area
    ok <- !is.na(dt_sub[[valcol]]) & !is.na(dt_sub$area_m2)
    w  <- dt_sub$area_m2[ok]
    x  <- dt_sub[[valcol]][ok]

    # area-weighted mean increase ET (mm/day)
    mean_aw <- if (length(x) == 0L) NA_real_ else sum(x * w) / sum(w)

    # area-weighted aggregate: mm/day * m^2 -> m^3/day
    # 1 mm = 0.001 m
    sum_m3day <- if (length(x) == 0L) NA_real_ else sum(x * w) * 1e-3

    data.table(resolution = res,
               k = k,
               n_cells_irrig_k = sum(dt_irrig$n_pos >= k, na.rm = TRUE),
               n_cells_join_k = nrow(dt_sub),
               n_cells_with_dET_k = sum(!is.na(dt_sub[[valcol]])),
               mean_dET_mmday = mean(dt_sub[[valcol]], na.rm = TRUE),
               sum_dET_mmday_unweighted = sum(dt_sub[[valcol]], na.rm = TRUE),
               mean_dET_mmday_areaweighted = mean_aw,
               sum_dET_m3day_areaweighted = sum_m3day,
               area_m2_with_dET_k = if (length(w) == 0L) 0 else sum(w)
    )
  }))
}))

out_k_all

# PLOT #########################################################################

# Condition ifesle -------------------------------------------------------------

total_col <- if ("sum_dET_m3day_areaweighted" %in% names(out_k_all)) {
  "sum_dET_m3day_areaweighted"
} else if ("sum_dET_mmday_unweighted" %in% names(out_k_all)) {
  "sum_dET_mmday_unweighted"
} else {
  stop("out_k_all has neither sum_dET_m3day_areaweighted nor sum_dET_mmday_unweighted")
}

# Convert m3/day to km3/day ----------------------------------------------------

out_plot <- copy(out_k_all)

if (total_col == "sum_dET_m3day_areaweighted") {
  out_plot[, total_km3day := get(total_col) / 1e9]
  total_ylab <- "Total irrigation-\ninduced ET (km³ d-¹)"

} else {

  out_plot[, total_km3day := get(total_col)]
  total_ylab <- "Total irrigation-induced \nET signal (sum of mm/day)"
}

# Keep only k = 5, ..., 10------------------------------------------------------

out_plot <- out_plot[k %in% 5:10]

# Plot mean increase in irrigation-induced ET ----------------------------------
# (average irrigation-induced evapotranspiration change per grid cell)

plot_mean_ET <- ggplot(out_plot, aes(k, mean_dET_mmday, colour = resolution,
                                     group = resolution)) +
  geom_line() +
  geom_point() +
  labs(x = "Datasets agreeing",
       y = expression(atop("Mean " * Delta * ET[irr], "(mm d"^{-1} * ")"))) +
  scale_x_reverse(breaks = 10:5) +
  scale_color_manual(values = res_palette) +
  theme_AP() +
  theme(legend.position = "none")

# Plot total irrigation induced ET ---------------------------------------------
# (global irrigation-induced evapotranspiration signal, obtained by integrating
# over all retained cells; e.g., total amount of water transferred to the
# atmosphere per day because of irrigation)

plot_total_ET <- ggplot(out_plot, aes(k, total_km3day, colour = resolution, group = resolution)) +
  geom_line() +
  geom_point() +
  labs(x = "Datasets agreeing", y = total_ylab) +
  scale_x_reverse(breaks = 10:5) +
  scale_color_manual(values = res_palette) +
  theme_AP() +
  theme(legend.position = "none")

# Plot area with irrigation-induced ET retained under k-------------------------

out_plot <- copy(out_k_all)[k %in% 5:10]

# convert m2 to Mha -----------------------------------
out_plot[, area_Mha:= area_m2_with_dET_k / 1e10]

# also show contraction relative to k=5 within each resolution
out_plot[, area_frac_vs_k5:= area_m2_with_dET_k / area_m2_with_dET_k[k == 5], by = resolution]

# Total land area of grid cells with reported irrigated areas
plot_area <- ggplot(out_plot, aes(k, area_Mha, colour = resolution, group = resolution)) +
  geom_line() +
  geom_point() +
  labs(x = "Datasets agreeing", y = "Area irrigation \nsignal (Mha)") +
  scale_x_reverse(breaks = 10:5) +
  scale_color_manual(values = res_palette) +
  theme_AP() +
  theme(legend.position = "none")

plot_area

# Also in case: fraction of area retained vs k=5
plot_area_fraction <- ggplot(out_plot, aes(k, area_frac_vs_k5, colour = resolution, group = resolution)) +
  geom_hline(yintercept = 1, linewidth = 0.4) +
  geom_line() +
  geom_point() +
  labs(x = "Datasets agreeing",
       y = expression("Area retaining "~Delta*ET~" (fraction of k=5)")) +
  scale_x_reverse(breaks = 10:5) +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_color_manual(values = res_palette) +
  theme_AP() +
  theme(legend.position = "none")

plot_area_fraction


## ----merge_ET_plots, dependson="irrigation_ET", fig.height=2, fig.width=4-----------------------

bottom <- plot_grid(plot_mean_ET, plot_total_ET, plot_area, labels = c("c", "d", "e"),
                    ncol = 3)
bottom


## ----merge_ET_plots2, dependson=c("merge_ET_plots", "merge_crop_breadbasket2"), fig.height=4.3, fig.width=5.5----

plot_grid(top, bottom, ncol = 1, rel_heights = c(0.685, 0.325))


## ----final_figure1, dependson=c("plot_uasa", "define_plots_taumax", "plot_map"), fig.width=5.5----

left <- plot_grid(plot_uncertainty1, plot_agreement1,
                  plot_stacked +
                    labs(y = "Fraction of cells"), ncol = 1, labels = "auto")
plot_grid(left, plot_raster, ncol = 2, rel_widths = c(0.3, 0.7), labels = c("", "d"))


## ----final_figure1_1, dependson=c("plot_uasa", "define_plots_taumax", "plot_map"), fig.width=6----

left <- plot_grid(plot_uncertainty1, plot_agreement1,
                  plot_regime, ncol = 1, labels = "auto")
middle <- plot_grid(p1 + facet_wrap(~resolution, ncol = 1),
                    plot_stacked +
                      labs(y = "Fraction of cells"),
                    ncol = 1, rel_heights = c(0.7, 0.3), labels = c("d", "e"))
plot_grid(left, middle, plot_raster, ncol = 3, rel_widths = c(0.26, 0.26, 0.48), labels = c("", "", "f"))


## ----final_figure2, dependson= c("plot_uasa", "uasa_on_mapping_paradigm", "uasa_on_weights"), fig.height=1.5, fig.width=5.5----

top <- plot_grid(plot_indices1 +
                   labs(x = "", y = "Fraction \n variance") +
                   scale_fill_discrete(name = "Sobol' \n indices") +
                   theme(legend.position = c(0.4, 0.7)),
                 a1 + labs(x = "exclusion", y = "Fraction \n disagree"),
                 b1, c1, ncol = 4, labels = c("a", "", "", ""),
                 rel_widths = c(0.23, 0.37, 0.2, 0.2))
top

middle <- plot_grid(plot_indices2 +
                   labs(x = "", y = "Fraction \n variance") +
                   scale_fill_discrete(name = "") +
                   theme(legend.position = "none"),
                 a2 + labs(x = "exclusion", y = "Fraction \n disagree"),
                 b2, c2, ncol = 4, labels = c("b", "", "", ""),
                 rel_widths = c(0.23, 0.37, 0.2, 0.2))

middle

bottom <- plot_grid(plot_indices3 +
                   labs(x = "", y = "Fraction \n variance") +
                   scale_fill_discrete(name = "") +
                   theme(legend.position = "none"),
                 a3 + labs(x = "exclusion", y = "Fraction \n disagree"),
                 b3, c3, d3, ncol = 5, labels = c("c", "", "", "", ""),
                 rel_widths = c(0.22, 0.34, 0.16, 0.15, 0.15))

bottom


## ----merge_all_sa, dependson="final_figure2"----------------------------------------------------

plot_grid(top, middle, bottom, ncol = 1)


## ----session_information------------------------------------------------------------------------

# SESSION INFORMATION ##########################################################

sessionInfo()

## Return the machine CPU -----------------------------------------------------

cat("Machine:     "); print(get_cpu()$model_name)

## Return number of true cores -------------------------------------------------

cat("Num cores:   "); print(parallel::detectCores(logical = FALSE))

## Return number of threads ---------------------------------------------------

cat("Num threads: "); print(parallel::detectCores(logical = FALSE))


















library(data.table)
library(ggplot2)
library(scales)

#--------------------------------------------
# 1. Settings
#--------------------------------------------

datasets <- c(
  "gaez_v4", "giam", "gmia", "gripc", "luh2",
  "meier", "mirca2000", "mirca_os", "nagaraj", "spam"
)

countries_sel <- c("Vietnam", "Bangladesh", "Thailand")
resolution_sel <- "0.2deg"

#--------------------------------------------
# 2. Keep only required rows
#--------------------------------------------

# dt_wide appears already to be at 0.2 deg
dtw <- copy(dt_wide)[country %in% countries_sel]

# dt_pres has multiple resolutions, so keep only the one of interest
dtp <- copy(dt_pres)[country %in% countries_sel & resolution == resolution_sel]

#--------------------------------------------
# 3. Melt to long format
#--------------------------------------------

dtw_long <- melt(
  dtw,
  id.vars = c("lon", "lat", "country", "code", "continent"),
  measure.vars = datasets,
  variable.name = "dataset",
  value.name = "irr_area"
)

dtp_long <- melt(
  dtp,
  id.vars = c("lon", "lat", "country", "code", "continent", "resolution"),
  measure.vars = datasets,
  variable.name = "dataset",
  value.name = "pres"
)

#--------------------------------------------
# 4. Join magnitude with presence mask
#--------------------------------------------

plot_dt <- merge(
  dtw_long,
  dtp_long[, .(lon, lat, country, dataset, pres)],
  by = c("lon", "lat", "country", "dataset"),
  all.x = TRUE
)

# If for some reason pres is missing after join, treat as absence
plot_dt[is.na(pres), pres := 0]

# Force cells with no irrigation presence to white (= 0)
plot_dt[, irr_plot := fifelse(pres == 1, irr_area, 0)]

#--------------------------------------------
# 5. Optional: nicer labels/order
#--------------------------------------------

plot_dt[, dataset := factor(
  dataset,
  levels = datasets,
  labels = c(
    "GAEZ v4", "GIAM", "GMIA", "GRIPC", "LUH2",
    "Meier", "MIRCA2000", "MIRCA-OS", "Nagaraj", "SPAM"
  )
)]

plot_dt[, country := factor(
  country,
  levels = c("Vietnam", "Bangladesh", "Thailand")
)]

#--------------------------------------------
# 6. Plot
#--------------------------------------------

max_val <- plot_dt[, max(irr_plot, na.rm = TRUE)]

dt_wide[, .(country, gripc, gmia)] %>%
  .[country %in% c("Vietnam")]

p <- ggplot(plot_dt, aes(x = lon, y = lat, fill = irr_plot)) +
  geom_raster() +
  coord_equal() +
  facet_grid(country~dataset) +
  scale_fill_gradient(
    low = "red",
    high = "darkgreen",
    limits = c(0, max_val),
    trans = "sqrt",
    name = "Irrigated area"
  )+
  labs(
    x = NULL,
    y = NULL,
    title = "Identification of irrigated areas across datasets",
    subtitle = "Red = no irrigation detected; darker green = larger irrigated area"
  ) +
  theme_AP() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(face = "bold")
  )

p

p


### STUDY OF THE DEFINITIONS ###################################################
################################################################################

# TEST 1: STUDY OF DEFINITION HETEROGENEITY ####################################

dt_definitions <- fread("./llm_irrigated_definitions/outputs/irrigation_definition_table.csv")

# Create crosswalk table -------------------------------------------------------

crosswalk <- data.table(
  source_file = c(
    "Grogan et al. - 2022 - Global gridded crop harvested area, production, yield, and monthly physical area data circa 2015.pdf",
    "Hurtt et al. - 2020 - Harmonization of global land use change and management for the period 850–2100 (LUH2) for CMIP6.pdf",
    "Meier et al_2018_A global approach to estimate irrigated areas.pdf",
    "Nagaraj et al_2021_A new dataset of global irrigation areas from 2001 to 2015.pdf",
    "Portmann et al. - 2010 - MIRCA2000-Global monthly irrigated and rainfed cro.pdf",
    "Salmon et al. - 2015 - Global rain-fed, irrigated, and paddy croplands A new high resolution map derived from remote sensi.pdf",
    "Siebert et al_2013_Update of the digital global map of irrigation areas to version 5.pdf",
    "Thenkabail et al_2009_Global irrigated area map (GIAM), derived from remote sensing, for the end of.pdf",
    "Yu et al. - 2020 - A cultivated planet in 2010 – Part 2 The global gridded agricultural-production maps.pdf",
    "s41597-024-04313-w-1.pdf"
  ),
  dataset = c(
    "gaez_v4",
    "luh2",
    "meier",
    "nagaraj",
    "mirca2000",
    "gripc",
    "gmia",
    "giam",
    "spam",
    "mirca_os"
  )
)

# Keep the columns needed for the audit ----------------------------------------

vars_def <- c("source_file",
              "paper_title",
              "irrigation_target",
              "representation_unit",
              "temporal_concept",
              "paddy_treatment",
              "permanent_crops_treatment",
              "irrigation_basis")

dt_definitions <- dt_definitions[, ..vars_def]

# Merge with dataset column ----------------------------------------------------

dt_definitions <- merge(dt_definitions, crosswalk, by = "source_file", all.x = TRUE)

# Convert empty strings to NA, then label as "unclear" ------------------------

for (v in setdiff(vars_def, c("source_file", "paper_title"))) {
  dt_definitions[get(v) == "" | is.na(get(v)), (v):= "unclear"]
}

# Make frequency tables --------------------------------------------------------

freq_list <- lapply(
  setdiff(vars_def, c("source_file", "paper_title")),
  function(v) make_freq_table(dt_definitions, v)
)

freq_dt <- rbindlist(freq_list, use.names = TRUE, fill = TRUE)
freq_dt

# Share of mixed / unclear datasets---------------------------------------------

summary_mixed_unclear <- rbindlist(lapply(
  setdiff(vars_def, c("source_file", "paper_title")),
  function(v) {
    tmp <- dt_definitions[, .(n_total = .N,
                              n_unclear = sum(get(v) == "unclear", na.rm = TRUE),
                              n_mixed = sum(get(v) %in% c("mixed", "mixed_or_multiple"), na.rm = TRUE))]
    tmp[, variable:= v]
    tmp[, share_unclear:= n_unclear / n_total]
    tmp[, share_mixed:= n_mixed / n_total]
    tmp[]
  }), fill = TRUE)

print(summary_mixed_unclear)


# Barplots of definitional heterogeneity ---------------------------------------

p1 <- plot_bar_var(freq_dt, "irrigation_target")
p2 <- plot_bar_var(freq_dt, "representation_unit")
p3 <- plot_bar_var(freq_dt, "temporal_concept")
p4 <- plot_bar_var(freq_dt, "paddy_treatment")
p5 <- plot_bar_var(freq_dt, "permanent_crops_treatment")
p6 <- plot_bar_var(freq_dt, "irrigation_basis")

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3, labels = "auto")

# Dataset × definition heatmap panel

dt_long <- melt(dt_definitions, id.vars = c("dataset"),
                measure.vars = c("irrigation_target",
                                 "representation_unit",
                                 "temporal_concept",
                                 "permanent_crops_treatment",
                                 "irrigation_basis"),
  variable.name = "definition_variable",
  value.name = "definition_value")

dataset_order <- dt_long[, .(
  score = sum(definition_value %in% c("unclear", "mixed", "mixed_or_multiple"))),
  by = dataset][order(-score, dataset)]$dataset

dt_long[, dataset:= factor(dataset, levels = rev(dataset_order))]

dt_long[, class_type:= fifelse(
  definition_value == "unclear", "unclear",
  fifelse(definition_value %in% c("mixed", "mixed_or_multiple"), "mixed", "specified")
)]

# Plot heatmap: Cells report the definition assigned to each dataset;
# colour highlights whether the definition is specified, mixed, or unclear

p_heat <- ggplot(dt_long, aes(definition_variable, dataset, fill = class_type)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = definition_value), size = 2.8) +
  scale_fill_manual(values = c("specified" = "#4daf4a", "mixed" = "#ff7f00",
                               "unclear" = "#bdbdbd")) +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "Definition heterogeneity audit across irrigation datasets",
       subtitle = "") +
  theme_AP() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top")

p_heat

# Summary metrics---------------------------------------------------------------

n_datasets <- uniqueN(dt_definitions$dataset)

summary_text_dt <- rbindlist(lapply(
  c("irrigation_target",
    "representation_unit",
    "temporal_concept",
    "paddy_treatment",
    "permanent_crops_treatment",
    "irrigation_basis"),
  function(v) {
    tmp <- dt_definitions[, .(
      n_unique_categories = uniqueN(get(v)),
      n_unclear = sum(get(v) == "unclear"),
      n_mixed = sum(get(v) %in% c("mixed", "mixed_or_multiple"))
    )]
    tmp[, variable:= v]
    tmp[, share_unclear:= round(n_unclear / n_datasets, 3)]
    tmp[, share_mixed:= round(n_mixed / n_datasets, 3)]
    tmp[]
  }
), fill = TRUE)

print(summary_text_dt)

# TEXT FOR THE PAPER:

# The literature audit reveals that global irrigation maps do not represent a single
# concept of irrigated area. Across products, irrigation is variously defined as
# equipped, actual, or harvested irrigation, while treatment of paddy rice, temporal
# reference, and irrigation basis also varies. Several datasets are mixed or
# conceptually unclear, indicating that definitional heterogeneity is real and must
# be explicitly accounted for before assessing whether it can explain the observed
# detectability crisis.

# TEST 2: PAIRWISE SAME DEFINITION VERSUS DIFFERENT DEFINITION TEST ------------

# choose one resolution for the pairwise test-----------------------------------

resolution_sel <- "0.2deg"

# PREPARE PRESENCE DATA --------------------------------------------------------

dtp <- copy(dt_pres)[resolution == resolution_sel]

# Dataset columns available in dt_pres -----------------------------------------

id_cols <- c("lon", "lat", "country", "code", "continent", "resolution", "n_pos")
dataset_cols <- setdiff(names(dtp), id_cols)

# Keep only datasets that are both in def_dt and dt_pres ------------------------

common_datasets <- intersect(dt_definitions$dataset, dataset_cols)

dt_definitions <- dt_definitions[dataset %in% common_datasets]
dtp <- dtp[, c(setdiff(id_cols, "n_pos"), common_datasets), with = FALSE]

cat("\nDatasets used in pairwise test:\n")
print(common_datasets)



pair_dt <- compute_pairwise_disagreement(dtp, common_datasets)


# ATTACH DEFINITIONAL ATTRIBUTES TO EACH PAIR ----------------------------------

# Attributes for dataset i -----------------------------------------------------

def_i <- copy(dt_definitions)
setnames(def_i,
         old = c("dataset", "irrigation_target", "representation_unit",
                 "temporal_concept", "paddy_treatment",
                 "permanent_crops_treatment", "irrigation_basis"),
         new = c("dataset_i", "irrigation_target_i", "representation_unit_i",
                 "temporal_concept_i", "paddy_treatment_i",
                 "permanent_crops_treatment_i", "irrigation_basis_i"))

# Attributes for dataset j------------------------------------------------------

def_j <- copy(dt_definitions)
setnames(def_j,
         old = c("dataset", "irrigation_target", "representation_unit",
                 "temporal_concept", "paddy_treatment",
                 "permanent_crops_treatment", "irrigation_basis"),
         new = c("dataset_j", "irrigation_target_j", "representation_unit_j",
                 "temporal_concept_j", "paddy_treatment_j",
                 "permanent_crops_treatment_j", "irrigation_basis_j"))

pair_dt <- merge(pair_dt, def_i[, .(
  dataset_i,
  irrigation_target_i,
  representation_unit_i,
  temporal_concept_i,
  paddy_treatment_i,
  permanent_crops_treatment_i,
  irrigation_basis_i
)], by = "dataset_i", all.x = TRUE)

pair_dt <- merge(pair_dt, def_j[, .(
  dataset_j,
  irrigation_target_j,
  representation_unit_j,
  temporal_concept_j,
  paddy_treatment_j,
  permanent_crops_treatment_j,
  irrigation_basis_j
)], by = "dataset_j", all.x = TRUE)

# same-definition indicators----------------------------------------------------

pair_dt[, same_irrigation_target:= irrigation_target_i == irrigation_target_j]
pair_dt[, same_representation_unit:= representation_unit_i == representation_unit_j]
pair_dt[, same_temporal_concept:= temporal_concept_i == temporal_concept_j]
pair_dt[, same_paddy_treatment:= paddy_treatment_i == paddy_treatment_j]
pair_dt[, same_permanent_crops_treatment:=permanent_crops_treatment_i == permanent_crops_treatment_j]
pair_dt[, same_irrigation_basis:= irrigation_basis_i == irrigation_basis_j]

# stricter indicators ----------------------------------------------------------

pair_dt[, same_core_definition:= same_irrigation_target & same_representation_unit &
          same_temporal_concept]

pair_dt[, same_extended_definition:= same_irrigation_target & same_representation_unit &
          same_temporal_concept & same_paddy_treatment & same_irrigation_basis]

# optional: exclude unclear/mixed comparisons-----------------------------------

pair_dt[, clear_irrigation_target:=
          !(irrigation_target_i %in% c("unclear", "mixed_or_multiple")) &
          !(irrigation_target_j %in% c("unclear", "mixed_or_multiple"))]

pair_dt[, clear_core_definition :=
          !(irrigation_target_i %in% c("unclear", "mixed_or_multiple")) &
          !(irrigation_target_j %in% c("unclear", "mixed_or_multiple")) &
          !(representation_unit_i %in% c("unclear", "mixed")) &
          !(representation_unit_j %in% c("unclear", "mixed")) &
          !(temporal_concept_i %in% c("unclear", "mixed")) &
          !(temporal_concept_j %in% c("unclear", "mixed"))]


# 5. DESCRIPTIVE COMPARISONS ---------------------------------------------------

# Mean disagreement by same/different target
desc_target <- pair_dt[, .(
  n_pairs = .N,
  mean_disagree_all = mean(share_disagree_all, na.rm = TRUE),
  median_disagree_all = median(share_disagree_all, na.rm = TRUE),
  mean_disagree_union = mean(share_disagree_union, na.rm = TRUE),
  median_disagree_union = median(share_disagree_union, na.rm = TRUE),
  mean_jaccard_dissimilarity = mean(jaccard_dissimilarity, na.rm = TRUE)
), by = same_irrigation_target]

print(desc_target)

# Mean disagreement by same/different extended definition
desc_extended <- pair_dt[, .(
  n_pairs = .N,
  mean_disagree_all = mean(share_disagree_all, na.rm = TRUE),
  median_disagree_all = median(share_disagree_all, na.rm = TRUE),
  mean_disagree_union = mean(share_disagree_union, na.rm = TRUE),
  median_disagree_union = median(share_disagree_union, na.rm = TRUE),
  mean_jaccard_dissimilarity = mean(jaccard_dissimilarity, na.rm = TRUE)
), by = same_extended_definition]

print(desc_extended)


# WILCOXON TESTS ---------------------------------------------------------------

# same vs different irrigation target -----------------
# Are the disagreement values in the two groups statistically different?
# It compares two groups of pairs of datasets:

# Group 1: pairs where same_irrigation_target == TRUE
# (both datasets claim to represent the same type of irrigation (e.g. both “actual”))
#
# Group 2: pairs where same_irrigation_target == FALSE
# (datasets represent different concepts (e.g. “actual” vs “equipped”))

# What is being compared?
# The variable: share_disagree_union
# (the fraction of cells where the two datasets disagree,
# conditional on at least one detecting irrigation.)
#
# What the test asks: Are the disagreement values in the two groups statistically
# different?

w_target <- wilcox.test(share_disagree_union ~ same_irrigation_target,
                        data = pair_dt)

print(w_target)
# Dataset pairs representing the same irrigation concept do not exhibit lower
# disagreement than those representing different concepts (Wilcoxon test, $p = 0.97$).
# This indicates that definitional alignment does not reduce disagreement.

# same vs different core definition-------------------

pair_dt[, core_similarity_score:= same_irrigation_target + same_representation_unit +
          same_temporal_concept]

cor.test(pair_dt$core_similarity_score, pair_dt$share_disagree_union, method = "spearman")

# Increasing definitional similarity does not reduce disagreement
# Even when datasets match across multiple definitional dimensions, they do not
# agree more spatially.

# PARAGRAPH: Disagreement does not decline with increasing definitional similarity
# (Spearman’s $\rho = -0.04$, $p = 0.78$), indicating that even datasets aligned
# across multiple conceptual dimensions fail to identify irrigation consistently.


# REGRESSION MODEL--------------------------------------------------------------

# Simple OLS on pairwise disagreement
# Dependent variable: disagreement conditional on at least one map detecting irrigation
mod1 <- lm(share_disagree_union ~
             same_irrigation_target +
             same_representation_unit +
             same_temporal_concept +
             same_paddy_treatment +
             same_irrigation_basis,
           data = pair_dt)

summary(mod1)
#
# None of the definitional dimensions significantly explains pairwise
# disagreement (all $p > 0.18$), and together they account for less than 8% of
# its variance. Definitional differences therefore fail to explain the widespread
# inconsistencies across irrigation maps.

# Matching definitions across five dimensions does not systematically reduce disagreement.

# PLOTS ------------------------------------------------------------------------


# Boxplot: same vs different irrigation target
p_target <- ggplot(pair_dt,
  aes(x = factor(same_irrigation_target,
                 levels = c(FALSE, TRUE),
                 labels = c("Different target", "Same target")),
      y = share_disagree_union)
) +
  geom_boxplot() +
  labs(
    x = NULL,
    y = "Pairwise disagreement (conditional on union)",
    title = "Pairwise disagreement for same-definition vs different-definition pairs",
    subtitle = "Comparison by irrigation target"
  ) +
  theme_AP()

p_target

# Boxplot: same vs different extended definition
p_ext <- ggplot(
  pair_dt,
  aes(x = factor(same_extended_definition,
                 levels = c(FALSE, TRUE),
                 labels = c("Different extended definition", "Same extended definition")),
      y = share_disagree_union)
) +
  geom_boxplot() +
  labs(
    x = NULL,
    y = "Pairwise disagreement (conditional on union)",
    title = "Pairwise disagreement for same-definition vs different-definition pairs",
    subtitle = "Comparison by extended definition"
  ) +
  theme_AP()

p_ext


# Global irrigation datasets do not represent a single, consistent concept of irrigated area.
# Across products, irrigation is variously defined as actual, equipped or mixed irrigation,
# and differs in representation unit (presence versus area), temporal reference
# (annual, seasonal or reference-period snapshots), and treatment of paddy rice and
# permanent crops (Fig. X). Several datasets combine multiple concepts or do not
# explicitly state their definition, indicating substantial conceptual heterogeneity
# across the ensemble.
#
# If definitional differences were the primary source of disagreement, datasets
# sharing the same conceptual definition should exhibit higher spatial agreement.
# However, this is not the case. Pairwise comparisons show that disagreement does
# not differ between dataset pairs that share the same irrigation target and those
# that do not (Wilcoxon test, p = 0.97), indicating that datasets aligned along
# several conceptual attributes do not agree more on where irrigation exists.
#
# We further tested whether definitional harmonisation reduces disagreement at the
# ensemble level by restricting the analysis to subsets of datasets sharing similar
# definitions. Although disagreement appears lower within these subsets, this
# reduction is entirely explained by the smaller number of datasets involved. When
# compared against randomly sampled subsets of equal size, definition-consistent
# subsets do not exhibit lower disagreement (Fig. X), demonstrating that the
# apparent reduction is a statistical artefact rather than a consequence of
# conceptual alignment.
#
# A multivariate regression including all definitional dimensions confirms these
# results. None of the variables significantly explains pairwise disagreement (all p > 0.18)
# and together they account for less than 8% of its variance (adjusted R2<0). Thus,
# even when considered jointly, differences in irrigation definitions fail to explain
# the widespread inconsistencies across datasets.
#
# Taken together, these results show that the detectability crisis cannot be attributed
# to definitional heterogeneity. Datasets do not disagree because they represent
# different types of irrigation; they disagree because they do not identify the
# same locations as irrigated.








# TEST NUMBER THREE:############################################################
##################

# Clean
for (v in c("irrigation_target", "representation_unit", "temporal_concept")) {
  def_dt[get(v) == "" | is.na(get(v)), (v) := "unclear"]
}

# Merge with your dataset names (same crosswalk as Test 2)
# assume def_dt already has column "dataset"

#--------------------------------------------------
# Select subsets
#--------------------------------------------------

subset_target <- def_dt[
  irrigation_target != "unclear",
  .(datasets = list(dataset)),
  by = irrigation_target
][order(-lengths(datasets))]

subset_unit <- def_dt[
  representation_unit != "unclear",
  .(datasets = list(dataset)),
  by = representation_unit
][order(-lengths(datasets))]

subset_time <- def_dt[
  temporal_concept != "unclear",
  .(datasets = list(dataset)),
  by = temporal_concept
][order(-lengths(datasets))]

print(subset_target)
print(subset_unit)
print(subset_time)


dt <- dcast(dt, lon + lat + resolution + country + code + continent ~ dataset,
      value.var = "mha") %>%
  .[resolution == "0.2deg"]

compute_disagreement_fraction <- function(dt, dataset_names) {

  A <- as.matrix(dt[, ..dataset_names])

  pres <- A > 0

  n_present <- rowSums(pres, na.rm = TRUE)
  n_datasets <- length(dataset_names)

  disagreement <- n_present > 0 & n_present < n_datasets

  frac_disagree <- mean(disagreement, na.rm = TRUE)

  return(frac_disagree)
}



# Use same resolution as main analysis
dtp <- dt_pres[resolution == "0.2deg"]

# All datasets
all_datasets <- intersect(names(dtp), def_dt$dataset)

frac_full <- compute_disagreement_fraction(dtp, all_datasets)

#--------------------------------------------------
# Loop over subsets
#--------------------------------------------------

results <- list()

for (i in seq_len(nrow(subset_target))) {

  ds <- subset_target$datasets[[i]]
  name <- subset_target$irrigation_target[i]

  ds <- intersect(ds, names(dtp))

  if (length(ds) >= 3) {

    frac <- compute_disagreement_fraction(dtp, ds)

    results[[length(results) + 1]] <- data.table(
      grouping = "irrigation_target",
      category = name,
      n_datasets = length(ds),
      frac_disagree = frac
    )
  }
}

for (i in seq_len(nrow(subset_unit))) {

  ds <- subset_unit$datasets[[i]]
  name <- subset_unit$representation_unit[i]

  ds <- intersect(ds, names(dtp))

  if (length(ds) >= 3) {

    frac <- compute_disagreement_fraction(dtp, ds)

    results[[length(results) + 1]] <- data.table(
      grouping = "representation_unit",
      category = name,
      n_datasets = length(ds),
      frac_disagree = frac
    )
  }
}

for (i in seq_len(nrow(subset_time))) {

  ds <- subset_time$datasets[[i]]
  name <- subset_time$temporal_concept[i]

  ds <- intersect(ds, names(dtp))

  if (length(ds) >= 3) {

    frac <- compute_disagreement_fraction(dtp, ds)

    results[[length(results) + 1]] <- data.table(
      grouping = "temporal_concept",
      category = name,
      n_datasets = length(ds),
      frac_disagree = frac
    )
  }
}

res_dt <- rbindlist(results)

# Add full ensemble
res_dt <- rbind(
  data.table(
    grouping = "full",
    category = "all datasets",
    n_datasets = length(all_datasets),
    frac_disagree = frac_full
  ),
  res_dt
)

print(res_dt)

fwrite(res_dt, "outputs/test3_definition_restriction_results.csv")

