[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19001232.svg)](https://doi.org/10.5281/zenodo.19001232)

# Where irrigation exists is globally contested

[Arnald Puy](https://www.arnaldpuy.com/), Olivia Richards, Seth N. Linga, Samuel
Flinders, Carmen Aguiló-Rivera, Josh Larsen

This study analyzes the disagreement across datasets in the identification of
global irrigated areas.

## Abstract

*Irrigation is one of the most extensive human modifications of the land surface and a key driver of food production, water use and regional climate. However, the location of irrigated land has never been systematically verified at global scale. Here we show that where irrigation exists is globally contested. Across ten datasets, 60–90% of cropland grid cells disagree on whether irrigation exists at all, remaining at 40–60% when irrigation is required to exceed 1% of the cell to confirm presence. Enforcing agreement across datasets reduces the irrigated area of irrigation hotspots such as Vietnam, China or Thailand by 70%, the share of global crop production associated with irrigation by 60-90% and the global irrigation-induced evapotranspiration signal by 67%. The consequences extend to greenhouse gas accounting as irrigated area misclassifications can underestimate rice methane emissions by 58%. This disagreement persists regardless of which datasets, methods or resolutions are considered and is consistent with irrigated areas being a non-identifiable variable. The way forward requires acknowledging that whether an area is irrigated is a judgement call that data alone cannot resolve.*

## Maps

The irrigated area maps used on our study are the following

* [Meier et al](https://hess.copernicus.org/articles/22/1119/2018/)  - The map by Meier et al.
* [GIAM](https://www.tandfonline.com/doi/full/10.1080/01431160802698919)  - The Global Irrigated Area Map.
* [GMIA](https://www.fao.org/aquastat/en/geospatial-information/global-maps-irrigated-areas/latest-version)  - The Global Map of Irrigated Areas.
* [GRIPC](https://www.sciencedirect.com/science/article/abs/pii/S0303243415000240?via%3Dihub)  - The Global Rain-fed, Irrigated and Paddy Croplands map.
* [Nagaraj](https://www.sciencedirect.com/science/article/abs/pii/S0309170821000658?via%3Dihub)  - The map by Nagaraj et al.
* [MIRCA-2000](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2008GB003435)  - The Monthly Irrigated and Rainfed Crop Areas map.
* [GAEZ+2015](https://www.nature.com/articles/s41597-021-01115-2)  - The Global Agro-Ecological Zones version 3 map.
* [SPAM2010](https://essd.copernicus.org/articles/12/3545/2020/)  - The 2010 Spatial Production Allocation Model map.
* [MIRCA-OS](https://www.nature.com/articles/s41597-024-04313-w)  - The Monthly Irrigated and Rainfed Crop Areas Open Source map.
* [LUH2](https://luh.umd.edu/)  - The Land-Use Harmonization$^2$ map.

## Software requirements

The workflow is implemented entirely in **R**. The following packages are required:

| Package | Purpose |
|---|---|
| `data.table` | Fast tabular data manipulation |
| `terra` | Raster reading, reprojection and aggregation |
| `tidyverse` | Data wrangling and plotting utilities |
| `ggplot2` | Visualisation (loaded via `tidyverse`) |
| `cowplot` | Multi-panel figure composition |
| `scales` | Axis and colour scale helpers |
| `here` | Portable file paths |
| `sf` | Vector spatial data (country polygons) |
| `rnaturalearth` | Country boundary data |
| `countrycode` | Country name and code standardisation |
| `sp` / `rworldmap` | Legacy spatial utilities (original dataset processing) |
| `readxl` | Reading the dataset-survey Excel file |
| `ncdf4` | Reading NetCDF climate model output |
| `wesanderson` | Colour palettes |
| `benchmarkme` | System information for reproducibility reporting |
| `parallel` | Parallel computation |
| `sensobol` | Uncertainty and sensitivity analyses |
| `magrittr` | Pipe operator |

Install all packages at once with:

```r
install.packages(c(
  "data.table", "terra", "tidyverse", "cowplot", "scales", "here", "sf",
  "rnaturalearth", "countrycode", "sp", "rworldmap", "readxl", "ncdf4",
  "wesanderson", "benchmarkme", "parallel", "sensobol", "magrittr"
))
```

A reproducibility seed of `123` is set at the start of each script.

## Replication

We provide all the functions needed to replicate our workflow in the `functions` folder.

### Functions

The `functions` folder contains all the custom functions coded for the analysis.
They are all sourced automatically from the `.R`, `.pdf` and `.Rmd` files and therefore the
user of the code does not need to source them separately.

The table below describes each function file and its role in the workflow:

| File | Description |
|---|---|
| `add_country_continent.R` | Assigns country name and continent to each grid cell via spatial join with Natural Earth polygons; processes data in configurable chunks to manage memory. |
| `compute_country_totals.R` | Sums irrigated area (Mha) per country for every dataset and computes within-dataset country rank. |
| `compute_loss_for_k.R` | For a given agreement threshold *k*, calculates what fraction of each dataset's irrigated cells (and their area) falls outside the *k*-consensus core. |
| `compute_rank_correlations.R` | Computes pairwise Kendall and Spearman rank correlations of country-level irrigation totals across all dataset pairs. |
| `compute_tau_fun.R` | Core detection-threshold analysis: given a minimum area threshold τ (Mha), classifies each cell as absent/present/disagree and returns a summary of agreement fractions, presence classes (None / Minority / Majority / All), and a marginal vs existential decomposition of disagreement. |
| `compute_tau_max_disagreement_fun.R` | Finds the largest τ at which existential disagreement persists in a cell (i.e. at least one dataset reports zero irrigation while another exceeds τ). |
| `compute_tau_weighted_fun.R` | Variant of `compute_tau_fun` that allows datasets to be weighted by user-defined scenarios (e.g. by publication year or methodological family). |
| `compute_topN_overlap.R` | Computes pairwise Jaccard similarity and intersection size of the top-N irrigation countries across all dataset pairs. |
| `country_to_continent.R` | Maps ISO country codes to continent labels. |
| `earthstat_functions.R` | Utilities for reading and processing EarthStat crop-harvested-area rasters used in the crop-production impact analysis. |
| `filter_long_dt_by_crop_mask_aggregated.R` | Applies a cropland mask to a long-format data table at a target resolution; supports two rules: `"fraction"` (cell passes if cropland fraction exceeds a threshold) and `"any"` (cell passes if any cropland pixel is present). |
| `functions_evapotranspiration.R` | Reads CESM climate model output (irrigation vs. no-irrigation runs), resamples to the three analysis grids, and computes irrigation-induced evapotranspiration (dET = ET_IRR − ET_NOI) in mm/day. |
| `lonlat_to_country.R` | Point-in-polygon lookup that returns the country name for a vector of longitude/latitude coordinates. |
| `non_identifiability_tests.R` | Stability and convergence tests: (1) leave-one-out and multi-tau perturbation analysis to measure cell-level classification flip rates; (2) random and exact subset-size sweeps to test whether agreement improves as more datasets are added. |
| `plot_attr.R` | Generates side-by-side boxplots comparing pairwise disagreement rates between datasets that share vs. differ in a given attribute (e.g. resolution, method family), with Wilcoxon test annotations. |
| `regridding_functions.R` | Creates explicit lon/lat and equal-area (Mollweide) template rasters; rasterises point datasets onto templates; reprojects and aggregates each irrigated-area map to 0.2°, 0.4° and 1° grids via a two-step equal-area intermediate; QA helpers for threshold-stability checks. |
| `run_test.R` | Labels dataset pairs by whether they share a given metadata attribute, runs a Wilcoxon rank-sum test and a permutation test for pairwise disagreement, and returns structured output for `plot_attr`. |

### Code

We offer the code in `.R`, `.pdf` and `.Rmd`. The full workflow is divided into the following scripts:

#### 1. `code_original_datasets`
Ingests each raw raster in its native format and converts it to a standardised long-format CSV with columns `lon`, `lat`, `country`, `continent`, `mha` and `dataset`. Area is expressed in **millions of hectares (Mha)** per grid cell.

The analyst must download all rasters from their respective sources. Most can be obtained by clicking the links in the section [Maps](#maps) above. Datasets that are under embargo or whose links are no longer active must be requested directly from the original authors.

Dataset-specific notes:

* **Meier et al.** — 1 km grid; pixels are classified 1–4 (irrigated) vs 0 (non-irrigated); cell area is computed from the geodetic formula for latitude-dependent pixel size.
* **GIAM** — 30 arc-sec grid; class-based Irrigation Area Fractions (IAFs) are used to derive Total Area Actually Irrigated (TAAI) in Mha.
* **GMIA, GRIPC, Nagaraj, MIRCA-2000, GAEZ+2015, SPAM2010, MIRCA-OS** — all regridded from their native resolutions.
* **LUH2** — the `irrig` fraction layer is multiplied by cell area to obtain Mha.

#### 2. `code_harmonization`
Harmonises all datasets to three common resolutions using a **two-step equal-area workflow**:

1. Each dataset is rasterised onto a global Mollweide equal-area grid (20 km resolution).
2. The Mollweide raster is reprojected back to geographic coordinates at 0.2°, 0.4° and 1° using bilinear resampling.

Template rasters are created with `make_template_ll()` for the three lon/lat grids and a fixed Mollweide extent. Zero-only or all-NA cells are dropped after harmonisation.

Output: a single long-format file `datasets/irrigated_areas/irrigated_areas.csv` with one row per dataset × grid cell combination, containing harmonised Mha values and administrative metadata.

#### 3. `code_main_analysis`
Core detectability analysis at the grid-cell level. The key steps are:

1. **Threshold sweep** (`compute_tau_fun`) — for minimum-area thresholds τ ranging from 0 to the 99th percentile of non-zero cell values, each grid cell is classified as:
   - *None*: no dataset detects irrigation (n = 0).
   - *Minority presence*: fewer than half of datasets detect irrigation.
   - *Majority presence*: more than half but not all datasets detect irrigation.
   - *All*: every dataset detects irrigation.
   - *Disagreement*: at least one dataset disagrees (1 ≤ n < 10).

2. **Disagreement decomposition** — disagreement is further split into:
   - *Marginal*: all datasets record some irrigation but not all exceed τ (intensity disagreement).
   - *Existential*: at least one dataset records zero irrigation while others exceed τ (presence/absence disagreement).

3. **Country-level analysis** (`compute_country_totals`, `compute_rank_correlations`, `compute_topN_overlap`) — quantifies how enforcing agreement (increasing *k*) shifts national irrigation rankings and totals.

4. **Crop-production impact** (`filter_long_dt_by_crop_mask_aggregated`, `earthstat_functions`) — measures how much of global crop production is attributed to irrigation under different *k* thresholds.

5. **ET impact** (`functions_evapotranspiration`) — quantifies how much irrigation-induced evapotranspiration (CESM model output) is associated with cells excluded at each *k* threshold.

6. **Stability tests** (`non_identifiability_tests`) — leave-one-out and multi-τ perturbation analyses confirm that disagreement is not an artefact of any single dataset or threshold choice.

## Supplementary material

* `irrigated_area_policy_SM.pdf` — *Policies using irrigated area statistics*. A review of policy documents, programmes 
and reports from major institutions (FAO Sustainable Development Goals, FAO Hand-in-Hand, 
World Bank, Asian Development Bank, etc.) that rely on irrigated area datasets to 
inform decisions on land degradation, food security, water scarcity and agricultural 
investment. The document catalogues which datasets each initiative cites and 
highlights that policy decisions are typically based on a single 
dataset (most often GMIA).

## Repository structure

```
code_detectability_irrigated_areas/
├── README.md
├── code_original_datasets.{R,Rmd,pdf}    # Ingest raw rasters → long-format CSV
├── code_harmonization.{R,Rmd,pdf}         # Regrid all datasets to 0.2°/0.4°/1°
├── code_main_analysis.{R,Rmd,pdf}         # Core detectability analysis
├── irrigated_area_policy_SM.pdf            # Supplementary material: policies that rely on irrigated-area datasets
├── functions/                              # Custom R functions (sourced automatically)
│   ├── add_country_continent.R
│   ├── compute_country_totals.R
│   ├── compute_loss_for_k.R
│   ├── compute_rank_correlations.R
│   ├── compute_tau_fun.R
│   ├── compute_tau_max_disagreement_fun.R
│   ├── compute_tau_weighted_fun.R
│   ├── compute_topN_overlap.R
│   ├── country_to_continent.R
│   ├── earthstat_functions.R
│   ├── filter_long_dt_by_crop_mask_aggregated.R
│   ├── functions_evapotranspiration.R
│   ├── lonlat_to_country.R
│   ├── non_identifiability_tests.R
│   ├── plot_attr.R
│   ├── regridding_functions.R
│   └── run_test.R
└── irrigation_masks/                       # Consensus irrigated-area masks
    ├── irrigation_mask_02.csv              # 0.2° resolution
    ├── irrigation_mask_04.csv              # 0.4° resolution
    └── irrigation_mask_1.csv              # 1° resolution
```

## Key analytical concepts

### Detection threshold τ
Rather than using a binary presence/absence definition, the analysis sweeps over a continuous minimum-area threshold τ (in Mha per cell). A cell is considered *irrigated* by a given dataset only if its reported area exceeds τ. This allows quantifying how disagreement depends on the minimum irrigable area that is considered meaningful. Cell-area nominal values used to convert percentage thresholds to Mha are: 0.05 Mha at 0.2°, 0.10 Mha at 0.4°, and 0.30 Mha at 1°.

### Agreement level *k*
The variable *k* counts how many of the ten datasets simultaneously classify a cell as irrigated. The masks in the `irrigation_masks/` folder encode *k* for every grid cell at each resolution. Enforcing *k* ≥ 2 (or higher) progressively restricts the area considered "confirmed" irrigated land.

### Instability and flip rate
The `non_identifiability_tests.R` functions measure how often a cell's irrigation classification flips when: (1) the dataset ensemble is perturbed by leaving one dataset out, or (2) the threshold τ is varied. The *flip rate* is defined as `1 − max(share_present, share_absent)` across all perturbation runs, so a flip rate of 0 means perfectly stable and 0.5 means maximally unstable.

## Irrigated area masks

We include the irrigation masks resulting from our study. The masks indicate how
many datasets agree that there is irrigation in a given grid cell (detection
limit set at τ > 1% of the grid cell), from (*k* = 0; 0 datasets) to (*k* = 10; 10 datasets).
They can be found in the `irrigation_masks/` folder, with one file per resolution:

* `irrigation_mask_02.csv` (0.2°).
* `irrigation_mask_04.csv` (0.4°).
* `irrigation_mask_1.csv` (1°).

Columns:

| Column | Description |
|---|---|
| `lon` | Longitude of the grid-cell centre (decimal degrees, WGS84) |
| `lat` | Latitude of the grid-cell centre (decimal degrees, WGS84) |
| `country` | Country name (Natural Earth, medium scale) |
| `code` | ISO numeric country code |
| `continent` | Continent name |
| `k` | Number of datasets (0–10) that classify the cell as irrigated |

Only cells where at least one dataset reports non-zero irrigation are included (all-zero cells are dropped). To obtain a consensus irrigated-area map for a chosen agreement level, filter by `k >= threshold`.

**Example in R:**

```r
library(data.table)
mask <- fread("irrigation_masks/irrigation_mask_02.csv")

# Cells where at least 5 of 10 datasets agree
consensus_5 <- mask[k >= 5]

# Country totals at majority agreement
consensus_5[, .N, country][order(-N)]
```

## Citation

If you use this workflow or the irrigated area masks, please cite:

A. Puy, O. Richards, S. N. Linga, S. Flinders, C. Aguiló-Rivera, J. Larsen. (2026).
Code and datasets of Where irrigation exists is globally contested. Zenodo. doi: 10.5281/zenodo.19001232.

## License

MIT License

Copyright (c) 2026 Arnald Puy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
