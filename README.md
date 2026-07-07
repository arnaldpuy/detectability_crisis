[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19001232.svg)](https://doi.org/10.5281/zenodo.19001232)

# Where irrigation exists is globally contested

[Arnald Puy](https://www.arnaldpuy.com/), Olivia Richards, Seth N. Linga, Samuel
Flinders, Carmen Aguiló-Rivera, Josh Larsen

This study analyzes the disagreement across datasets in the identification of
global irrigated areas.

## Abstract

*Irrigation is one of the most extensive human modifications of the land surface and a key driver of food production, water use and regional climate. However, the location of irrigated land has never been systematically verified at global scale. Here we show that where irrigation exists is globally contested. Across ten datasets, 60–90\% of cropland grid cells disagree on whether irrigation exists at all, remaining at 40–60\% when irrigation is required to exceed 1\% of the cell to confirm presence. Enforcing agreement across datasets reduces the irrigated area of irrigation hotspots such as Vietnam, China or Thailand by 30-80\%, the share of global crop production associated with irrigation by 60-90\% and the global irrigation-induced evapotranspiration signal by $\sim67$\%. The consequences extend to water governance, where discrepancies across datasets reclassify the Sustainable Development Goal 6.4.2 water-stress status of 42\% of countries. The disagreement persists regardless of which datasets, methods or resolutions are considered and is consistent with irrigated areas being a non-identifiable variable. The way forward requires acknowledging that whether an area is irrigated is a judgement call that data alone cannot resolve.*

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

## 1. System requirements

### Operating systems
The workflow is implemented entirely in **R** and uses no operating-system-specific
calls. It was developed and tested on:

* **macOS 26.5.1** (Apple Silicon, `arm64`; build 25F80).

It is expected to run on any recent Linux or Windows system with a working R
installation, but it has only been tested on the macOS configuration above.

### Software dependencies
* **R ≥ 4.5** (tested on **R 4.5.2**, `aarch64-apple-darwin20`).
* The R packages below. `parallel` ships with base R; all others are on CRAN. The
  right-hand column lists the exact versions the workflow was tested with.

| Package | Purpose | Tested version |
|---|---|---|
| `data.table` | Fast tabular data manipulation | 1.18.4 |
| `terra` | Raster reading, reprojection and aggregation | 1.8.93 |
| `sf` | Vector spatial data (country polygons) | 1.1.0 |
| `tidyverse` | Data wrangling and plotting utilities | 2.0.0 |
| `ggplot2` | Visualisation (loaded via `tidyverse`) | 4.0.3 |
| `cowplot` | Multi-panel figure composition | 1.2.0 |
| `scales` | Axis and colour scale helpers | 1.4.0 |
| `here` | Portable file paths | 1.0.2 |
| `rnaturalearth` | Country boundary data | 1.1.0 |
| `countrycode` | Country name and code standardisation | 1.6.1 |
| `sp` | Legacy spatial utilities (original dataset processing) | 2.2.1 |
| `rworldmap` | Legacy spatial utilities (original dataset processing) | 1.3.8 |
| `readxl` | Reading the dataset-survey Excel file | 1.4.5 |
| `ncdf4` | Reading NetCDF climate model output | 1.24 |
| `wesanderson` | Colour palettes | 0.3.7 |
| `benchmarkme` | System information for reproducibility reporting | 1.0.8 |
| `sensobol` | Uncertainty and sensitivity analyses | 1.2.0 |
| `jsonlite` | Parsing the UNSD / World Bank API responses | 2.0.0 |
| `magrittr` | Pipe operator | 2.0.5 |
| `parallel` | Parallel computation | base R |

(`ggplot2` and `cowplot` were tested with current development builds of 4.0.3 and
1.2.0; the released versions behave identically for this workflow. A reproducibility
seed of `123` is set at the start of each script.)

### Hardware
No non-standard hardware is required; the workflow runs on a standard desktop or
laptop. Because the full pipeline ingests large global rasters and a harmonised
ensemble of ~130,000 grid cells, **≥16 GB RAM is recommended** for the
dataset-ingestion and harmonisation steps (the test machine had 64 GB). The demo
and the consensus masks run comfortably in under 2 GB.

## 2. Installation guide

1. Install **R ≥ 4.5** from [CRAN](https://cran.r-project.org/) (optionally with
   [RStudio](https://posit.co/download/rstudio-desktop/)).
2. Download or clone this repository:
   ```bash
   git clone https://github.com/arnaldpuy/detectability_crisis.git
   ```
   or download the archive from [Zenodo](https://doi.org/10.5281/zenodo.19001232).
3. Install the R dependencies:
   ```r
   install.packages(c(
     "data.table", "terra", "tidyverse", "cowplot", "scales", "here", "sf",
     "rnaturalearth", "countrycode", "sp", "rworldmap", "readxl", "ncdf4",
     "wesanderson", "benchmarkme", "sensobol", "magrittr", "jsonlite"
   ))
   ```

**Typical install time.** On a normal desktop with a broadband connection,
installing the full package set takes roughly **10–30 minutes**. Most packages
install in seconds from pre-compiled binaries; the time is dominated by `terra`,
`sf` and `sensobol`, which may compile from source and depend on the system
libraries GDAL, GEOS and PROJ (normally bundled with the binary builds).

## 3. Demo

The consensus irrigation masks shipped in `irrigation_masks/` provide a
self-contained demo that needs no external data. The snippet below counts the grid
cells where at least five of the ten datasets agree on irrigation presence and
lists the leading countries:

```r
library(data.table)
mask <- fread("irrigation_masks/irrigation_mask_02.csv")   # 0.2 deg, 129,628 cells

consensus_5 <- mask[k >= 5]              # cells where >= 5 of 10 datasets agree
nrow(consensus_5)                        # -> 46920

consensus_5[, .N, country][order(-N)][1:6]
```

**Expected output:**

```
> nrow(consensus_5)
[1] 46920

         country     N
1:         China  9390
2:         India  5595
3: United States  5269
4:          Iran  1766
5:        Brazil  1675
6:        Mexico  1276
```

**Expected run time:** under 1 second on a normal desktop.

For a larger, figure-producing demo, knit `code_reply_to_reviewers.Rmd`: it
reproduces every quantitative result in our reply to the reviewers from the
harmonised 0.2° ensemble (run time ~1–2 minutes; requires the `datasets/` folder
from the Zenodo archive).

## 4. Instructions for use

### Running the masks on your own data

The consensus masks in `irrigation_masks/` are the main reusable product of the
study. To use them in your own analysis, read the mask at the desired resolution
and keep the cells that meet your agreement criterion *k*. For example,
`mask[k >= 5]` retains only cells where at least half of the ten datasets agree
that irrigation is present; join on `lon`/`lat` (grid-cell centres, WGS84) to mask
your own gridded data. See [Irrigated area masks](#irrigated-area-masks) below for
the full column definitions and a worked example.

### Reproduction instructions

The full study is reproduced by running the scripts in the order below. Each is
provided as `.R`, `.Rmd` and a rendered `.pdf`. The raw rasters and the
intermediate `datasets/` folder are archived on
[Zenodo](https://doi.org/10.5281/zenodo.19001232); the consensus masks ship with
this repository. All functions in `functions/` are sourced automatically, so they
do not need to be loaded by hand.

| Order | Script | Approximate run time |
|---|---|---|
| 1 | `code_original_datasets` | hours (ingests and converts all raw rasters) |
| 2 | `code_harmonization` | ~1–2 hours |
| 3 | `code_main_analysis` | several hours (full τ sweep, perturbation and ET analyses) |
| 4 | `code_sdg_analysis` | ~1–2 minutes |
| 5 | `code_reply_to_reviewers` | ~1–2 minutes |

Run times are indicative for the test machine (Apple Silicon, 64 GB RAM). Steps 1–3
are the heavy raster-processing stages; steps 4–5 run from the cached harmonised
ensemble and are fast. The two descriptive blocks below detail the custom functions
and each script.

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

#### 4. `code_sdg_analysis`
Propagates the detectability crisis into the two FAO water-governance indicators **SDG 6.4.1** (change in water-use efficiency) and **SDG 6.4.2** (level of water stress). Rather than reimplementing AQUASTAT or GlobWat, it takes FAO's *published* indicator values and **swaps only the irrigation-derived input**: the GMIA-derived agricultural-withdrawal term is replaced by each of the ten maps (and by the *k*-of-10 consensus masks), holding all other inputs fixed and using the approximately linear scaling of agricultural water withdrawal with irrigated area documented by [Puy et al](https://www.nature.com/articles/s41467-021-24508-8). Official inputs are retrieved from the UNSD SDG API and the World Bank; the five-class water-stress thresholds follow the SDG 6.4.2 metadata (FAO / UN Statistics Division). The script reports how many countries change their official water-stress band across the ensemble, contrasts this with the (band-less) efficiency indicator, and produces the SDG propagation figures.

#### 5. `code_reply_to_reviewers`
Reproduces every quantitative figure cited in our reply to the reviewers, organised by reviewer. Working from the harmonized 0.2° ensemble restricted to cropland cells, it reports: the detectability threshold *τ* as a fraction of the grid-cell area at each resolution (Reviewer 1); the temporal-heterogeneity tests showing that maps representing the same nominal year disagree on the existence of irrigation as much as maps spanning 1999–2015, and that the most recent maps record no irrigation where older maps detect it (Reviewer 1); the existential disagreement among products that map the same variable — including that the five census-based products are the single most-agreeing five-map subset, that remote-sensing and census products are not a subset of one another, and which variable is extracted for each product (Reviewer 2); and the existential and extreme disagreement by country, showing it does not track irrigation practice (Reviewer 3).

The notebook closes with three additional analyses supporting the non-identifiability claim (first round of revisions): an underdetermination count, a latent-class test of the "identifiable truth plus noise" hypothesis and a comparison of the best US regional products (LANID, MIrAD-US, AIM-HPA). The rationale and interpretation are given in the manuscript and its Supplementary Materials. The regional products aggregated to the 0.2° grid ship in `us_regional_products/us_regional_products_02.csv`, so the US analysis reproduces without downloading the raw rasters.

The CSV is a derivative (aggregation to 0.2°) of third-party data redistributed under their own terms: LANID ([Xie & Lark 2021](https://doi.org/10.5281/zenodo.5548555), CC-BY 4.0) and AIM-HPA ([Deines et al. 2019](https://doi.org/10.4211/hs.a371fd69d41b4232806d81e17fe4efcb), CC-BY 4.0), both modified by aggregation; and MIrAD-US ([USGS](https://doi.org/10.5066/P9NA3EO8), US public domain).

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
├── code_sdg_analysis.{Rmd,pdf}            # Propagation into SDG 6.4.1 / 6.4.2 indicators
├── code_reply_to_reviewers.{Rmd,pdf}      # Analyses reproducing the reply to the reviewers
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
├── irrigation_masks/                       # Consensus irrigated-area masks
│   ├── irrigation_mask_02.csv              # 0.2° resolution
│   ├── irrigation_mask_04.csv              # 0.4° resolution
│   └── irrigation_mask_1.csv               # 1° resolution
└── us_regional_products/                   # US regional products at 0.2° (reply to reviewers)
    └── us_regional_products_02.csv         # LANID, MIrAD-US, AIM-HPA (2012, 2017) at US cropland cells
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
