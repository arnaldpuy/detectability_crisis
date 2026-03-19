[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19001232.svg)](https://doi.org/10.5281/zenodo.19001232)
# Where irrigation exists is globally contested

[Arnald Puy](https://www.arnaldpuy.com/), Olivia Richards, Seth N. Linga, Samuel
Flinders, Carmen Aguiló-Rivera

This study analyzes the disagreement across datasets in the identification of 
global irrigated areas.

## Abstract

*Irrigation is one of the most extensive human modifications of the land surface and a key driver of food production, water use and regional climate. However, the global location of irrigated land has never been empirically verified. Here we show that irrigation presence is not observationally identifiable at planetary scale. Across ten global datasets, 60-90% of cropland grid cells disagree on whether irrigation exists at all, and these contradictions persist across resolutions and detection thresholds. Restricting analysis to areas with confirmed irrigation by all datasets makes irrigation hotspots such as Vietnam, Bangladesh or Thailand drop by more than 20 positions in global rankings and reduces global crop production and irrigation-induced evapotranspiration estimates by up to 65%. Global policies and food, water and climate estimates rely on spatial representations of irrigated areas whose very existence is unresolved.*

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

## Replication

We provide all the functions needed to replicate our workflow in the "functions" folder. 

#### Generated data

The `datasets/irrigated_areas_regridded` folder contains the data produced and analysed in this study. 

* `irrigated_areas_regridded` folder: It contains the ten maps mentioned above
harmonized to the following resolutions: 

  - `irrigated_areas_regridded_02.csv`: at a 0.2º resolution.
  - `irrigated_areas_regridded_04.csv`: at a 0.4º resolution.
  - `irrigated_areas_regridded_10.csv`: at a 1º resolution.
  
The datasets are in long format. To turn them into wide format, the user can
run the following code snippet in `R` (example with the dataset at 0.2º):
  
```r
library(data.table)

dt <- fread("irrigated_areas_regridded_02.csv")

dt_wide <- dcast(dt, lon + lat + country + code + continent + resolution ~ dataset, value.var = "mha")

```

### Functions

The `functions` folder contains all the custom functions coded for the analysis.
They are all sourced from the `.R`, `.pdf` and `.Rmd` files and therefore the 
user of the code does not need to source them separately.

### Code

We offer the code in `.R`, `.pdf` and `.Rmd`. There are four main analyses:

* 1. `code_original_datasets`: workflow to transform all rasters from `tif` to tabular form.
* 2. `code_harmonization`: harmonization of all irrigated area datasets at 0.2º, 0.4º and 1º.
* 3. `code_crop_maps`: extraction of global irrigated crop production from [GAEZ+2015](https://www.nature.com/articles/s41597-021-01115-2).
* 4. `code_main_analysis`: analysis of the detectability of irrigated areas at the grid cell level.

