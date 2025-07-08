# Anaemia Malaria Fitting

This repository contains R scripts and vignettes for modeling the impact of malaria on haemoglobin levels in pregnancy. The vignettes demonstrate how to generate synthetic data and how to run the inferential framework on both censored and uncensored datasets.

## Requirements

Make sure you have R (version 4.0 or later) with the packages listed in `R/load_libraries.R` installed. The main packages used are `drjacoby`, `ggplot2`, `dplyr`, `patchwork` and their dependencies. You can install them from CRAN:

```r
install.packages(c("drjacoby", "ggplot2", "dplyr", "patchwork", "easypackages", "tidyr", "gridExtra", "ggnewscale", "ggpubr", "viridis", "matrixStats", "readxl", "grid", "DescTools", "ggtext"))
```

## Running the vignettes

1. Clone this repository and open the project in RStudio or your preferred R environment.
2. Source the helper scripts so that functions are available:

```r
source("R/load_libraries.R")
source("R/cubic_spline.R")
source("R/immune_functions.R")
source("R/prob_censor.R")
source("R/data_generation_functions.R")
source("R/define_params_functions.R")
source("R/log_likeli.R")
source("R/log_prior.R")
source("R/plot_functions.R")
```

3. Knit the vignettes using `rmarkdown::render` or by clicking *Knit* in RStudio.

```r
rmarkdown::render("vignette/data_generation.Rmd")
rmarkdown::render("vignette/malaria_hb_vignette.Rmd")
```

Rendered HTML versions (`*.html`) are provided in the `vignette/` directory for quick browsing.

## Example data

Pre-computed example outputs and infection history data used by the vignettes are stored in `inst/extdata/` and loaded via `readRDS` within the R Markdown files.

## License

This project is released under the MIT License. See `LICENSE` for details.
