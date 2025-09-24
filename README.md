# dsc: Dynamic Synthetic Control for Time Series with Heterogeneous Adjustment Speeds

[![GitHub License](https://img.shields.io/github/license/conflictlab/dsc)](https://github.com/conflictlab/dsc/blob/main/LICENSE) 
![GitHub Repo stars](https://img.shields.io/github/stars/conflictlab/dsc?style=social) 

Implements the Dynamic Synthetic Control (DSC) method as described in Cao and Chadefaux (2024). DSC extends the synthetic control method (Abadie and Gardeazabal, 2003; Abadie, Diamond, and Hainmueller, 2010, 2011, 2014) by accounting for varying response speeds across units and time. Using dynamic time warping (DTW), DSC aligns donor units with the treated unit in the pre-treatment period, mitigating bias caused by differential adaptation rates. It provides a data-driven procedure to construct synthetic control units based on a weighted combination of comparison units, improving counterfactual estimation in comparative case studies. The package includes functions for DSC implementation, placebo tests, and visualization tools for treatment effect inference.

Authors: Jian Cao and Thomas Chadefaux

Maintainer: Jian Cao. caoj@tcd.ie

License: GPL-3

Depends: R (>= 3.5.0)

Imports: 

URL: https://github.com/conflictlab/dsc


## Installation

You can install `dsc` directly from GitHub with:

```R
devtools::install_github("conflictlab/dsc")
```

## Example Usage

Here's a basic example to get you started:

```R
# Load the dsc package
library(dsc)

# Load the Basque dataset from the Synth package
data(basque, package = "Synth")
data <- basque

# Rename relevant columns for clarity
colnames(data)[1:4] <- c("id", "unit", "time", "value")

# Compute additional variables
data$invest_ratio <- data$invest / data$value
data$value_raw <- data$value

# Define special predictors for the model
special_preds <- expression(list(
  list(dep.var, 1960:1969, c("mean")),
  list("invest_ratio", 1964:1969, c("mean")),
  list("popdens", 1969, c("mean")),
  list("sec.agriculture", 1961:1969, c("mean")),
  list("sec.energy", 1961:1969, c("mean")),
  list("sec.industry", 1961:1969, c("mean")),
  list("sec.construction", 1961:1969, c("mean")),
  list("sec.services.venta", 1961:1969, c("mean")),
  list("sec.services.nonventa", 1961:1969, c("mean")),
  list("school.illit", 1964:1969, c("mean")),
  list("school.prim", 1964:1969, c("mean")),
  list("school.med", 1964:1969, c("mean")),
  list("school.high", 1964:1969, c("mean")),
  list("school.post.high", 1964:1969, c("mean"))
))

# Execute the DSC analysis
result <- dsc(
  data = data,
  start.time = 1955,
  end.time = 1997,
  treat.time = 1970,
  dependent = "Basque Country (Pais Vasco)",
  predictors = NULL,
  parallel = TRUE,
  special.predictors = special_preds,
  time.predictors.prior = 1955:1969,
  time.optimize.ssr = 1955:1969
)
```

## Citation
To acknowledge the use of this package in your work, please use the citation format provided below:

Jian Cao (2023). dsc: Dynamic Synthetic Control Method for Comparative Case Studies. Version 0.1.0. URL: github.com/conflictlab/dsc




## Feedback and Contributions

We welcome your feedback and contributions! Please open an issue on [GitHub](https://github.com/conflictlab/dsc/issues) or submit a pull request if you want to propose changes.
