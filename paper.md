---
title: 'dsc: Dynamic Synthetic Control for Time Series with Heterogeneous Adjustment
  Speeds'
authors:
- given-names: Jian
  surname: Cao
  affiliation: 1
- given-names: Thomas
  surname: Chadefaux
  affiliation: 1
date: "14 May 2025"
output: pdf_document
affiliations:
- index: 1
  name: Department of Political Science, Trinity College Dublin, Ireland
---


# Summary

The dsc package (Cao and Chadefaux, 2025) implements the Dynamic Synthetic Control (DSC) method, an extension of the synthetic control framework designed for comparative case studies involving time series data. Traditional synthetic control assumes that treated and donor units respond to shocks at the same rate, which can lead to poor pre-treatment fit and biased treatment effect estimates when this assumption is violated.

The DSC method addresses this issue by aligning time series data from donor units with the treated unit using Dynamic Time Warping (DTW), a technique that accounts for differences in the speed of adjustment across units. The adjustment is only made to the misalignment that originates from the pre-treatment period, preserving the integrity of treatment effects post-intervention. This correction improves the counterfactual estimation, as shown in both empirical examples and simulations.

The method is described in detail in Cao and Chadefaux (2024), where it is applied to comparative political economy questions in a peer-reviewed context. The dsc package provides a clean and replicable implementation of the methodology introduced in that work, extending its accessibility to empirical researchers.

The dsc package is implemented in R and integrates easily into existing workflows. It supports parallel computation, placebo tests, and multiple diagnostic plots. Its interface is user-friendly and familiar to users of the Synth package. DSC is particularly useful in political science, economics, and public policy, where variation in adjustment dynamics is common.

# Statement of Need

Synthetic control methods are now standard tools in causal inference, especially for policy evaluation. However, the standard approach assumes all units respond to shocks or treatments at a uniform rate. This is rarely true in practice: countries, regions, or institutions often differ in how quickly they adjust to change. If this heterogeneity is not accounted for, it may lead to incorrect counterfactuals and misleading conclusions.

The dsc package fills this gap by incorporating time-series alignment using DTW before the construction of the synthetic control. This ensures that donor units are temporally matched with the treated unit, even if their dynamics differ. Unlike some other synthetic control extensions that address poor pre-treatment fit by removing problematic donors or augmenting the model, DSC fixes the problem at the alignment stage without discarding information.

No other R package currently offers DTW-based alignment within the synthetic control framework in a fully automated and documented way. Thus, dsc represents a substantial contribution to applied research workflows.

# Model Overview

The synthetic control estimator builds counterfactuals for treated units using weighted combinations of untreated donor units. Let $y_{1t}$ denote the treated unit, and $y_{jt}$ denote donors $j = 2, ..., J+1$. The goal is to find weights $w_j$ such that:

$$
y_{1t} \approx \sum_{j=2}^{J+1} w_j y_{jt}
$$

for the pre-treatment period $t < T$.

However, if donor series differ in response speed, direct comparison misaligns their dynamics. The DSC method introduces DTW-based warping functions $g_j(t)$ such that the warped donor series $y_{jt}^{(w)} = y_{j, g_j(t)}$ is aligned with the treated unit's trajectory in the pre-treatment window. Synthetic weights are then computed using these aligned donor series.

The warping is applied only to misalignments that originate from the pre-treatment period. In the post-treatment window, the same transformation is propagated, which ensures comparability and interpretability.

# Implementation

The core of the `dsc` method is a three-step process:

1. **Warping Pre-Treatment Series**: Use DTW to align each donor series $y_j$ to the treated unit $y_1$ during the pre-treatment period.
2. **Propagate Speed Alignment**: Apply the inferred warping path to the post-treatment period of each donor unit.
3. **Construct Synthetic Control**: Estimate weights $w_j$ to best fit the warped donor series $y^w_j$ to $y_1$ before treatment.

This preserves any speed differences introduced by the treatment itself, while eliminating those inherited from structural or institutional differences.

The package also includes:

Diagnostic plots (e.g., observed vs synthetic, treatment gap)

Placebo tests

Parallel processing for faster estimation

# Code Example

The `dsc` package (Cao and Chadefaux, 2024) can be installed from GitHub using:

```r
devtools::install_github("conflictlab/dsc")
```

Here’s a representative use case based on the Basque Country dataset:

```r
library(dsc)
library(Synth)

# Load dataset
data(basque, package = "Synth")
data <- basque

# Prepare data
colnames(data)[1:4] <- c("id", "unit", "time", "value")
data$invest_ratio <- data$invest / data$value

# Specify special predictors
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

# Run DSC
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
  time.optimize.ssr = 1955:1969,
  plot.figures=TRUE
)

```

![Observed GDP of the Basque Country compared with its DSC-based synthetic control.](example.pdf){#fig:basque-dsc width=100%}


# Empirical Applications

## Terrorism and GDP in the Basque Country

We replicate Abadie and Gardeazabal (2003), estimating the effect of terrorism on GDP using DSC. Compared to traditional synthetic control, DSC shows a closer match for placebo units and reduced mean squared error (Figure @fig:basque-dsc).


## Proposition 99: Tobacco Control in California

We revisit the effect of California’s anti-smoking policy, Proposition 99. DSC estimates a larger reduction in cigarette consumption and outperforms the original model on placebo test sharpness.

## German Reunification

We assess the impact of reunification on West Germany’s GDP. Again, DSC improves the counterfactual fit for placebo countries and yields more precise treatment estimates.

# Monte Carlo Evaluation

To validate performance, we simulate data where units react to shocks at varying speeds. Across 100 replications, DSC consistently produces treatment effect estimates with lower variance and bias compared to standard synthetic control.

We define the relative improvement as:

$$
r = \log \left( \frac{\text{MSE}_{DSC}}{\text{MSE}_{SC}} \right)
$$

The average $r$ is negative across all scenarios, indicating that DSC yields lower mean squared errors.

# Discussion and Limitations

The `dsc` method improves synthetic control estimation by accounting for reaction speed heterogeneity. Extensions to multi-treatment cases or endogenously determined timing remain areas for future work.

# Acknowledgements

This research was supported by the European Research Council (ERC) under the European Union's Horizon 2020 programme (grant agreement No. 101002240).

# References

Abadie, A., & Gardeazabal, J. (2003). *The economic costs of conflict: A case study of the Basque Country*. American Economic Review, 93(1), 113–132.

Abadie, A., Diamond, A., & Hainmueller, J. (2010). *Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program*. Journal of the American Statistical Association, 105(490), 493–505.

Abadie, A., Diamond, A., & Hainmueller, J. (2015). *Comparative politics and the synthetic control method*. American Journal of Political Science, 59(2), 495–510.

Cao, J., & Chadefaux, T. (2024). Dynamic Synthetic Control. Political Analysis, 32(1), 1–18. https://doi.org/10.1017/pan.2023.50

Vintsyuk, T. K. (1968). *Speech discrimination by dynamic programming*. Cybernetics, 4(1), 52–57.

Cao, J., & Chadefaux, T. (2025). dsc: Dynamic Synthetic Control for Time Series with Heterogeneous Adjustment Speeds [Computer software]. GitHub. https://github.com/conflictlab/dsc
