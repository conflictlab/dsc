---

title: "dsc: Dynamic Synthetic Control for Time Series with Heterogeneous Adjustment Speeds"
authors:
  - given-names: Jian
    surname: Cao
    affiliation: 1
  - given-names: Thomas
    surname: Chadefaux
    affiliation: 1


---

# Summary

The `dsc` package introduces Dynamic Synthetic Control, a new approach for comparative case studies in time series settings. Synthetic control methods are widely used to estimate causal effects, but they often fail when treated and donor units react to shocks at different speeds. The `dsc` package addresses this by incorporating Dynamic Time Warping (DTW) to account for heterogeneous adjustment speeds across units, improving counterfactual estimation.

Implemented in R, `dsc` aligns donor units to the treated unit using speed-adjusted time series before estimating synthetic weights. This innovation helps avoid biases arising from asynchronous reactions and supports more accurate estimation of treatment effects. The package is useful for applications in economics, public policy, and political science.

# Statement of Need

Standard synthetic control techniques assume that all units react to shocks or policies at the same speed. This can result in poor fits and misleading conclusions when donor units respond more slowly or more quickly than treated ones. For example, institutional inertia might delay reactions in one country relative to another, even if underlying economic mechanisms are similar.

The `dsc` package provides a principled solution by using DTW to synchronize pre-treatment time series between treated and donor units. This synchronization reduces mean squared error in treatment effect estimation by up to 70% in simulations and improves placebo test performance in real-world datasets. It fills a gap in the causal inference toolkit by allowing for varying speeds of adjustment, a common real-world phenomenon that existing packages ignore.

# Model Overview

Synthetic control methods construct counterfactuals for treated units using weighted combinations of untreated donor units. Let \$y\_{1t}\$ denote the treated unit, and \$y\_{jt}\$ denote donors \$j = 2, ..., J+1\$. The goal is to find weights \$w\_j\$ such that:

$$
y_{1t} \approx \sum_{j=2}^{J+1} w_j y_{jt}
$$

for the pre-treatment period \$t < T\$.

However, if donor units respond to latent shocks \$z\_t\$ with lags, then the pre-treatment series are not aligned in time. The `dsc` package addresses this by warping \$y\_{jt}\$ to align with \$y\_{1t}\$ using DTW. The warped donor series \$y^w\_{jt}\$ is then used in synthetic control estimation.

# Implementation

The core of the `dsc` method is a three-step process:

1. **Warping Pre-Treatment Series**: Use DTW to align each donor series \$y\_j\$ to the treated unit \$y\_1\$ during the pre-treatment period.
2. **Propagate Speed Alignment**: Apply the inferred warping path to the post-treatment period of each donor unit.
3. **Construct Synthetic Control**: Estimate weights \$w\_j\$ to best fit the warped donor series \$y^w\_j\$ to \$y\_1\$ before treatment.

This preserves any speed differences introduced by the treatment itself, while eliminating those inherited from structural or institutional differences.

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
  time.optimize.ssr = 1955:1969
)

# Visualize results
plot(result)
```

# Empirical Applications

## Terrorism and GDP in the Basque Country

We replicate Abadie and Gardeazabal (2003), estimating the effect of terrorism on GDP using DSC. Compared to traditional synthetic control, DSC shows a closer match for placebo units and reduced mean squared error.

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

The average \$r\$ is negative across all scenarios, indicating that DSC yields lower mean squared errors.

# Discussion and Limitations

The `dsc` method improves synthetic control estimation by accounting for reaction speed heterogeneity. However, it assumes that speed differences are stable in the pre-treatment period and that no spillover effects contaminate donor units. Extensions to multi-treatment cases or endogenously determined timing remain areas for future work.

# Acknowledgements

This research was supported by the European Research Council (ERC) under the European Union's Horizon 2020 programme (grant agreement No. 101002240).

# References

```bibtex
@article{abadie2003economic,
  title={The economic costs of conflict: A case study of the Basque Country},
  author={Abadie, Alberto and Gardeazabal, Javier},
  journal={American Economic Review},
  volume={93},
  number={1},
  pages={113--132},
  year={2003}
}

@article{abadie2010synthetic,
  title={Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program},
  author={Abadie, Alberto and Diamond, Alexis and Hainmueller, Jens},
  journal={Journal of the American Statistical Association},
  volume={105},
  number={490},
  pages={493--505},
  year={2010}
}

@article{abadie2015comparative,
  title={Comparative politics and the synthetic control method},
  author={Abadie, Alberto and Diamond, Alexis and Hainmueller, Jens},
  journal={American Journal of Political Science},
  volume={59},
  number={2},
  pages={495--510},
  year={2015}
}

@article{vintsyuk1968dtw,
  title={Speech discrimination by dynamic programming},
  author={Vintsyuk, T. K.},
  journal={Cybernetics},
  volume={4},
  number={1},
  pages={52--57},
  year={1968}
}
```
