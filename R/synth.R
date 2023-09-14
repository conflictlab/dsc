#' Apply the Synth Method to Time-Series Data (Internal)
#'
#' This function performs the synthetic control method on provided data. It prepares the data
#' and fits a synthetic control model using the `synth` package. This is an internal helper
#' function and is not intended to be directly used by package users.
#'
#' @param df A data frame containing the data for analysis.
#' @param dep.var A string indicating the dependent variable column name.
#' @param dependent.id An identifier for the treated unit or case.
#' @param predictors A list or vector of predictor variable names.
#' @param special.predictors A list or vector of special predictor names.
#' @param time.predictors.prior Time periods for which to use the predictors.
#' @param time.optimize.ssr Time period to optimize the SSR.
#'
#' @return A list containing the actual values (`value`), average values (`average`),
#'         and the synthetic control predictions (`synthetic`).
#'
#' @keywords internal
do.synth = function(df, dep.var, dependent.id, predictors,
                    special.predictors, time.predictors.prior,
                    time.optimize.ssr){
  # find v
  dataprep.out <-
    Synth::dataprep(
      foo = df,
      predictors    = eval(predictors),
      dependent     = dep.var,
      unit.variable = 1,
      time.variable = 3,
      special.predictors = eval(special.predictors),
      treatment.identifier = dependent.id,
      controls.identifier = setdiff(unique(df$id), dependent.id),
      time.predictors.prior = time.predictors.prior,
      time.optimize.ssr = time.optimize.ssr,
      unit.names.variable = 2,
      time.plot = sort(unique(df$time))
    )

  # fit training model
  synth.out <-
    Synth::synth(
      data.prep.obj=dataprep.out,
      Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
    )

  value = df %>% filter(id == dependent.id) %>% `$`(value_raw)
  average = df %>% filter(id != dependent.id) %>% group_by(time) %>%
    summarise(average = mean(!!sym(dep.var), na.rm = TRUE)) %>% `$`(average)
  synthetic = dataprep.out$Y0%*%synth.out$solution.w %>% as.numeric

  return(list(value = value,
              average = average,
              synthetic = synthetic))
}
