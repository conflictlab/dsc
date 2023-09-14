#' Dynamic Synthetic Control (DSC) method
#'
#' The DSC method combines ideas from the synthetic control method and dynamic time warping.
#' This function processes the given data, applies time-series preprocessing, and then
#' computes the synthetic control using dynamic time warping.
#'
#' @param data A data frame containing the observational data.
#' @param start.time Starting time for the analysis.
#' @param end.time Ending time for the analysis.
#' @param treat.time Treatment time.
#' @param dependent The dependent variable name in the dataset.
#' @param k Integer, number of control units used for dynamic time warping. Default is 4.
#' @param filter.width Integer, width of the filter. Default is 5.
#' @param buffer Integer, buffer for time series alignment. Default is 0.
#' @param norm.method Method for normalization. Default is "t".
#' @param match.method Method for matching. Default is "fixed".
#' @param step.pattern1 Step pattern for the DTW. Default is dtw::symmetricP1.
#' @param step.pattern2 Alternative step pattern for DTW. Default is dtw::asymmetricP2.
#' @param plot.figures Logical, if TRUE plots will be generated. Default is FALSE.
#' @param n.burn Integer, number of initial time periods to disregard. Default is 3.
#' @param ma Integer, moving average length. Default is 3.
#' @param ma.na Method to handle missing values in moving average. Default is "original".
#' @param dist.quant Numeric, quantile for distance measure. Default is 1.
#' @param n.IQR Numeric, factor for IQR in outlier detection. Default is 3.
#' @param window.type Type of window for DTW. Default is "none".
#' @param default.margin Default margin size. Default is 3.
#' @param n.q Integer, number of synthetic controls to use. Default is 1.
#' @param n.r Integer, number of predictors to use. Default is 1.
#' @param parallel Logical, if TRUE parallel processing will be enabled. Default is TRUE.
#' @param rescale Logical, if TRUE data will be rescaled. Default is TRUE.
#' @param dependent.id Numeric, ID of the dependent unit.
#' @param predictors List, names of predictor variables.
#' @param special.predictors List, names of special predictor variables.
#' @param time.predictors.prior List, names of time predictor variables for the prior period.
#' @param time.optimize.ssr List, names of time predictor variables for the SSR optimization.
#'
#' @return A list containing results of the synthetic control analysis.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(...)
#' result <- dsc(data, ...)
#' }
#'
#' @export
dsc <- function(data, start.time, end.time, treat.time,
                dependent,
                k = 4, filter.width = 5,
                buffer = 0, norm.method = "t",
                match.method = "fixed",
                step.pattern1 = dtw::symmetricP1,
                step.pattern2 = dtw::asymmetricP2,
                plot.figures = FALSE, n.burn = 3,
                ma = 3, ma.na = "original",
                dist.quant = 1, n.IQR = 3,
                window.type = "none",
                default.margin = 3,
                n.q = 1, n.r = 1,
                parallel = TRUE, rescale = TRUE,
                dependent.id, predictors,
                special.predictors, time.predictors.prior,
                time.optimize.ssr){

  # Preprocess data
  if (rescale) {
    rescale_data <- data %>%
      filter(time <= treat.time) %>%
      group_by(unit) %>%
      summarize(
        value.min = min(value),
        value.max = max(value)) %>%
      ungroup()

    mean.diff = mean(rescale_data$value.max - rescale_data$value.min)
    rescale_data = rescale_data %>%
      mutate(
        multiplier = mean.diff/(value.max - value.min)
      )

    data <- left_join(data, rescale_data, by = "unit") %>%
      mutate(
        value.bak = value_raw,
        value_raw = (value_raw - value.min) * multiplier,
        value = value_raw
      )
  }

  # Filter data
  if (!is.null(filter.width)) {
    data = preprocessing(data, filter.width)
  }

  # Prepare data for TFDTW
  t.treat <- (treat.time - start.time) + 1

  y_values <- data %>%
    filter(unit == dependent & time <= end.time) %>%
    select(value_raw, value)

  y.raw <- y_values[["value_raw"]]
  y.processed <- y_values[["value"]]

  x.list <- data %>%
    filter(unit != dependent & time <= end.time) %>%
    select(unit, time, value, value_raw) %>%
    group_by(unit) %>%
    group_split(.keep = TRUE)

  # TFDTW
  fun_map <- ifelse(parallel, furrr::future_map, purrr::map)
  results <- x.list %>%
    set_names(sapply(x.list, function(df) df[[1,1]])) %>%
    fun_map(
      ~{
        item <- .
        x.raw <- item[["value_raw"]]
        unit <- item[[1, 1]]

        res <- TFDTW(x = item[["value"]],
                     y = y.processed,
                     k = k,
                     t.treat = t.treat,
                     buffer = buffer,
                     norm.method = norm.method,
                     match.method = match.method,
                     step.pattern1 = step.pattern1,
                     step.pattern2 = step.pattern2,
                     plot.figures = plot.figures,
                     n.burn = n.burn,
                     ma = ma, ma.na = ma.na,
                     dist.quant = dist.quant, n.IQR = n.IQR,
                     window.type = window.type,
                     default.margin = default.margin,
                     n.q = n.q, n.r = n.r)

        x.warped <- c(
          warpWITHweight(x.raw[1:res$cutoff], res$weight.a)[1:t.treat],
          warpWITHweight(x.raw[res$cutoff:length(x.raw)], res$avg.weight)[-1]
        )

        df.synth <- data.frame(
          time = 1:length(y.raw) + start.time - 1,
          unit = unit,
          value_warped = NA,
          stringsAsFactors = FALSE
        )
        df.synth$value_warped <- x.warped[1:length(y.raw)]

        list(unit = unit,
             x.warped = x.warped,
             df.synth = df.synth)
      }
    )

  # Prepare data for Synth
  df.synth <- lapply(results, "[[", "df.synth") %>%
    bind_rows() %>%
    add_row(time = seq(start.time, length(y.raw) + start.time - 1),
            unit = dependent,
            value_warped = y.raw)

  df.synth = right_join(data, df.synth, by = c("unit", "time"))
  dependent.id = df.synth$id[which.max(df.synth$unit == dependent)]

  # Synth
  res.synth = do.synth(df = df.synth,
                        dep.var = "value_warped",
                        dependent.id = dependent.id,
                        predictors = predictors,
                        special.predictors = special.predictors,
                        time.predictors.prior = time.predictors.prior,
                        time.optimize.ssr = time.optimize.ssr)

  return(res.synth)
}

