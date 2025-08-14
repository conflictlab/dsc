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
#' @importFrom dplyr group_by summarize ungroup mutate select left_join
#' @importFrom dplyr right_join group_split bind_rows mutate_at mutate_all
#' @importFrom ggplot2 ggplot aes geom_line theme_bw
#' @importFrom magrittr `%>%`
#' @importFrom stats quantile sd
#' @importFrom tibble add_row
#' @importFrom purrr set_names
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' library(dsc)
#'
#' # Load the Basque dataset from the Synth package
#' data(basque, package = "Synth")
#' data <- basque
#'
#' # Rename relevant columns for clarity
#' colnames(data)[1:4] <- c("id", "unit", "time", "value")
#'
#' # Compute additional variables
#' data$invest_ratio <- data$invest / data$value
#' data$value_raw <- data$value
#'
#' # Define special predictors for the model
#' special_preds <- expression(list(
#'   list(dep.var, 1960:1969, c("mean")),
#'   list("invest_ratio", 1964:1969, c("mean")),
#'   list("popdens", 1969, c("mean")),
#'   list("sec.agriculture", 1961:1969, c("mean")),
#'   list("sec.energy", 1961:1969, c("mean")),
#'   list("sec.industry", 1961:1969, c("mean")),
#'   list("sec.construction", 1961:1969, c("mean")),
#'   list("sec.services.venta", 1961:1969, c("mean")),
#'   list("sec.services.nonventa", 1961:1969, c("mean")),
#'   list("school.illit", 1964:1969, c("mean")),
#'   list("school.prim", 1964:1969, c("mean")),
#'   list("school.med", 1964:1969, c("mean")),
#'   list("school.high", 1964:1969, c("mean")),
#'   list("school.post.high", 1964:1969, c("mean"))
#' ))
#'
#' # Execute the DSC analysis
#' result <- dsc(
#'   data = data,
#'   start.time = 1955,
#'   end.time = 1997,
#'   treat.time = 1970,
#'   dependent = "Basque Country (Pais Vasco)",
#'   predictors = NULL,
#'   parallel = TRUE,
#'   special.predictors = special_preds,
#'   time.predictors.prior = 1955:1969,
#'   time.optimize.ssr = 1955:1969
#' )
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
      dplyr::filter(.data$time <= treat.time) %>%
      group_by(.data$unit) %>%
      summarize(
        value.min = min(.data$value),
        value.max = max(.data$value)) %>%
      ungroup()

    mean.diff = mean(rescale_data$value.max - rescale_data$value.min)
    rescale_data = rescale_data %>%
      mutate(
        multiplier = mean.diff/(.data$value.max - .data$value.min)
      )

    data <- left_join(data, rescale_data, by = "unit") %>%
      mutate(
        value.bak = .data$value_raw,
        value_raw = (.data$value_raw - .data$value.min) * .data$multiplier,
        value = .data$value_raw
      )
  }

  # Filter data
  if (!is.null(filter.width)) {
    data = preprocessing(data, filter.width)
  }

  # Prepare data for TFDTW
  t.treat <- (treat.time - start.time) + 1

  y_values <- data %>%
    dplyr::filter(.data$unit == dependent & .data$time <= end.time) %>%
    select(.data$value_raw, .data$value)

  y.raw <- y_values[["value_raw"]]
  y.processed <- y_values[["value"]]

  x.list <- data %>%
    dplyr::filter(.data$unit != dependent & .data$time <= end.time) %>%
    select(.data$unit, .data$time, .data$value, .data$value_raw) %>%
    group_by(.data$unit) %>%
    group_split(.keep = TRUE)

  # TFDTW
  if (parallel) {
    fun_map <- furrr::future_map
    future::plan(future::multisession, workers = parallel::detectCores() - 1)
  }else{
    fun_map <- purrr::map
  }
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
  dependent.id = df.synth[["id"]][which.max(df.synth$unit == dependent)]

  # Synth
  res.synth = do.synth(df = df.synth,
                        dep.var = "value_warped",
                        dependent.id = dependent.id,
                        predictors = predictors,
                        special.predictors = special.predictors,
                        time.predictors.prior = time.predictors.prior,
                        time.optimize.ssr = time.optimize.ssr)
  
  # Plot figure
  if (plot.figures) {
    # Prepare data frame
    df.plot <- data.frame(
      value = c(res.synth$value, res.synth$synthetic),
      series = c(rep("Treated Unit", length(res.synth$value)),
                 rep("DSC", length(res.synth$value))),
      x = rep(start.time:end.time, 2)
    )
    
    # Plot with ggplot
    ggplot(df.plot, aes(x = x, y = value, color = series)) +
      geom_line(size = 1) +
      geom_vline(xintercept = treat.time, linetype = "dashed") +
      labs(
        title = dependent,
        x = "Time",
        y = "Value",
        color = NULL
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
  }

  return(res.synth)
}

