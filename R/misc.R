#' @title Standard Z-score Normalization
#' @description Normalizes the data using z-score.
#' @param data A numeric vector to be normalized.
#' @param reference An optional reference for normalization.
#' @return Normalized data.
t.normalize <- function(data, reference = NULL) {
  if (is.null(reference)) {
    mu = mean(data)
    sigma = sd(data)
  }else{
    mu = mean(reference)
    sigma = sd(reference)
  }
  res = (data - mu)/sigma
  return(res)
}

#' @title Min-Max Normalization
#' @description Normalizes the data using min-max scaling.
#' @param data A numeric vector to be normalized.
#' @param reference An optional reference for normalization.
#' @return Normalized data.
minmax.normalize <- function(data, reference = NULL) {
  if (is.null(reference)) {
    minimum = min(data)
    maximum = max(data)
  }else{
    minimum = min(reference)
    maximum = max(reference)
  }
  res = (data - minimum)/(maximum - minimum)
  return(res)
}

#' @title General Normalization Function
#' @description Applies the specified normalization to the data.
#' @param data A numeric vector to be normalized.
#' @param norm.method A string that specifies the normalization method.
#' @param reference An optional reference for normalization.
#' @return Normalized data.
normalize <- function(data, norm.method, reference = NULL) {
  if (norm.method == "minmax") {
    data = minmax.normalize(data, reference)
  }else if (norm.method == "t") {
    data = t.normalize(data, reference)
  }
  return(data)
}

#' @title Add Buffer to Time Series
#' @description Adds a buffer to a time series using auto.arima.
#' @param TS Time series data.
#' @param n Buffer size.
#' @return Time series with added buffer.
add.buffer <- function(TS, n) {
  model.right = forecast::auto.arima(TS)
  right <- as.numeric(forecast::forecast(model.right, h = n)$mean)
  model.left = forecast::auto.arima(rev(TS))
  left <- rev(as.numeric(forecast::forecast(model.left, h = n)$mean))

  return(c(left, TS, right))
}

#' @title Pre-process Data
#' @description Applies several pre-processing steps to the data.
#' @param data The data to be processed.
#' @param filter.width Width of the filter.
#' @param norm.method Normalization method.
#' @param n.poly Degree of the polynomial.
#' @param n.deri Order of the derivative.
#' @param plot.data Logical, indicating whether to plot the data.
#' @return Processed data.
preprocessing <- function(data, filter.width = 5, norm.method = "t",
                          n.poly = 3, n.deri = 2, plot.data = FALSE) {
  # transform
  values = reshape2::dcast(data %>% select(c("unit", "time", "value_raw")),
                           time ~ unit, value.var = "value_raw")

  # normalize
  values = values %>% mutate_at(setdiff(colnames(values), "time"),
                                ~normalize(., norm.method))

  # add buffer
  n.buffer = (filter.width - 1)/2
  values.w.buffer = sapply(values %>% select(-.data$time),
                           add.buffer, n = n.buffer) %>%
    data.frame()

  # derivative
  values.w.buffer = values.w.buffer %>%
    mutate_all(~signal::sgolayfilt(., n.poly, filter.width, n.deri))
  values[-1] = values.w.buffer[(n.buffer + 1):(n.buffer + nrow(values)),]

  # join
  df <- reshape2::melt(values, id.vars = 'time', variable.name = 'unit')
  data = right_join(df, data %>% select(-.data$value), by = c("time", "unit"))
  data = data[c("id", "unit", "time", "value", colnames(data)[-(1:4)])]

  # plot
  if (plot.data) {
    ggplot(data, aes(x = .data$time, y = .data$value, color = .data$unit)) +
      geom_line() +
      theme_bw()
  }

  return(data)
}

#' @title Convert Warping Path to Weight
#' @description Transforms a warping path to weights.
#' @param W The warping path.
#' @return Weights.
warp2weight <- function(W) {
  w = as.matrix(W)
  count = rep(1/colSums(w), nrow(w)) %>%
    matrix(.,
           nrow = ncol(w),
           ncol = nrow(w)) %>%
    t()
  weight = rowSums(w * count)

  return(weight)
}

#' @title Warp Time Series with Weights
#' @description Warps a time series using the provided weights.
#' @param ts Time series to be warped.
#' @param weight Weights for warping.
#' @return Warped time series.
warpWITHweight <- function(ts, weight) {
  n.ts = length(ts)
  n.weight = length(weight)
  if (n.ts != n.weight) {
    print("Lengths of ts and weight are different.")
    return()
  }

  # form tranformation matrix
  count = 1
  residual = 0
  w = matrix(0, nrow = n.weight, ncol = 2*n.weight)
  for (i in 1:n.weight) {
    z = weight[i]
    while (z > 0) {
      if (z + residual >= 1) {
        now = 1 - residual
        residual = 0
        z = z - now
        w[i, count] = now
        count = count + 1
      }else{
        residual = residual + z
        now = z
        z = 0
        w[i, count] = now
      }
    }
  }

  col.sum = colSums(w)
  ind.max = which(col.sum > 0) %>%
    max
  w = w[, 1:ind.max]
  w[n.weight, ind.max] = 1 - sum(w[-n.weight, ind.max])

  # warp time series
  ts.warped = t(w) %*% matrix(ts, ncol = 1)

  return(ts.warped[,1])
}

#' @title Check if Reference is Short in DTW
#' @description Checks if the reference series is too short for dynamic time warping.
#' @param query The query series.
#' @param reference The reference series.
#' @param step.pattern Step pattern for DTW.
#' @param window.type Type of window for DTW.
#' @param window.size Size of the window for DTW.
#' @return Logical indicating if reference is too short.
RefTooShort <- function(query, reference,
                        step.pattern = dtw::symmetricP2,
                        window.type = "none",
                        window.size = NULL) {
  alignment = tryCatch(dtw::dtw(reference, query,
                                step.pattern = step.pattern,
                                window.type = window.type,
                                window.size = window.size,
                                open.end = TRUE),
                       error = function(e) return(NULL))
  if (is.null(alignment)) {
    return(FALSE)
  }
  if (length(unique(alignment$index2)) == 1) {
    return(TRUE)
  }
  wq = suppressWarnings(dtw::warp(alignment, index.reference = FALSE))
  return(round(max(wq)) < length(query))
}

#' @title Remove Outliers in Weight Matrix
#' @description Removes outliers from the data based on interquartile range.
#' @param data Numeric data.
#' @param n.IQR Multiplier for the interquartile range.
#' @return Data with outliers removed.
RemoveOutliers <- function(data, n.IQR = 3) {
  Q1 = stats::quantile(data, 0.25, na.rm = TRUE)
  Q3 = stats::quantile(data, 0.75, na.rm = TRUE)
  IQR = Q3 - Q1
  upper = Q3 + n.IQR*IQR
  lower = Q1 - n.IQR*IQR
  data[data > upper] = NaN
  data[data < lower] = NaN
  return(data)
}

