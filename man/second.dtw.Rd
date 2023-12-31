% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TFDTW.R
\name{second.dtw}
\alias{second.dtw}
\title{Dynamic Time Warping: Second Phase}
\usage{
second.dtw(
  x.post,
  x.pre,
  k,
  weight.a,
  norm.method = "t",
  default.margin = 3,
  n.q = 1,
  n.r = 1,
  step.pattern = dtw::asymmetricP2,
  window.type = "none",
  window.size = NULL,
  dist.quant = 1,
  n.IQR = 3,
  ...
)
}
\arguments{
\item{x.post}{Post values of x.}

\item{x.pre}{Pre values of x.}

\item{k}{An integer.}

\item{weight.a}{Weight values for a.}

\item{norm.method}{Normalization method. Default is "t".}

\item{default.margin}{Default margin size. Default is 3.}

\item{n.q}{An integer, default 1.}

\item{n.r}{An integer, default 1.}

\item{step.pattern}{Step pattern to use. Default is \code{dtw::asymmetricP2}.}

\item{window.type}{Type of window. Default is "none".}

\item{window.size}{Optional size for the window. Default is NULL.}

\item{dist.quant}{Distance quantile. Default is 1.}

\item{n.IQR}{Number of IQR. Default is 3.}

\item{...}{Additional arguments passed to other methods.}
}
\value{
A list containing various outputs from the second DTW process.
}
\description{
This function calculates the dynamic time warping for the second phase.
}
\keyword{internal}
