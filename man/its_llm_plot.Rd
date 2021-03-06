% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/its_llm_plot.R
\name{its_llm_plot}
\alias{its_llm_plot}
\title{Plot the Estimation and Placebo Sample}
\usage{
its_llm_plot(df, rvar, outcome, bw, donut = 0)
}
\arguments{
\item{df}{(required) \code{data.frame} containing all variables}

\item{rvar}{(required) the name of the running variable in \code{df}}

\item{outcome}{(required) the name of the outcome variable in \code{df}}

\item{bw}{either a scalar or a vector of length 2 defining the bandwidth to the left (right) of the cut-point on the scale of \code{rvar}}

\item{donut}{either a scalar or a vector of length 2 defining the length of the period to the left (right) of the cut-point for which the data are dropped (on the scale of \code{rvar}).}

\item{trend}{include a linear term ('lin'), a quadratic term ('quad') or no trend at all ('none')?}
}
\description{
\code{its_llm_plot} helps to visualize which data points are used to estimate the intercept shift 
at the cut-point and the placebo estimates.
}
\examples{
 \dontrun{
 
  N <- 21
  time <- seq(-1,1,length.out=N)
  treat <- as.numeric(time >= 0)
  y <- 1 + time + treat*1 + rnorm(N,0,0.25)
  df <- data.frame(y=y, time=time)
  its_llm_plot(df, rvar="time", outcome="y", bw=0.25)

}

}
\references{

}
\seealso{
\code{\link{its_llm_placebo}}, \code{\link{its_llm}}.
}

