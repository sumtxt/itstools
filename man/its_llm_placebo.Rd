% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/its_llm_placebo.R
\name{its_llm_placebo}
\alias{its_llm_placebo}
\title{Placebo Estimates via Local Linear Regression with user-selected bandwidth}
\usage{
its_llm_placebo(df, rvar, outcome, trend = "none", bw, donut = 0,
  nsim = 200)
}
\arguments{
\item{df}{(required) \code{data.frame} containing all variables}

\item{rvar}{(required) the name of the running variable in \code{df}}

\item{outcome}{(required) the name of the outcome variable in \code{df}}

\item{trend}{include a linear term ('lin'), a quadratic term ('quad') or no trend at all ('none')?}

\item{bw}{either a scalar or a vector of length 2 defining the bandwidth to the left (right) of the cut-point on the scale of \code{rvar}}

\item{donut}{either a scalar or a vector of length 2 defining the length of the period to the left (right) of the cut-point for which the data are dropped (on the scale of \code{rvar}).}

\item{nsim}{if the number of potential placebo estimates is larger than \code{nsim}, 
only a random sample of \code{nsim} estimates is formed.}
}
\value{
\code{numeric} vector of placebo estimates.
}
\description{
\code{its_llm} estimates placebo intercept shift of a time series permuting the cut-point.
}
\details{
Permutes the cut-point of a time series within the data points that are observed before the cut-point. The resulting 
distribution can be thought of as a distribution of typical shifts in the time series conditional on a selected  
local linear regression specification and might used as a reference distribution to make inferences about the actual estimate
at the cut-point.  

The function only forms placebo estimates for values of the running variable that are sufficiently enough away 
from the left boundary of the time series (\code{min(rvar)+bwL}) and the cut-point (\code{0-bwR-donut})
such that a local linear regression with the specified bandwidth can be estimated. Use \code{\link{its_plot_samples}} to 
visualize the implied interval from which placebo estimates are formed.
}
\examples{
\dontrun{

  N <- 21
  time <- seq(-1,1,length.out=N)
  treat <- as.numeric(time >= 0)
  y <- 1 + time + treat*1 + rnorm(N,0,0.25)
  df <- data.frame(y=y, time=time)
  its_llm_placebo(df, rvar="time", outcome="y", bw=0.25)

}

}
\seealso{
\code{\link{its_plot_samples}}, \code{\link{its_llm}}.
}

