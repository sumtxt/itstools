#' Placebo Estimates via Local Linear Regression with user-selected bandwidth
#'
#' \code{its_llm} estimates placebo intercept shift of a time series permuting the cut-point.
#'
#'  
#' @param df (required) \code{data.frame} containing all variables 
#' @param rvar (required) the name of the running variable in \code{df}
#' @param outcome (required) the name of the outcome variable in \code{df} 
#' @param trend include a linear term ('lin'), a quadratic term ('quad') or no trend at all ('none')?
#' @param bw either a scalar or a vector of length 2 defining the bandwidth to the left (right) of the cut-point on the scale of \code{rvar}
#' @param donut either a scalar or a vector of length 2 defining the length of the period to the left (right) of the cut-point for which the data are dropped (on the scale of \code{rvar}).
#' @param nsim if the number of potential placebo estimates is larger than \code{nsim}, 
#' 	only a random sample of \code{nsim} estimates is formed. 
#'
#' 
#' @details 
#' Permutes the cut-point of a time series within the data points that are observed before the cut-point. The resulting 
#' distribution can be thought of as a distribution of typical shifts in the time series conditional on a selected  
#' local linear regression specification and might used as a reference distribution to make inferences about the actual estimate
#' at the cut-point.  
#' 
#' The function only forms placebo estimates for values of the running variable that are sufficiently enough away 
#' from the left boundary of the time series (\code{min(rvar)+bwL}) and the cut-point (\code{0-bwR-donut})
#' such that a local linear regression with the specified bandwidth can be estimated. Use \code{\link{its_plot_samples}} to 
#' visualize the implied interval from which placebo estimates are formed. 
#' 
#' @examples 
#' \dontrun{
#'
#'   N <- 21
#'   time <- seq(-1,1,length.out=N)
#'   treat <- as.numeric(time >= 0)
#'   y <- 1 + time + treat*1 + rnorm(N,0,0.25)
#'   df <- data.frame(y=y, time=time)
#'   its_llm_placebo(df, rvar="time", outcome="y", bw=0.25)
#' 
#' }
#' 
#' @return \code{numeric} vector of placebo estimates. 
#' 
#' @seealso \code{\link{its_plot_samples}}, \code{\link{its_llm}}.
#' 
#' @export
its_llm_placebo <- function(df, rvar, outcome, 
	trend="none", bw, donut=0, nsim=200){
	list[bwL,bwR] <- parseLR(bw)
	list[donutL,donutR] <- parseLR(donut)
	df <- as.data.frame(df)
	pset <- make_permut_set(df, rvar=rvar, bwL=bwL, bwR=bwR, donutL=donutL, donutR=donutR)
	pset <- sort(unique(pset))
	if ( length(pset) > nsim ) pset <- sample(pset, nsim)
	m <- sapply(pset, function(value){
		dat <- df
		dat[,rvar] <- (dat[,rvar] - value)
		fit <- its_llm(dat, rvar=rvar, 
			outcome=outcome, trend=trend, 
			bw=c(bwL,bwR), donut=c(donutL,donutR) )
		return(fit['est'])
		})
	return(as.vector(do.call(c, m)))
	}

