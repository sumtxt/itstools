#' ITS Placebo Estimates 
#'
#' \code{its_llm_placebo} gives a vector of placebo estimates
#'
#'  
#' @param df (required) \code{data.frame} containing all variables 
#' @param rvar (required) a running variable in \code{df}
#' @param outcome (required) the outcome variable in \code{df}
#' @param trend can be 'lin' or 'quad' (linear or quadratic)
#' @param bw bandwidth selection on the scale of the \code{rvar}
#' @param donut length of the period in after the cutpoint that is not used for the estimation 
#' @param nsim number of simulations to run 
#' 
#' @details 
#' tbd.
#' 
#' @examples 
#'  \dontrun{
#'
#' }
#' 
#' @references 
#'	
#' 
#' 
#' @export
its_llm_placebo <- function(df, rvar, outcome, 
	trend=NULL, bw, donut=0, nsim=200){
	pset <- make_permut_set(df, rvar, bw,donut)
	if ( length(pset) > nsim ) pset <- sample(pset, nsim)
	m <- sapply(pset, function(value){
		dat <- df
		dat[,rvar] <- (dat[,rvar] - value)
		fit <- its_llm(dat, rvar=rvar, 
			outcome=outcome, trend=trend, 
			bw=bw, donut=donut)
		return(fit['est'])
		})
	return(as.vector(do.call(c, m)))
	}
