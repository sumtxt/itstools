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
#'   	eventmonth <- -12:12
#'   	treat <- as.numeric(eventmonth >= 0)
#'   	y <- 1 + (1 * eventmonth) + treat*2 + rnorm(length(eventmonth))
#'   	df <- data.frame(y=y, eventmonth=eventmonth)
#'   	its_llm_placebo(df, rvar="eventmonth", outcome="y", bwL=2,bwR=2, donut=0) 		
#' 
#' }
#' 
#' @references 
#'	
#' 
#' 
#' @export
its_llm_placebo <- function(df, rvar, outcome, 
	trend="none", bwL, bwR, donut=0, nsim=200){
	df <- as.data.frame(df)
	pset <- sort(unique(make_permut_set(df, rvar, bwL, bwR, donut)))
	if ( length(pset) > nsim ) pset <- sample(pset, nsim)
	m <- sapply(pset, function(value){
		dat <- df
		dat[,rvar] <- (dat[,rvar] - value)
		fit <- its_llm(dat, rvar=rvar, 
			outcome=outcome, trend=trend, 
			bwL=bwL, bwR=bwR, donut=donut)
		return(fit['est'])
		})
	return(as.vector(do.call(c, m)))
	}

