#' Plot the Estimation and Placebo Sample
#'
#' \code{its_llm_plot} helps to visualize which data points are used to estimate the intercept shift 
#' at the cut-point and the placebo estimates.  
#'
#' @param df (required) \code{data.frame} containing all variables 
#' @param rvar (required) the name of the running variable in \code{df}
#' @param outcome (required) the name of the outcome variable in \code{df} 
#' @param trend include a linear term ('lin'), a quadratic term ('quad') or no trend at all ('none')?
#' @param bwL bandwidth selection to the left of the cut-point on the scale of \code{rvar}
#' @param bwR bandwidth selection to the right of the cut-point on the scale of \code{rvar}
#' @param donut length of the period after the cut-point for which the data are dropped (on the scale of \code{rvar}).
#' 
#' 
#' 
#' @examples 
#'  \dontrun{
#' 
#'   eventmonth <- -12:12
#'   treat <- as.numeric(eventmonth >= 0)
#'   y <- 1 + (1 * eventmonth) + treat*2 + rnorm(length(eventmonth))
#'   df <- data.frame(y=y, eventmonth=eventmonth)
#'   its_llm_plot(df, rvar="eventmonth", outcome="y", bwL=2,bwR=2, donut=0)
#' 
#' }
#' 
#' @references 
#'	
#' 
#' 
#' @seealso \code{\link{its_llm_placebo}}, \code{\link{its_llm}}.
#' 
#' @export
its_llm_plot <- function(df, rvar, outcome, bwL, bwR, donut){
	df <- as.data.frame(df)
	pset <- make_permut_set(df, rvar, bwL, bwR, donut)
	df_ <- make_donut_df(df, rvar=rvar, donut=donut)
	sset <- df_[ -bwL <= df_[,rvar] & df_[,rvar] <= (donut+bwR),rvar]
	plot(df[,rvar],df[,outcome], type='n', xlab=rvar, ylab=outcome)
	abline(v=0, col='blue')
	abline(v=bwR+donut,lty=2, col='red')
	abline(v=donut,lty=2, col='blue') 
	abline(v=-bwL,lty=2, col='red')
	points(df[,rvar],df[,outcome], , pch=20, col='grey')
	points(df[df[,rvar] %in% pset,rvar], 
		df[df[,rvar] %in% pset,outcome], col='blue', pch=4)
	points(df[df[,rvar] %in% sset,rvar], 
		df[df[,rvar] %in% sset,outcome], col='green', pch=8)
	legend('topleft', c("Estimation Sample","Placebo Events"), col=c("green", 'blue'), 
			pch=c(4,8), bty = "n")
	}
