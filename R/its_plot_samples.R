#' Plots ITS Samples 
#'
#' \code{its_plot_samples} ... 
#'
#'  
#' @param df (required) \code{data.frame} containing all variables 
#' @param rvar (required) a running variable in \code{df}
#' @param outcome (required) the outcome variable in \code{df}
#' @param bw bandwidth selection on the scale of the \code{rvar}
#' @param donut length of the period in after the cutpoint that is not used for the estimation 
#' 
#' @details 
#' tbd.
#' 
#' @examples 
#'  \dontrun{
#' 
#'   eventmonth <- -12:12
#'   treat <- as.numeric(eventmonth >= 0)
#'   y <- 1 + (1 * eventmonth) + treat*2 + rnorm(length(eventmonth))
#'   df <- data.frame(y=y, eventmonth=eventmonth)
#'   its_plot_samples(df, rvar="eventmonth", outcome="y", bwL=2,bwR=2, donut=0)
#' 
#' }
#' 
#' @references 
#'	
#' 
#' 
#' @export
its_plot_samples <- function(df, rvar, outcome, bwL, bwR, donut){
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