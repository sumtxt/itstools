#' Plot the Estimation and Placebo Sample
#'
#' \code{its_llm_plot} helps to visualize which data points are used to estimate the intercept shift 
#' at the cut-point and the placebo estimates.  
#'
#' @param df (required) \code{data.frame} containing all variables 
#' @param rvar (required) the name of the running variable in \code{df}
#' @param outcome (required) the name of the outcome variable in \code{df} 
#' @param trend include a linear term ('lin'), a quadratic term ('quad') or no trend at all ('none')?
#' @param bw either a scalar or a vector of length 2 defining the bandwidth to the left (right) of the cut-point on the scale of \code{rvar}
#' @param donut either a scalar or a vector of length 2 defining the length of the period to the left (right) of the cut-point for which the data are dropped (on the scale of \code{rvar}).
#' 
#' 
#' 
#' @examples 
#'  \dontrun{
#'  
#'   N <- 21
#'   time <- seq(-1,1,length.out=N)
#'   treat <- as.numeric(time >= 0)
#'   y <- 1 + time + treat*1 + rnorm(N,0,0.25)
#'   df <- data.frame(y=y, time=time)
#'   its_llm_plot(df, rvar="time", outcome="y", bw=0.25)
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
its_llm_plot <- function(df, rvar, outcome, bw, donut=0){
	list[bwL,bwR] <- parseLR(bw)
	list[donutL,donutR] <- parseLR(donut)
	df <- as.data.frame(df)
	df <- make_inSampleVar(df,rvar=rvar, bwL=bwL, bwR=bwR, donutL=donutL, donutR=donutR)
	sset <- df[df[,'inSample']==TRUE,rvar]
	pset <- make_permut_set(df, rvar=rvar, bwL=bwL, bwR=bwR, donutL=donutL, donutR=donutR)
	plot(df[,rvar],df[,outcome], type='n', xlab=rvar, ylab=outcome)
	abline(v=0, col='blue')
	abline(v=bwR+donutR,lty=2, col='red')
	abline(v=donutR,lty=2, col='blue') 
	abline(v=-(bwL+donutL),lty=2, col='red')
	abline(v=-donutL,lty=2, col='blue') 
	points(df[,rvar],df[,outcome], , pch=20, col='grey')
	points(df[df[,rvar] %in% pset,rvar], 
		df[df[,rvar] %in% pset,outcome], col='blue', pch=4)
	points(df[df[,rvar] %in% sset,rvar], 
		df[df[,rvar] %in% sset,outcome], col='green', pch=8)
	legend('topleft', c("Estimation Sample","Placebo Events"), col=c("green", 'blue'), 
			pch=c(4,8), bty = "n")
	legend('bottomright', c("Bandwidth", "Donut"), col=c('red', 'blue'), lty=c(2,2), bty='n')
	}

