#' Estimates via Local Linear Regression with user-selected bandwidth
#'
#' \code{its_llm} estimates the intercept shift of a time series at a cut-point. 
#'
#'  
#' @param df (required) \code{data.frame} containing all variables 
#' @param rvar (required) the name of the running variable in \code{df}
#' @param outcome (required) the name of the outcome variable in \code{df} 
#' @param trend include a linear term ('lin'), a quadratic term ('quad') or no trend at all ('none')?
#' @param bwL bandwidth selection to the left of the cut-point on the scale of \code{rvar}
#' @param bwR bandwidth selection to the right of the cut-point on the scale of \code{rvar}
#' @param donut length of the period after the cut-point for which the data are dropped (on the scale of \code{rvar}).
#' 
#' @details 
#' Estimates the size of the intercept shift of a time series at at cut-point (at zero) using 
#' a linear regression model with separate trends for the running variable to both sides
#' of the cut-point and within the neighborhood as defined by the bandwidth parameters (\code{bwL},\code{bwR}).
#' Standard errors are calculated based on the heteroskedasticity-consistent covariance matrix (HC3) from the \link{sandwich} package.
#' 
#' Use \code{\link{its_plot_samples}} to understand which data points are included when choosing different values for \code{bwL}, 
#' \code{bwR} and \code{donut}. 
#' 
#' When no values for \code{trend}, \code{bwL}, \code{bwR} and \code{donut} are supplied, the functions defaults to estimating the 
#' difference in means pooling all available data to the left/right of the cut-point. 
#' 
#' 
#' @examples 
#'  \dontrun{
#'		
#'   	eventmonth <- -12:12
#'   	treat <- as.numeric(eventmonth >= 0)
#'   	y <- 1 + (1 * eventmonth) + treat*2 + rnorm(length(eventmonth))
#'   	df <- data.frame(y=y, eventmonth=eventmonth)
#'   	its_llm(df, rvar="eventmonth", outcome="y", trend='lin', bwL=4,bwR=4, donut=0) 		
#' 
#' 		}
#' 
#'	@return \code{data.frame} with a single row and entries for the point estimate (\code{est}), 95\% confidence interval (\code{lo,hi}), 
#'	standard error (\code{se}), p-value (\code{pval}) and the number of data points to the left/right of the cut-point used in the 
#'  estimation (\code{Nleft,Nright}).
#' 
#' @seealso \code{\link{its_plot_samples}}, \code{\link{its_llm_placebo}}.
#' 
#' @export
its_llm <- function(df, rvar, outcome, trend="none", bwL=Inf, bwR=Inf, donut=0){
	df <- as.data.frame(df)
	df <- df[order(df[,rvar]),]
	make_donut_df(df, rvar=rvar, donut=donut)
	df$treat <- as.numeric(df[,rvar]>=0) 
	df$inSample <- (-bwL <= df[,rvar]) & (df[,rvar] <= (donut+bwR))
	if ( !(trend %in% c("none", "lin", "quad")) ) stop("Parameter 'trend' must be either 'none', 'lin' or 'quad'.")
	if ( trend=="none" ) spec <- paste(outcome, "~ treat",sep="")
	if ( trend=="lin" ) spec <- paste(outcome, " ~ treat * ",rvar, sep="")
	if ( trend=="quad" ) spec <- paste(outcome, " ~ treat * (", rvar, " + I(", rvar, "^2))", sep="")
	NL <- sum(df[df[,'inSample']==TRUE,'treat']==0)
	NR <- sum(df[df[,'inSample']==TRUE,'treat']==1)
	if (NL==0 | NR==0) return(data.frame(est=NA, lo=NA, hi=NA, se=NA, pval=NA, Nleft=NL, Nright=NR))
	m <- lmHC(formula=as.formula(spec), data=df[df[,'inSample']==TRUE,])
	m$Nleft <- NL 
	m$Nright <- NR
	return(m)
	}
