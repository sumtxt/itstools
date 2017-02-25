#' ITS Estimates via Local Linear Regression 
#'
#' \code{its_llm} Estimates change in time series via local linear regression 
#'
#'  
#' @param df (required) \code{data.frame} containing all variables 
#' @param rvar (required) a running variable in \code{df}
#' @param outcome (required) the outcome variable in \code{df}
#' @param trend can be 'lin' or 'quad' (linear or quadratic)
#' @param bw bandwidth selection on the scale of the \code{rvar}
#' @param donut length of the period in after the cutpoint that is not used for the estimation 
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
#'   	its_llm(df, rvar="eventmonth", outcome="y", bwL=2,bwR=2, donut=0) 		
#' 
#' 		}
#' 
#' @references 
#'	
#' 
#' 
#' @export
its_llm <- function(df, rvar, outcome, trend="none", bwL=Inf, bwR=Inf, donut=0){
	df <- as.data.frame(df)
	df <- df[order(df[,rvar]),]
	make_donut_df(df, rvar=rvar, donut=donut)
	df$treat <- as.numeric(df[,rvar]>=0) 
	df$inSample <- (-bwL <= df[,rvar]) & (df[,rvar] <= (donut+bwR))
	spec <- paste(outcome, "~ treat",sep="")
	if ( trend=="lin" ) spec <- paste(outcome, " ~ treat * ",rvar, sep="")
	if ( trend=="quad" ) spec <- paste(outcome, " ~ treat * (", rvar, " + I(", rvar, ")^2)", sep="")
	NL <- sum(df[df[,'inSample']==TRUE,'treat']==0)
	NR <- sum(df[df[,'inSample']==TRUE,'treat']==1)
	if (NL==0 | NR==0) return(data.frame(est=NA, lo=NA, hi=NA, se=NA, pval=NA, Nleft=NL, Nright=NR))
	m <- lmHC(formula=as.formula(spec), data=df[df[,'inSample']==TRUE,])
	m$Nleft <- NL 
	m$Nright <- NR
	return(m)
	}
