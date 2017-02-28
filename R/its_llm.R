#' Estimates via Local Linear Regression with user-selected bandwidth
#'
#' \code{its_llm} estimates the intercept shift of a time series at a cut-point. 
#'
#'  
#' @param df (required) \code{data.frame} containing all variables 
#' @param rvar (required) the name of the running variable in \code{df}
#' @param outcome (required) the name of the outcome variable in \code{df} 
#' @param trend include a linear term ('lin'), a quadratic term ('quad') or no trend at all ('none')?
#' @param bw either a scalar or a vector of length 2 defining the bandwidth to the left (right) of the cut-point on the scale of \code{rvar}
#' @param donut either a scalar or a vector of length 2 defining the length of the period to the left (right) of the cut-point for which the data are dropped (on the scale of \code{rvar}).
#' @param verbose set to any value other than zero to show which data points are included in the estimation
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
#' \dontrun{
#'		
#'   N <- 21
#'   time <- seq(-1,1,length.out=N)
#'   treat <- as.numeric(time >= 0)
#'   y <- 1 + time + treat*1 + rnorm(N,0,0.25)
#'   df <- data.frame(y=y, time=time)
#'   its_llm(df, rvar="time", outcome="y", bw=0.25)
#' 
#' }
#' 
#' @return \code{data.frame} with a single row and entries for the point estimate (\code{est}), 95\% confidence interval (\code{lo,hi}), 
#' standard error (\code{se}), p-value (\code{pval}) and the number of data points to the left/right of the cut-point used in the 
#' estimation (\code{Nleft,Nright}).
#' 
#' @seealso \code{\link{its_plot_samples}}, \code{\link{its_llm_placebo}}.
#' 
#' @export
its_llm <- function(df, rvar, outcome, trend="none", bw=Inf, donut=0, verbose=0){
	list[bwL,bwR] <- parseLR(bw)
	list[donutL,donutR] <- parseLR(donut)
	df <- as.data.frame(df)
	df <- df[order(df[,rvar]),]
	df <- make_inSampleVar(df,rvar=rvar, bwL=bwL, bwR=bwR, donutL=donutL, donutR=donutR)
	df$treat <- as.numeric(df[,rvar]>=0) 
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
	if (verbose!=0){
		df_ <- df[df[,'inSample']==TRUE,]
		cat("\nData Points: \n")
		cat("to the left of 0:\n")
		cat("- smallest value:", min(df_[df_[,rvar]<0,rvar]), "\n") 
		cat("- largest value:", max(df_[df_[,rvar]<0,rvar]), "\n" ) 
		cat("to the right of 0:\n")
		cat("- smallest value:", min(df_[df_[,rvar]>=0,rvar]), "\n") 
		cat("- largest value:", max(df_[df_[,rvar]>=0,rvar]), "\n\n" ) 
	}
	return(m)
	} 


