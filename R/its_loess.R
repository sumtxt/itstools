#' ITS Estimates via Local Polynomial Regression 
#'
#' \code{its_llm} Estimates change in time series via local linear regression 
#'
#'  
#' @param df (required) \code{data.frame} containing all variables 
#' @param rvar (required) a running variable in \code{df}
#' @param outcome (required) the outcome variable in \code{df}
#' @param span the span passed to \code{loess}
#' @param degree the degree of the polynomial passed to \code{loess}
#' @param othervar vector of variable names to be passed into the returned \code{data.frame}
#' @param only_estimate if \code{TRUE} only returns the estimates at the cutpoint (default)
#' 
#' @details 
#' This function estimates two loess functions to the right and left of the cutoff and uses the predicted 
#' values from the fitted functions to form an estimate about the shift. 
#' 
#' @examples 
#'  \dontrun{
#'		
#'   	eventmonth <- -12:12
#'   	treat <- as.numeric(eventmonth >= 0)
#'   	y <- 1 + (1 * eventmonth) + treat*2 + rnorm(length(eventmonth))
#'   	df <- data.frame(y=y, eventmonth=eventmonth)
#'   	its_loess(df, rvar="eventmonth", outcome="y") 		
#' 
#' 		}
#' 
#' @references 
#'	
#' 
#' 
#' @export
its_loess <- function(df,rvar,outcome,span=0.8,degree=1,othervar=NULL, only_estimates=TRUE){ 
	
	df <- as.data.frame(df)	
	dat <- df[,c(rvar, outcome, othervar)]
	dat$treat <- as.numeric(dat[,rvar] >= 0) 
	dat[,'span'] <- span
	dat[,'degree'] <- degree

	spec <- paste(outcome, " ~ ",rvar, sep="")
	param <- loess.control(surface='direct', statistics='exact')

	set.seed(42)
	m0 <- loess(spec, data=dat[dat[,'treat']==0,], control=param, degree=degree, span=span)
	m1 <- loess(spec, data=dat[dat[,'treat']==1,], control=param, degree=degree, span=span)
	
	yhat0 <- predict(m0, se=TRUE, newdata=dat[dat[,rvar]<=0,rvar] )
	yhat1 <- predict(m1, se=TRUE, newdata=dat[dat[,rvar]>=0,rvar] )

	dat[dat[,rvar]<=0, 'yhat0'] <- yhat0$fit
	dat[dat[,rvar]>=0, 'yhat1'] <- yhat1$fit

	dat[dat[,rvar]<=0, 'yhat0.se'] <- yhat0$se.fit
	dat[dat[,rvar]>=0, 'yhat1.se'] <- yhat1$se.fit

	dat[dat[,rvar]<=0, 'yhat0.lo'] <- yhat0$se.fit * qnorm(0.025) + yhat0$fit
	dat[dat[,rvar]>=0, 'yhat1.lo'] <- yhat1$se.fit * qnorm(0.025) + yhat1$fit

	dat[dat[,rvar]<=0, 'yhat0.hi'] <- yhat0$se.fit * qnorm(0.975) + yhat0$fit
	dat[dat[,rvar]>=0, 'yhat1.hi'] <- yhat1$se.fit * qnorm(0.975) + yhat1$fit

	res <- welch_test_2sided(mx=dat[, 'yhat1'], my=dat[, 'yhat0'], 
					  stderrx=dat[, 'yhat1.se'], stderry=dat[, 'yhat0.se'],
					  nx=yhat1$df, ny=yhat0$df)
	
	if (only_estimates==TRUE) {
		res <- na.omit(res)
		res$degree <- degree
		res$span <- span
		rownames(res) <- NULL
		return(res)
		}
	else return(cbind(dat, res))

	}