#' Estimates via Local Polynomial Regression 
#'
#' \code{its_loess} estimates the intercept shift of a time series at a cut-point.
#'
#' @param df (required) \code{data.frame} containing all variables 
#' @param rvar (required) the name of the running variable in \code{df}
#' @param outcome (required) the name of the outcome variable in \code{df} 
#' @param span either a scalar or a vector of length 2 for the span to be used to the left and right of the cut-point.
#' @param degree either a scalar or a vector of length 2 defining the degree of the polynomial to the left and right of the cut-point.
#' @param donut either a scalar or a vector of length 2 defining the length of the period to the left (right) of the cut-point for which the data are dropped (on the scale of \code{rvar}).
#' @param othervar vector of variable names that will be included into the returned \code{data.frame} (if \code{only_estimate=FALSE}).
#' @param only_estimate if \code{TRUE} only returns the estimates at the cut-point (default)
#' @param na.action parameter passed on to \code{loess} about how to handle missing values
#' 
#' @details 
#' This function estimates two local polynomial regressions via \code{\link{loess}} to the left and 
#' right of the cut-point at zero. The predicted values from the regressions are used to form an 
#' estimate of the intercept shift at the cut-point. The reported t-statistic and p-value is based 
#' on a standard Welch's t-test with unequal variances. 
#' 
#' As span increases, the regression becomes a linear regression with a linear (degree=1) or quadratic trend term (degree=2).
#' The loess parameters are set such that the surface and statistics are computed exactly, i.e. \code{loess.control(surface='direct', statistics='exact')}.
#' 
#' @return If \code{only_estimate=TRUE} a \code{data.frame} with a single row and entries for the point estimate of the intercept shift (\code{diff.est}), 
#' standard error (\code{diff.se}), t-statistic (\code{diff.tstat}), p-value (\code{diff.pval}) as well as the degree and span used in the
#' estimation. If \code{only_estimate=FALSE}, a \code{data.frame} that contains the variables as named in \code{rvar,outcome,othervar} as 
#' well as the predicted values (\code{yhat*}), standard errors (\code{yhat*.se}) and upper/lower bound of the 95\% confidence interval (based on a normal distribution)
#' (\code{yhat*.lo},\code{yhat*hi}) from the \code{\link{loess}} regressions to the left (*=0) and to the 
#' right (*=1) of the cut-point.
#' 	
#' @seealso \code{\link{loess}}.
#' 
#' @examples 
#'  \dontrun{
#'		
#'   	eventmonth <- -12:12
#'   	treat <- as.numeric(eventmonth >= 0)
#'   	y <- 1 + (1 * eventmonth) + treat*10 + rnorm(length(eventmonth))
#'   	df <- data.frame(y=y, eventmonth=eventmonth)
#'   	its_loess(df, rvar="eventmonth", outcome="y") 		
#' 
#' 		m <- its_loess(df, rvar="eventmonth", outcome="y", only_estimate=FALSE)
#' 		with(m, plot(eventmonth, y)) 		
#' 	 	with(m, lines(eventmonth, yhat0, col='blue'))
#' 		with(m, lines(eventmonth, yhat1, col='blue'))
#' 
#' 		}
#' 
#' @export
its_loess <- function(df,rvar,outcome, span=0.8, degree=1, donut=0, 
	othervar=NULL, only_estimates=TRUE, na.action="na.omit"){ 

	list[spanL,spanR] <- parseLR(span)
	list[degreeL,degreeR] <- parseLR(degree)
	list[donutL,donutR] <- parseLR(donut)

	df <- as.data.frame(df)	
	df <- make_inSampleVar(df,rvar=rvar, bwL=Inf, bwR=Inf, donutL=donutL, donutR=donutR)

	dat <- df[,c(rvar, outcome, othervar,'inSample')]
	dat$treat <- as.numeric(dat[,rvar] >= 0) 

	dat[,'span'] <- paste(span, collapse="|")
	dat[,'degree'] <- paste(degree, collapse="|")
	dat[,'donut'] <- paste(donut, collapse="|")

	spec <- paste(outcome, " ~ ",rvar, sep="")
	param <- loess.control(surface='direct', statistics='exact')

	set.seed(42)
	m0 <- loess(spec, data=dat[dat[,'treat']==0 & dat[,'inSample']==1,], control=param, 
		degree=degreeL, span=spanL, na.action=na.action)
	m1 <- loess(spec, data=dat[dat[,'treat']==1 & dat[,'inSample']==1,], control=param, 
		degree=degreeR, span=spanR, na.action=na.action)
	
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
		res$degree <- paste(degree, collapse="|")
		res$span <- paste(span, collapse="|")
		res$donut <- paste(donut, collapse="|")
		rownames(res) <- NULL
		return(res)
		}
	else return(cbind(dat, res))

	}

