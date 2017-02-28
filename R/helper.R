# Source: https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
   args <- as.list(match.call())
   args <- args[-c(1:2,length(args))]
   length(value) <- length(args)
   for(i in seq(along=args)) {
     a <- args[[i]]
     if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
   }
   x
}

parseLR <- function(x){
	varname <- deparse(substitute(span))
	if ( !(length(x) %in% c(1,2)) ) {
		stop(paste("Parameter ",varname," must be either a scalar or a vector of length 2.",sep=""))
	} 
	if (length(x)==1){
		xL <- xR <- x
		}
	if (length(x)==2){
		xL <- x[1]
		xR <- x[2]
		}
	return(list(xL=xL, xR=xR))
	}

lmHC <- function(..., HC="HC3", var='treat', estimates=TRUE){
	m <- lm(...)
	m <- coeftest(m, vcov = vcovHC(m, HC))
	if( estimates==TRUE) return(get_lmHC_est(m,var=var))
	else return(m)
	}

get_lmHC_est <- function(m,var){
	est <- m[rownames(m)==var,'Estimate'] 
	se <- m[rownames(m)==var,'Std. Error'] 
	pval <- m[rownames(m)==var,4] 
	low <- est + (se * qnorm(0.975))
	hig <- est + (se * qnorm(0.025))
	return(data.frame(est=est, lo=low, hi=hig, se=se, pval=pval))
	}


make_inSampleVar <- function(df, rvar, bwL, bwR, donutL, donutR){
	outer <-  (-(bwL+donutL) <= df[,rvar] ) & ( df[,rvar] <= (donutR+bwR) )
	inner <-  ( df[,rvar] < -donutL ) | ( donutR <= df[,rvar] )
	df[,'inSample'] <- as.numeric( outer & inner )
	return(df)	
	}

make_permut_set <- function(df, rvar, bwL, bwR, donutL, donutR){
	minleft <- min(df[,rvar])
	placeboL <- minleft + (bwL+donutL)
	placeboR <- 0-bwR-donutR
	permutset <- df[,rvar][(df[,rvar] >= placeboL) & (df[,rvar] < placeboR)]
	return(permutset)
	}

# Modified version of: t.test(..., var.equal=FALSE) 
welch_test_2sided <- function(mx,my,stderrx,stderry,nx,ny){
    mu <- 0
    stderr <- sqrt(stderrx^2 + stderry^2)
    df <- stderr^4/(stderrx^4/(nx-1) + stderry^4/(ny-1))	
    tstat <- (mx - my - mu)/stderr
	pval <- 2 * pt(-abs(tstat), df)
    res <- data.frame(diff=mx - my, diff.se=stderr, diff.pval=pval, diff.tstat=tstat)
    return(res)
    }