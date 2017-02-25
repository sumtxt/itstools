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

make_donut_df <- function(df, rvar, donut){
	df <- df[order(df[,rvar]),]
	inSample <- df[,rvar] < 0 | df[,rvar] >= donut
	df <- df[inSample==TRUE,]
	df[,rvar] <- df[,rvar]
	return(df)
	}

make_permut_set <- function(df, rvar, bwL, bwR, donut){
	minleft <- min(df[,rvar])
	placeboL <- minleft + bwL 
	placeboR <- 0-bwR-donut
	permutset <- df[,rvar][(df[,rvar] >= placeboL) & (df[,rvar] < placeboR)]
	return(permutset)
	}

# Identical to the Welch's t-test implement via t.test(..., var.equal=FALSE)
welch_test_2sided <- function(mx,my,stderrx,stderry,nx,ny){
    mu <- 0
    stderr <- sqrt(stderrx^2 + stderry^2)
    df <- stderr^4/(stderrx^4/(nx-1) + stderry^4/(ny-1))	
    tstat <- (mx - my - mu)/stderr
	pval <- 2 * pt(-abs(tstat), df)
    res <- data.frame(diff=mx - my, diff.se=stderr, diff.pval=pval, diff.tstat=tstat)
    return(res)
    }
