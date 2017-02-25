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


make_permut_set <- function(df, rvar, bw, donut){
	minleft <- min(df[,rvar])
	placeboL <- minleft + bw 
	placeboR <- 0-bw-donut
	permutset <- df[,rvar][(df[,rvar] >= placeboL) & (df[,rvar] < placeboR)]
	return(permutset)
	}
