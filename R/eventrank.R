#' Creates a vector of ranks from the left and right of a cutoff 
#' 
#' @param x a vector to be ranked
#' @param cutoff a scalar number used as cutoff (zero by default)
#' 
#' @return vector of ranks. If an element in x is equal to the cutoff, it is assigned '0'. 
#' 
#' @examples 
#'  \dontrun{
#' 
#'  
#' 
#'  x <- sample(0:10, 10, replace=TRUE)
#'  data.frame(x=x, r=eventrank(x, 5) ) %>% arrange(x)
#' 
#'  } 
#' 
#' 
#' 
#' @export
eventrank <- function(x, cutoff=0){
	N <- length(x)
	left <- dense_rank(x[x < cutoff]*(-1))*(-1)
	right <- dense_rank(x[x > cutoff])
	ranks <- rep(NA, N)
	ranks[x < cutoff] <- left
	ranks[x > cutoff] <- right
	ranks[x == cutoff] <- 0
	return(ranks)
	}
