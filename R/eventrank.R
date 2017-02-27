#' Creates a vector of ranks from the left and right of a cutpoint 
#' 
#' @param x a vector to be ranked
#' @param cutpoint a scalar number used as cut-point (zero by default)
#' 
#' @return vector of ranks. If an element in x is equal to the cut-point, it is assigned '0'. 
#' 
#' @examples 
#'  \dontrun{
#'
#'  x <- sample(0:10, 10, replace=TRUE)
#'  data.frame(x=x, r=eventrank(x, 5) ) 
#' 
#'  } 
#' 
#' 
#' 
#' @export
eventrank <- function(x, cutpoint=0){
	N <- length(x)
	left <- dense_rank(x[x < cutpoint]*(-1))*(-1)
	right <- dense_rank(x[x > cutpoint])
	ranks <- rep(NA, N)
	ranks[x < cutpoint] <- left
	ranks[x > cutpoint] <- right
	ranks[x == cutpoint] <- 0
	return(ranks)
	}
