% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eventrank.R
\name{eventrank}
\alias{eventrank}
\title{Creates a vector of ranks from the left and right of a cutpoint}
\usage{
eventrank(x, cutpoint = 0)
}
\arguments{
\item{x}{a vector to be ranked}

\item{cutpoint}{a scalar number used as cut-point (zero by default)}
}
\value{
vector of ranks. If an element in x is equal to the cut-point, it is assigned '0'.
}
\description{
Creates a vector of ranks from the left and right of a cutpoint
}
\examples{
 \dontrun{

 x <- sample(0:10, 10, replace=TRUE)
 data.frame(x=x, r=eventrank(x, 5) ) 

 } 



}

