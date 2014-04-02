\name{climbingTools-class}
\Rdversion{1.1}
\docType{class}
\alias{climbingTools-class}
\alias{allocationProfit,climbingTools-method}
\alias{findNeighborhood,climbingTools-method}
\alias{standardError,climbingTools-method}

\title{Class \code{"climbingTools"}}
\description{
Object of class "climbingTools". This will contain parameters that will be used by methods that are instrumental to the different hill climbing models, such as the methods to calculate profit or to find the neighborhoods.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("climbingTools", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{Coeff_inputvars}:}{Object of class \code{"numeric"}. Coefficient vector of the input variables }
    \item{\code{Coeff_squares}:}{Object of class \code{"numeric"}. Coefficient vector of the squares  }
    \item{\code{Coeff_crossprod}:}{Object of class \code{"numeric"}. Coefficient vector of the cross product }
    \item{\code{X1}:}{Object of class \code{"numeric"}. A vector representing an allocation or just profit values. }
    \item{\code{Num_Connections}:}{Object of class \code{"numeric"}. Number of connections. }
  }
}
\section{Methods}{
  \describe{
    \item{allocationProfit}{\code{signature(object = "climbingTools")}: returns the profit of an allocation }
    \item{findNeighborhood}{\code{signature(object = "climbingTools")}: returns the neighborhoods of an allocation }
    \item{standardError}{\code{signature(object = "climbingTools")}: returns the standard error of a dataset }
	 }
}
\author{
Kanishka Malik
}
\seealso{
\code{\link{findNeighborhood}}, \code{\link{allocationProfit}},
\code{\link{standardError}}
}
\examples{
showClass("climbingTools")
}