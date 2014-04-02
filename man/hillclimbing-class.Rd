\name{hillclimbing-class}
\Rdversion{1.1}
\docType{class}
\alias{hillclimbing-class}
\alias{hillclimb, hillclimbing-method}
\alias{hillclimbboxplot, hillclimbing-method}
\alias{hillclimbtable, hillclimbing-method}
\alias{landscapePlot, hillclimbing-method}
\alias{twoD.simplex, hillclimbing-method}

\title{Class \code{"hillclimbing"}}
\description{
Class "hillclimbing" will contain parameters which will be used by the different hillclimbing methods set as generic.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("hillclimbing", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{approach}:}{Object of class \code{"character"}: the ascent method; either "steepascent", "medianascent" & "leastascent". }
    \item{\code{Budget}:}{Object of class \code{"numeric"} Budget of the firms. }
    \item{\code{Num_Inputs}:}{Object of class \code{"numeric"} Number of Inputs where this budget will be allocated. }
    \item{\code{Num_Connections}:}{Object of class \code{"numeric"}: Number of connections between inputs. }
    \item{\code{Num_Landscapes}:}{Object of class \code{"numeric"}: Number of landscapes across which we would like to run simulations. }
    \item{\code{Upperbound_inputvars}:}{.Object of class \code{"numeric"}: Upper bound of coefficients of input variables}
    \item{\code{Lowerbound_inputvars}:}{Object of class \code{"numeric"}: Lower bound of coefficients of input variables}
    \item{\code{Upperbound_squares}:}{Object of class \code{"numeric"}: upper bound of coefficients of squares}
    \item{\code{Lowerbound_squares}:}{Object of class \code{"numeric"}: lower bound of coefficients of squares. }
    \item{\code{Upperbound_crossprod}:}{Object of class \code{"numeric"}: upper bound of coefficients of cross product. }
    \item{\code{Lowerbound_crossprod}:}{Object of class \code{"numeric"}: lower bound of coefficients of cross product. }
  }
}
\details{
For an object of class "hillclimbing" to be valid, the number of connections must be less than the number of inputs. And the approach name must also be either "medianascent", "leastascent" or "steepascent". 
}

\section{Methods}{
  \describe{
    \item{hillclimb}{\code{signature(object = "hillclimbing")}: returns standard error and mean of normalized profits and steps taken to reach the maxima. }
    \item{hillclimbboxplot}{\code{signature(object = "hillclimbing")}: multiple boxplots of mean normalized profits for all ascent methods across all connections from one to five. }
    \item{hillclimbtable}{\code{signature(object = "hillclimbing")}: table containing mean normalized profits for all ascent methods for all connections from one to five. }
    \item{landscapePlot}{\code{signature(object = "hillclimbing")}: method which inherits the budget from the object to create a randomly generated landscape for three inputs. }
    \item{twoD.simplex}{\code{signature(object = "hillclimbing")}: method which inherits the budget from the object and returns a 2D simplex scaled to hundred }
	 }
}
\author{
Kanishka Malik
}
\seealso{
\code{\link{hillclimb}}, \code{\link{hillclimbboxplot}}, 
\code{\link{hillclimbtable}}, \code{\link{landscapePlot}},
\code{\link{hillclimb}}, \code{\link{twoD.simplex}}
}
\examples{
showClass("hillclimbing")
}