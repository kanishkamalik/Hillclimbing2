\name{hillclimbboxplot}
\alias{hillclimbboxplot}
\title{
hillclimbboxplot()
}
\description{
Method which returns multiple boxplots using functions from package "ggplot2" to represent the mean normalized profits for all ascent approaches, for all connections from 1 to 5
}
\usage{
hillclimbboxplot(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ object of class "hillclimbing" from which this method inherits the Budget, Number of landscapes, inputs and the bounds of the coefficients}
}
\value{
Multiple Boxplots
}
\author{
Kanishka Malik
}
\examples{
##object of class "hillclimbing"
hc <- new("hillclimbing", approach="steepascent",
    Budget=20, Num_Inputs=7, 
    Num_Connections=6, Num_Landscapes=100,
    Upperbound_inputvars = 1, Lowerbound_inputvars = -1, 
    Upperbound_squares = 1, Lowerbound_squares = -1, 
    Upperbound_crossprod = 1, Lowerbound_crossprod = -1)
\dontrun{
hillclimbboxplot(hc)
}
}