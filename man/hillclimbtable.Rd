\name{hillclimbtable}
\alias{hillclimbtable}
\title{
Table Summarizing Ascent Methods
}
\description{
Table containing mean normalized profits for all ascent methods for all connections from one to five.
}
\usage{
hillclimbtable(object)
}
\arguments{
  \item{object}{
object of class "hillclimbing" from which this method will inherit the Budget, Number of Inputs, Landscapes, and bounds of the coefficient vectors.
}
}
\value{
Table containing mean normalized profits for all three ascent methods and all connections from 1 to 5.
}
\author{
Kanishka Malik
}
\examples{
##object of class "hillclimbing"
hf <- new("hillclimbing", approach="steepascent", 
  Budget=20, Num_Inputs=6, Num_Connections=5, 
  Num_Landscapes=100, Upperbound_inputvars = 1, 
  Lowerbound_inputvars = -1, Upperbound_squares = 1, 
  Lowerbound_squares = -1, Upperbound_crossprod = 1, 
  Lowerbound_crossprod = -1)
\dontrun{
hillclimbtable(hf)
}
}