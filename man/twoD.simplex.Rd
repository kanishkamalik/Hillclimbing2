\name{twoD.simplex}
\alias{twoD.simplex}
\title{
2D Simplex
}
\description{
Creates a 2D simplex representing all allocations generated through 3 inputs with a budget inherited from an object of class "hillclimbing".
}
\usage{
twoD.simplex(object)
}
\arguments{
  \item{object}{
      object of class "hillclimbing" which will contain the Budget.
}
}
\value{
2D simplex where each axis represents an input scaled to 100.
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
twoD.simplex(hf)
}
}