\name{landscapePlot}
\alias{landscapePlot}
\title{
Random Landscape over 3 Inputs
}
\description{
A randomly generated landscape created using plot3d() 
}
\usage{
landscapePlot(object)
}
\arguments{
  \item{object}{
object of class "hillclimbing" from which it will inherit the budget
}
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
landscapePlot(hf)
}
}