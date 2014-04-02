\name{hillclimb}
\alias{hillclimb}
\title{
hillclimb()
}
\description{
Method to find mean normalized profits and mean steps to reach the local maxima using a given ascent method, over a given number of landscapes
}
\usage{
hillclimb(object)
}
\arguments{
  \item{object}{
object of class "hillclimbing" whose slots will contain the Budget, Number of Connections, Landscapes, Inputs, ascent method and bounds of the coefficients.
}
}
\details{
The paper evaluates three different approaches for ascending a landscape. This package contains a method hillclimb() which inherits an object of class "hillclimbing" which contains the ascent method, budget, inputs, connections, landscapes and bounds defined in an object of class “hillclimbing”. It returns the mean normalized profits and the mean number of steps taken to reach the local maxima over the required number of landscapes.
}
\value{
Mean & Standard Error of the normalized Profit and steps to reach the landscape. 
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
hillclimb(hc)
}
}