\name{findNeighborhood}
\alias{findNeighborhood}
\title{
Neighborhoods of an Allocation
}
\description{
Takes an allocation and the number of connections to find the neighborhoods, based on the definition of a neighborhood in the paper.
}
\usage{
findNeighborhood(object)
}
\arguments{
  \item{object}{ object of class "climbingTools". }
}
\details{
In the paper, a neighborhood is any point such that “all the individual allocations are identical except for two” and that “each input is “partnered” with one other input to which it may transfer dollars”. The number of inputs with which one input is partnered is the number of connections. 

Based on my interpretation, if we have 10 inputs and 3 connections, there will be 2 subsets of 4 inputs that will be transferring dollars from each other. Naturally, the remaining two inputs would be constant as if they are also partnered, there will be at least two other inputs that will end up having more than three connections, thus violating the condition. 

Therefore if we generate a neighborhood of allocations with N inputs and C connections, the number of neighborhoods = choose(c+1,2)*2*(number of subsets)
}
\value{
Matrix where rows represent neighborhoods.
}
\author{
Kanishka Malik
}
\examples{
object <- new("climbingTools", X1=c(1,2,3), Coeff_inputvars=c(1,2,3), Coeff_squares=c(1,2,3), Coeff_crossprod=c(1,2,3), Num_Connections=2)
\dontrun{
findNeighborhood(object)
}
}