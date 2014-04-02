\name{allocationProfit}
\alias{allocationProfit}
\title{
allocationProfit()
}
\description{
Inherits the allocation and coefficient vectors from an object of class "climbingTools" to return a profit value using the profit function defined in the paper.
}
\usage{
allocationProfit(object)
}
\arguments{
  \item{object}{ Object of class "climbingTools" that contains an allocation of a given number of inputs and the coefficient vectors for the input variables, squares and the cross product }
}
\details{
To calculate the elements of the cross product, this package uses a special matrix containing 1s and -1s, generated using a function contrmat (from package "multcomp").
When every row of this matrix is multiplied with the allocation, it generates a vector containing the elements of the cross product.
This technique is used in the method allocationProfit() which accepts an allocation, and the coefficient vectors to calculate the profit using the given profit function.
}
\value{
Profit
}
\author{
Kanishka Malik
}
\examples{
##object of class "climbing tools"
object <- new("climbingTools", X1=c(1:3), Coeff_inputvars=c(1:3), Coeff_squares=c(1:3), Coeff_crossprod=c(1:3), Num_Connections=2)
\dontrun{
allocationProfit(object)
}
}