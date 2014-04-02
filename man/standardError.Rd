\name{standardError}
\alias{standardError}
\title{
Standard Error of a Vector
}
\description{
Calculates the standard error of a vector inherited from object of class "climbingTools"
}
\usage{
standardError(object)
}
\arguments{
  \item{object}{
      object of class "climbingTools" from which the vector X1 will be inherited
}
}
\value{
standard error of the vector if the vector is non-empty}
\author{
Kanishka Malik
}
\examples{
##object of class "climbing tools"
object <- new("climbingTools", X1=c(1:3), Coeff_inputvars=c(1:3), Coeff_squares=c(1:3), Coeff_crossprod=c(1:3), Num_Connections=2)
\dontrun{
standardError(object)
}
}