\name{hillclimbing-package}
\alias{hillclimbing-package}
\alias{hillclimbing}
\docType{package}
\title{
Hill Climbing Models
}
\description{
This package implements the ideas presented in the paper "Local Hillclimbing on an Economic Landscape". The package allows us to determine the optimal strategy a firm must apply towards allocating its budget, in order to maximize profits. We can individually test the different ascent methods a firm can use, and set the parameters of the economic landscapes across which it will climb. It also contains other useful tools that not only complement the methods of this package, but which will enhance our understanding of this paper, such as methods to find neighborhoods or graphical representations of simulations.
}
\details{
\tabular{ll}{
Package: \tab hillclimbing\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-01\cr
License: \tab GPL-2\cr
Depends: \tab methods, ggplot2, mvtnorm, splines, survival, TH.data, ggtern, EMT, multcomp, rgl\cr
}
}
\author{
Kanishka Malik <kanishka.college@gmail.com>

Maintainer: Kanishka Malik <kanishka.college@gmail.com>
}
\references{
Genolini, Christophe, A (Not So) Short Introduction to S4: Object Oriented Programming in R. August 20, 2008

Leisch, Friedrich, Creating R Packages: A Tutorial. Department of Statistics, Ludwig-Maximilians-Universitat Munchen, and R Development Core Team. September 14, 2009

Kane, David, Local Hillclimbing on an Economic Landscape, John F. Kennedy School of Government Harvard University, Cambridge, 1996.

Writing R extensions, R Development Core Team, 2013

Pages, Herve, Implementing S4 objects in your package, Fred Hutchinson Cancer Research Center, Seattle, 2011.

}

\seealso{
\code{\link{hillclimb}}, \code{\link{hillclimbboxplot}}, 
\code{\link{hillclimbtable}}, \code{\link{landscapePlot}},
\code{\link{findNeighborhood}}, \code{\link{allocationProfit}},
\code{\link{hillclimb}}, \code{\link{twoD.simplex}}
}

\examples{
##object of class "hillclimbing"
hf <- new("hillclimbing", approach="steepascent", 
  Budget=20, Num_Inputs=6, Num_Connections=5, 
  Num_Landscapes=100, Upperbound_inputvars = 1, 
  Lowerbound_inputvars = -1, Upperbound_squares = 1, 
  Lowerbound_squares = -1, Upperbound_crossprod = 1, 
  Lowerbound_crossprod = -1)
##method hillclimb() which returns the mean normalized profits and steps for a given ascent method
\dontrun{
hillclimb(object)
}
##method hillclimbboxplot() which returns mutliple boxplots containing the mean normalized profits for all methods, over all connections from one to five
\dontrun{
hillclimbboxplot(object)
}
}
