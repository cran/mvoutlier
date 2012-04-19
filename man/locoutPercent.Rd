\name{locoutPercent}
\alias{locoutPercent}
\title{Diagnostic plot for identifying local outliers with fixed size of neighborhood}
\description{
Computes global and pairwise Mahalanobis distances for visualizing global and
local multivariate outliers. The size of the neighborhood (number of neighbors) is
fixed, but the fraction of neighbors is varying.
}
\usage{
locoutPercent(dat, X, Y, dist = NULL, k = 10, chisqqu = 0.975, sortup = 10, sortlow = 90, 
  nlinesup = 10, nlineslow = 10, indices = NULL, xlab = "(Sorted) Index", 
  ylab = "Distance to neighbor", col = gray(0.7), ...)
}
\arguments{
  \item{dat}{multivariate data set (without coordinates)}
  \item{X}{X coordinates of the data points}
  \item{Y}{Y coordinates of the data points}
  \item{dist}{maximum distance to search for neighbors; if nothing is provided, k for kNN is used}
  \item{k}{number of nearest neighbors to search - not taken if a value for dist is provided}
  \item{chisqqu}{quantile of the chisquare distribution for splitting the plot}
  \item{sortup}{sort local outliers accorting to given percentage}
  \item{sortlow}{sort local inliers accorting to given percentage}
  \item{nlinesup}{number of lines to be plotted for upper part}
  \item{nlineslow}{number of lines to be plotted for lower part}
  \item{indices}{if this is not NULL, these should be indices of observations to be highlighted}
  \item{xlab}{x-axis label for plot}
  \item{ylab}{y-axis label for plot}
  \item{col}{color for lines}
  \item{\dots}{additional parameters for plotting}
}
\details{
For this diagnostic tool, the number of neighbors is fixed, but propneighb (called beta) is 
varied. For each observation we compute the degree of isolation from a fraction of 1-beta of 
its neighbors. Neighborhood can be defined either via the Euclidean distance or by k-Nearest-Neighbors. 
The critical value for outliers is the quantile chisqqu of the chisquare distribution.
One can also provide indices of observations (for indices). Then the corresponding lines in
the plots will be highlighted.
}
\value{
  \item{ret}{list containing indices of regular and outlying observations}
}
\references{
P. Filzmoser, A. Ruiz-Gazen, and C. Thomas-Agnan: Identification of local multivariate
outliers. Submitted for publication, 2012.
}
\author{Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\seealso{ \code{\link{locoutNeighbor}}, \code{\link{locoutSort}} }
\examples{
# use data from illustrative example in paper:
data(X)
data(Y)
data(dat)
res <- locoutPercent(dat,X,Y,k=10,chisqqu=0.975, indices=c(1,11,24,36))
}
\keyword{multivariate}
\keyword{robust}
