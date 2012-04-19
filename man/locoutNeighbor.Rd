\name{locoutNeighbor}
\alias{locoutNeighbor}
\title{Diagnostic plot for identifying local outliers with varying size of neighborhood}
\description{
Computes global and pairwise Mahalanobis distances for visualizing global and
local multivariate outliers. The size of the neighborhood (number of neighbors) is
varying, but the fraction of neighbors is fixed. 
}
\usage{
locoutNeighbor(dat, X, Y, propneighb = 0.1, variant = c("dist", "knn"), usemax = 1/3, 
   npoints = 50, chisqqu = 0.975, indices = NULL, xlab = NULL, ylab = NULL, colall = gray(0.7), 
   colsel = 1, ...)
}
\arguments{
  \item{dat}{multivariate data set (without coordinates)}
  \item{X}{X coordinates of the data points}
  \item{Y}{Y coordinates of the data points}
  \item{propneighb}{proportion of neighbors to be included in tolerance ellipse}
  \item{variant}{either search for neighbors according to the Eucl.Distance, or according to kNN}
  \item{usemax}{for either variant: give fraction of points (max Dist) that is used for the plot}
  \item{npoints}{computation is done at most at npoints points}
  \item{chisqqu}{quantile of the chisquare distribution for splitting the plot}
  \item{indices}{if this is not NULL, these should be indices of observations to be highlighted}
  \item{xlab}{x-axis label for plot}
  \item{ylab}{y-axis label for plot}
  \item{colall}{color for lines if indices is NULL}
  \item{colsel}{color for lines if indices are selected}
  \item{\dots}{additional parameters for plotting}
}
\details{
For this diagnostic tool, the number of neighbors is varied up to a fraction of usemax observations.
Then propneighb (called beta) is fixed, and for each observation we compute the degree of
isolation from a fraction of 1-beta of its neighbors. Neighborhood can be defined either
via the Euclidean distance or by k-Nearest-Neighbors. For computational reasons, all computations
are done at most for npoints points. The critical value for outliers is the quantile chisqqu
of the chisquare distribution.
One can also provide indices of observations (for indices). Then the corresponding lines in
the plots will be highlighted.
}
\value{
  \item{indices.reg}{indices of the (selected) observations being regular observations}
  \item{indices.out}{indices of the (selected) observations being golbal outliers}
}
\references{
P. Filzmoser, A. Ruiz-Gazen, and C. Thomas-Agnan: Identification of local multivariate
outliers. Submitted for publication, 2012.
}
\author{Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\seealso{ \code{\link{locoutPercent}}, \code{\link{locoutSort}} }
\examples{
# use data from illustrative example in paper:
data(X)
data(Y)
data(dat)
res <- locoutNeighbor(dat,X,Y,variant="knn",usemax=1,chisqqu=0.975,indices=c(1,11,24,36),
              propneighb=0.1,npoints=100)
}
\keyword{multivariate}
\keyword{robust}
