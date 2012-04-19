\name{locoutSort}
\alias{locoutSort}
\title{Interactive diagnostic plot for identifying local outliers}
\description{
Computes global and pairwise Mahalanobis distances for visualizing global and
local multivariate outliers. The plot is split into regular (left) and global 
(right) outliers, and points can be selected interactively. In a second plot,
these points are shown by spatial coordinates.
}
\usage{
locoutSort(dat, X, Y, distc = NULL, k = 10, propneighb = 0.1, chisqqu = 0.975, sel = NULL, ...)
}
\arguments{
  \item{dat}{multivariate data set (without coordinates)}
  \item{X}{X coordinates of the data points}
  \item{Y}{Y coordinates of the data points}
  \item{distc}{maximum distance to search for neighbors; if nothing is provided, k for kNN is used}
  \item{k}{number of nearest neighbors to search - not taken if a value for dist is provided}
  \item{propneighb}{proportion of neighbors to be included in tolerance ellipse}
  \item{chisqqu}{quantile of the chisquare distribution for splitting the plot}
  \item{sel}{optional list with x and y, i.e. coordinates with selected polygon}
  \item{\dots}{additional parameters for plotting}
}
\details{
For this diagnostic tool, the number of neighbors is fixed, and propneighb (called beta) is also fixed.
For each observation we compute the degree of isolation from a fraction of 1-beta of its neighbors. 
The observations are sorted according to this degree of isolation, and this sorted index forms
the x-axis of the left plot. This plot is also split into regular (left) and global (right) outliers.
Then one can select with the mouse a region in this plot, meaning an observation and (some of)
its neighbors. Alternatively, this region can be supplied by sel. The selected observations
are then shown in the right plot. Links to the neighbors are also shown.
}
\value{
list(sel=sel,index.regular=res$indices.regular,index.outliers=res$indices.outliers)
  \item{sel}{plot coordinates of the selected region}
  \item{indices.reg}{indices of the bservations being regular observations}
  \item{indices.out}{indices of the observations being golbal outliers}
}
\references{
P. Filzmoser, A. Ruiz-Gazen, and C. Thomas-Agnan: Identification of local multivariate
outliers. Submitted for publication, 2012.
}
\author{Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\seealso{ \code{\link{locoutPercent}}, \code{\link{locoutNeighbor}} }
\examples{
# use data from illustrative example in paper:
data(X)
data(Y)
data(dat)
sel <- locoutSort(dat,X,Y,k=10,propneighb=0.1,chisqqu=0.975,
    sel=list(x=c(87.5,87.5,89.3,89.3),y=c(4.3,0.7,0.7,4.3)))
}
\keyword{multivariate}
\keyword{robust}
