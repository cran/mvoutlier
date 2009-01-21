\name{sign2}
\alias{sign2}
\title{ Sign Method for Outlier Identification in High Dimensions - Sophisticated Version}
\description{
  Fast algorithm for identifying multivariate outliers in high-dimensional
  and/or large datasets, using spatial signs, see Filzmoser, Maronna, and
  Werner (CSDA, 2007). The computation of the distances is based on
  principal components.
}
\usage{
sign2(x, makeplot = FALSE, explvar = 0.99, qcrit = 0.975, ...)
}
\arguments{
  \item{x}{ a numeric matrix or data frame which provides the data for
  outlier detection }
  \item{makeplot}{ a logical value indicating whether a diagnostic plot
  should be generated (default to FALSE) }
  \item{explvar}{ a numeric value between 0 and 1 indicating how much variance
  should be covered by the robust PCs (default to 0.99) }
  \item{qcrit}{ a numeric value between 0 and 1 indicating the quantile to
  be used as critical value for outlier detection (default to 0.975) }
  \item{\dots}{ additional plot parameters, see help(par) }
}
\details{
  Based on the robustly sphered and normed data, robust principal components
  are computed which are needed for determining distances for each
  observation. The distances are transformed to approach chi-square
  distribution, and a critical value is then used as outlier cutoff.
}
\value{
  \item{wfinal01 }{0/1 vector with final weights for each observation;
  weight 0 indicates potential multivariate outliers.}
  \item{x.dist }{numeric vector with  distances used for outlier detection.}
  \item{const }{ outlier  cutoff value.}
}
\references{
  P. Filzmoser, R. Maronna, M. Werner.
  Outlier identification in high dimensions,
  \emph{Computational Statistics and Data Analysis}, 52, 1694--1711, 2008.

  N. Locantore, J. Marron, D. Simpson, N. Tripoli, J. Zhang, and K. Cohen.
  Robust principal components for functional data,
  \emph{Test} 8, 1-73, 1999.
}
\author{Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\seealso{ \code{\link{pcout}}, \code{\link{sign1}} }
\examples{
# geochemical data from northern Europe
data(bsstop)
x=bsstop[,5:14]
# identify multivariate outliers
x.out=sign2(x,makeplot=FALSE)
# visualize multivariate outliers in the map
op <- par(mfrow=c(1,2))
data(bss.background)
pbb(asp=1)
points(bsstop$XCOO,bsstop$YCOO,pch=16,col=x.out$wfinal01+2)
title("Outlier detection based on signout")
legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))

# compare with outlier detection based on MCD:
require(robustbase)
x.mcd=covMcd(x)
pbb(asp=1)
points(bsstop$XCOO,bsstop$YCOO,pch=16,col=x.mcd$mcd.wt+2)
title("Outlier detection based on MCD")
legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))
par(op)
}
\keyword{multivariate}
\keyword{robust}
