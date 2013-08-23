\name{pcout}
\alias{pcout}
\title{ PCOut Method for Outlier Identification in High Dimensions }
\description{
  Fast algorithm for identifying multivariate outliers in high-dimensional
  and/or large datasets, using the algorithm of Filzmoser, Maronna, and
  Werner (CSDA, 2007).
}
\usage{
pcout(x, makeplot = FALSE, explvar = 0.99, crit.M1 = 1/3, crit.c1 = 2.5, 
   crit.M2 = 1/4, crit.c2 = 0.99, cs = 0.25, outbound = 0.25, ...)
}
\arguments{
  \item{x}{ a numeric matrix or data frame which provides the data for
  outlier detection }
  \item{makeplot}{ a logical value indicating whether a diagnostic plot
  should be generated (default to FALSE) }
  \item{explvar}{ a numeric value between 0 and 1 indicating how much variance
  should be covered by the robust PCs (default to 0.99) }
  \item{crit.M1}{ a numeric value between 0 and 1 indicating the quantile to
  be used as lower boundary for location outlier detection (default to 1/3) }
  \item{crit.c1}{ a positive numeric value used for determining the upper
  boundary for location outlier detection (default to 2.5) }
  \item{crit.M2}{ a numeric value between 0 and 1 indicating the quantile to
  be used as lower boundary for scatter outlier detection (default to 1/4) }
  \item{crit.c2}{ a numeric value between 0 and 1 indicating the quantile to
  be used as upper boundary for scatter outlier detection (default to 0.99) }
  \item{cs}{ a numeric value indicating the scaling constant for combined 
  location and scatter weights (default to 0.25) }
  \item{outbound}{ a numeric value between 0 and 1 indicating the outlier
  boundary for defining values as final outliers (default to 0.25) }
  \item{\dots}{ additional plot parameters, see help(par) }
}
\details{
  Based on the robustly sphered data, semi-robust principal components
  are computed which are needed for determining distances for each
  observation. Separate weights for location and scatter outliers
  are computed based on these distances. The combined weights are
  used for outlier identification.
}
\value{
  \item{wfinal01 }{0/1 vector with final weights for each observation; 
  weight 0 indicates potential multivariate outliers.}
  \item{wfinal }{numeric vector with final weights for each observation; 
  small values indicate potential multivariate outliers.}
  \item{wloc }{numeric vector with weights for each observation; 
  small values indicate potential location outliers.}
  \item{wscat }{numeric vector with weights for each observation; 
  small values indicate potential scatter outliers.}
  \item{x.dist1 }{numeric vector with  distances for location outlier detection.}
  \item{x.dist2 }{numeric vector with  distances for scatter outlier detection.}
  \item{M1 }{upper boundary for assigning weight 1 in location outlier detection.}
  \item{const1 }{lower boundary for assigning weight 0 in location outlier detection.}
  \item{M2 }{upper boundary for assigning weight 1 in scatter outlier detection.}
  \item{const2 }{lower boundary for assigning weight 0 in scatter outlier detection.}
}
\references{
  P. Filzmoser, R. Maronna, M. Werner.
  Outlier identification in high dimensions,
  \emph{Computational Statistics and Data Analysis}, 52, 1694-1711, 2008.
}
\author{Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\seealso{ \code{\link{sign1}}, \code{\link{sign2}} }
\examples{
# geochemical data from northern Europe
data(bsstop)
x=bsstop[,5:14]
# identify multivariate outliers
x.out=pcout(x,makeplot=FALSE)
# visualize multivariate outliers in the map
op <- par(mfrow=c(1,2))
data(bss.background)
pbb(asp=1)
points(bsstop$XCOO,bsstop$YCOO,pch=16,col=x.out$wfinal01+2)
title("Outlier detection based on pcout")
legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))

# compare with outlier detection based on MCD:
x.mcd <- robustbase::covMcd(x)
pbb(asp=1)
points(bsstop$XCOO,bsstop$YCOO,pch=16,col=x.mcd$mcd.wt+2)
title("Outlier detection based on MCD")
legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))
par(op)
}
\keyword{multivariate}
\keyword{robust}
