\name{symbol.plot}
\alias{symbol.plot}
\title{Symbol Plot}
\description{
The function symbol.plot plots the (two-dimensional) data using different symbols according to
the robust mahalanobis distance based on the mcd estimator with adjustment.}
\usage{symbol.plot(x, quan=1/2, alpha=0.025, ...)}
\arguments{
\item{x}{two dimensional matrix or data.frame containing the data.}
\item{quan}{amount of observations which are used for mcd estimations. 
has to be between 0.5 and 1, default ist 0.5}
\item{alpha}{amount of observations used for calculating the adjusted quantile (see function arw).}
\item{...}{additional graphical parameters}
}
\value{
\item{outliers}{boolean vector of outliers}
\item{md}{robust mahalanobis distances of the data}
}
\details{The function symbol.plot plots the (two-dimensional) data using different symbols.
In addition a legend and four ellipsoids are drawn, on which mahalanobis distances are constant. 
As the legend shows, these constant values correspond to the 25\%, 50\%, 75\% and adjusted 
(see function arw) quantiles of the chi-square distribution.}
\author{Moritz Gschwandtner <\email{e0125439@student.tuwien.ac.at}>
\cr Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://cstat.tuwien.ac.at/filz/}}
\references{P. Filzmoser, R.G. Garrett, and C. Reimann.
Multivariate outlier detection in exploration geochemistry.
\emph{Computers & Geosciences}, 31:579-587, 2005.}
\seealso{\code{\link{dd.plot}}, \code{\link{color.plot}}, \code{\link{arw}}}
\examples{
# create data:
x <- cbind(rnorm(100), rnorm(100))
y <- cbind(rnorm(10, 5, 1), rnorm(10, 5, 1))
z <- rbind(x,y)
# execute:
symbol.plot(z, quan=0.75)
}
\keyword{dplot}
