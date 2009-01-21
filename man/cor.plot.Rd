\name{cor.plot}
\alias{cor.plot}
\title{Correlation Plot: robust versus classical bivariate correlation}
\description{
The function cor.plot plots the (two-dimensional) data and adds two correlation ellipsoids, 
based on classical and robust estimation of location and scatter. Robust estimation 
can be thought of as estimating the mean and covariance of the 'good' part of the data.}
\usage{cor.plot(x, y, quan=1/2, alpha=0.025, ...)}
\arguments{
\item{x}{vector to be plotted against y and of which the correlation with y is calculated.}
\item{y}{vector to be plotted against x and of which the correlation with x is calculated.}
\item{quan}{amount of observations which are used for mcd estimations. 
has to be between 0.5 and 1, default ist 0.5}
\item{alpha}{Determines the size of the ellipsoids. An observation will be outside of the 
ellipsoid if its mahalanobis distance exceeds the 1-alpha quantile of the chi-squared 
distribution.}
\item{...}{additional graphical parameters}
}
\value{
\item{cor.cla}{correlation between x and y based on classical estimation of location and 
scatter}
\item{cor.rob}{correlation between x and y based on robust estimation of location and 
scatter}
}
\author{Moritz Gschwandtner <\email{e0125439@student.tuwien.ac.at}>
\cr Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\seealso{\code{\link[robustbase]{covMcd}}}
\examples{
# create data:
x <- cbind(rnorm(100), rnorm(100))
y <- cbind(rnorm(10, 3, 1), rnorm(10, 3, 1))
z <- rbind(x,y)
# execute:
cor.plot(z[,1], z[,2])
}
\keyword{dplot}
