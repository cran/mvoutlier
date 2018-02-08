\name{color.plot}
\alias{color.plot}
\title{Color Plot}
\description{
The function color.plot plots the (two-dimensional) data using different symbols according to
the robust mahalanobis distance based on the mcd estimator with adjustment and using different colors 
according to the euclidean distances of the observations.
}
\usage{color.plot(x, quan=1/2, alpha=0.025, ...)}
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
\item{euclidean}{euclidean distances of the observations according to the minimum of the data.}
}
\details{
The function color.plot plots the (two-dimensional) data using different symbols 
(see function symbol.plot)
according to the robust mahalanobis distance based on the mcd estimator with adjustment
and using different 
colors according to the euclidean distances of the observations. Blue is typical for a 
little distance, whereas red is the opposite. In addition four ellipsoids 
are drawn, on which mahalanobis distances are constant. These constant values correspond to 
the 25\%, 50\%, 75\% and adjusted quantiles (see function arw) of the chi-square distribution
(see Filzmoser et al., 2005).
}
\author{Moritz Gschwandtner <\email{e0125439@student.tuwien.ac.at}>
\cr Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://cstat.tuwien.ac.at/filz/}}
\references{P. Filzmoser, R.G. Garrett, and C. Reimann.
Multivariate outlier detection in exploration geochemistry.
\emph{Computers & Geosciences}, 31:579-587, 2005.}
\seealso{\code{\link{symbol.plot}}, \code{\link{dd.plot}}, \code{\link{arw}}}
\examples{
# create data:
x <- cbind(rnorm(100), rnorm(100))
y <- cbind(rnorm(10, 5, 1), rnorm(10, 5, 1))
z <- rbind(x,y)
# execute:
color.plot(z, quan=0.75)
}
\keyword{dplot}
