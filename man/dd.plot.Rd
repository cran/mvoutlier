\name{dd.plot}
\alias{dd.plot}
\title{Distance-Distance Plot}
\description{
The function dd.plot plots the classical mahalanobis distance of the data against the robust 
mahalanobis distance based on the mcd estimator. Different symbols (see function symbol.plot)
and colours (see function color.plot) are used depending 
on the mahalanobis and euclidean distance of the observations (see Filzmoser et al., 2005).}
\usage{dd.plot(x, quan=1/2, alpha=0.025, ...)}
\arguments{
\item{x}{matrix or data frame containing the data}
\item{quan}{amount of observations which are used for mcd estimations. 
has to be between 0.5 and 1, default ist 0.5}
\item{alpha}{amount of observations used for calculating the adjusted quantile (see function arw).}
\item{...}{additional graphical parameters}
}
\value{
\item{outliers}{boolean vector of outliers}
\item{md.cla}{mahalanobis distances of the observations based on classical estimators of
location and scatter.}
\item{md.rob}{mahalanobis distances of the observations based on robust estimators of
location and scatter (mcd).}
}
\author{Moritz Gschwandtner <\email{e0125439@student.tuwien.ac.at}>
\cr Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\references{P. Filzmoser, R.G. Garrett, and C. Reimann (2005).
Multivariate outlier detection in exploration geochemistry.
\emph{Computers & Geosciences}. To appear.}
\seealso{\code{\link{symbol.plot}}, \code{\link{color.plot}}, \code{\link{arw}}, \code{\link[rrcov]{covPlot}}}
\examples{
# create data:
x <- cbind(rnorm(100), rnorm(100))
y <- cbind(rnorm(10, 3, 1), rnorm(10, 3, 1))
z <- rbind(x,y)
# execute:
dd.plot(z)
#
# Identify multivariate outliers for Co-Cu-Ni in humus layer of Kola data:
data(humus)
dd.plot(log(humus[,c("Co","Cu","Ni")]))
}
\keyword{dplot}
