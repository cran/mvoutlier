\name{uni.plot}
\alias{uni.plot}
\title{Univariate Presentation of Multivariate Outliers}
\description{
The function uni.plot plots each variable of x parallel in a one-dimensional 
scatter plot and in addition marks multivariate outliers.}
\usage{uni.plot(x, symb=FALSE, quan=1/2, alpha=0.025, ...)}
\arguments{
\item{x}{matrix or data.frame containing the data.}
\item{symb}{logical. if FALSE, only two colors and no special symbols are used. outliers are
marked red. 
if TRUE different symbols (cross means big value, circle means little value) according to 
the robust mahalanobis distance based on the mcd estimator and different colors (red means 
big value, blue means little value) according to the euclidean distances of the observations 
are used.}
\item{quan}{amount of observations which are used for mcd estimations. 
has to be between 0.5 and 1, default ist 0.5}
\item{alpha}{amount of observations used for calculating the adjusted quantile (see function arw).}
\item{...}{additional graphical parameters}
}
\value{
\item{outliers}{boolean vector of outliers}
\item{md}{robust multivariate mahalanobis distances of the data}
\item{euclidean}{(only if symb=TRUE) multivariate euclidean distances of the observations 
according to the minimum of the data.}
}
\details{The function uni.plot shows the mutlivariate outliers in the single variables
by one-dimensional scatter plots. 
If symb=FALSE (default), only two colors and no special symbols are used to
mark multivariate outliers (the outliers are marked red).
If symb=TRUE different symbols and colors are used.
The symbols (cross means big value, circle means little value) are selected according to
the robust mahalanobis distance based on the adjusted mcd estimator (see function symbol.plot)
Different colors (red means big value, blue means little value)
according to the euclidean distances of the observations (see function color.plot)
are used. For details see Filzmoser et al. (2005).}
\author{Moritz Gschwandtner <\email{e0125439@student.tuwien.ac.at}>
\cr Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\references{P. Filzmoser, R.G. Garrett, and C. Reimann (2005).
Multivariate outlier detection in exploration geochemistry.
\emph{Computers & Geosciences}. To appear.}
\seealso{\code{\link{map.plot}}, \code{\link{symbol.plot}}, \code{\link{color.plot}}, \code{\link{arw}}}
\examples{
data(swiss)
uni.plot(swiss)
#
# Geostatistical data:
data(humus) # Load humus data
uni.plot(log(humus[, c("As", "Cd", "Co", "Cu", "Mg", "Pb", "Zn")]),symb=TRUE)
}
\keyword{dplot}
