\name{map.plot}
\alias{map.plot}
\title{Plot Multivariate Outliers in a Map}
\description{
The function map.plot creates a map using geographical (x,y)-coordinates. 
This is thought for spatially dependent data of which coordinates are available. 
Multivariate outliers are marked.}
\usage{map.plot(coord, data, quan=1/2, alpha=0.025, symb=FALSE, plotmap=TRUE, 
    map="kola.background",which.map=c(1,2,3,4),map.col=c(5,1,3,4),
    map.lwd=c(2,1,2,1), ... )}
\arguments{
\item{coord}{(x,y)-coordinates of the data}
\item{data}{matrix or data.frame containing the data.}
\item{quan}{amount of observations which are used for mcd estimations. 
has to be between 0.5 and 1, default ist 0.5}
\item{alpha}{amount of observations used for calculating the adjusted quantile (see function arw).}
\item{symb}{logical for plotting special symbols (see details).}
\item{plotmap}{logical for plotting the background map.}
\item{map}{see plot.kola.background()}
\item{which.map}{see plot.kola.background()}
\item{map.col}{see plot.kola.background()}
\item{map.lwd}{see plot.kola.background()}
\item{...}{additional graphical parameters}
}
\value{
\item{outliers}{boolean vector of outliers}
\item{md}{robust mahalanobis distances of the data}
\item{euclidean}{(only if symb=TRUE) euclidean distances of the observations according to 
the minimum of the data.}
}
\details{The function map.plot shows mutlivariate outliers in a map.
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
\url{http://cstat.tuwien.ac.at/filz/}}
\references{P. Filzmoser, R.G. Garrett, and C. Reimann. 
Multivariate outlier detection in exploration geochemistry.
\emph{Computers & Geosciences}, 31:579-587, 2005.}
\seealso{\code{\link{symbol.plot}}, \code{\link{color.plot}}, \code{\link{arw}}}
\examples{
data(humus) # Load humus data
xy <- humus[,c("XCOO","YCOO")] # X and Y Coordinates
myhumus <- log(humus[, c("As", "Cd", "Co", "Cu", "Mg", "Pb", "Zn")])
map.plot(xy, myhumus, symb=TRUE)
}
\keyword{dplot}
