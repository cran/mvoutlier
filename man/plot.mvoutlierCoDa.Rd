\name{plot.mvoutlierCoDa}
\alias{plot.mvoutlierCoDa}
\title{Plots for interpreting multivatiate outliers of CoDa}
\description{
Plots the computed information by \code{mvoutlier.CoDa} for supporting the interpretation
of multivariate outliers in case of compositional data.}
\usage{
\method{plot}{mvoutlierCoDa}(x, ..., which = c("biplot", "map", "uni", "parallel"), 
   choice = 1:2, coord = NULL, map = NULL, onlyout = TRUE, bw = FALSE, symb = TRUE, 
   symbtxt = FALSE, col = NULL, pch = NULL, obj.cex = NULL, transp = 1)
}
\arguments{
  \item{x}{resulting object from function \code{mvoutlier.CoDa} }
  \item{\dots}{further plotting arguments }
  \item{which}{type of plot that should be made }
  \item{choice}{select the pair of PCs used for the biplot }
  \item{coord}{coordinates for the presentation in a map }
  \item{map}{coordinates for the background map; see details below }
  \item{onlyout}{if TRUE only the outliers are shown in the plot }
  \item{bw}{if TRUE symbold will be in grey scale rather than in color }
  \item{symb}{if TRUE special symbols are used according to outlyingness }
  \item{symbtxt}{if TRUE text labels are used for plotting }
  \item{col}{define colors to be used for outliers and non-outliers }
  \item{pch}{define plotting symbols to be used for outliers and non-outliers }
  \item{obj.cex}{define symbol size for outliers and non-outliers  }
  \item{transp}{define transparancy for parallel coordinate plot }
}
\details{
The function \code{mvoutlier.CoDa} prepares the information needed for this plot function:
In a first step, the raw compositional data set in transformed by the isometric logratio
(ilr) transformation to the usual Euclidean space. Then adaptive outlier detection is
perfomed: Starting from a quantile 1-alpha of the chisquare distribution, one looks for the
supremum of the differences between the chisquare distribution and the empirical distribution
of the squared Mahalanobis distances. The latter are derived from the MCD estimator using
the proportion quan of the data. The supremum is the outlier cutoff, and certain colors
and symbols for the outliers are computed: The colors should reflect the magnitude of the
median element concentration of the observations, which is done by computing for each observation
along the single ilr variables the distances to the medians. The mediab of all distances
determines the color (or grey scale): a high value, resulting in a red (or dark) symbol,
means that most univariate parts have higher values than the average, and a low value (blue
or light symbol) refers to an observation with mainly low values. The symbols are according
to the cut-points from the quantiles 0.25, 0.5, 0.75, and the outlier cutoff of the
squared Mahalanobis distances. This plot function then allows to visualize the information.

The optional background map for the representation of the outliers in a map can be
included using the argument \code{map}. This should consist of one or more polygons
representing the geographical x- and y-coordinates of the background map. Of course,
this map should be represented in the same coordinate system as the coordinates for
the sample locations provided by \code{coord}. The structure of \code{map} is as follows:
It should consist of 2 columns, one for the x-, one for the y-coordinates. If a polygon
ends, a row with 2 entries NA should follow. At the end two NA rows are needed.
See also examples below.
}
\value{
A plot is drawn.
}
\references{
P. Filzmoser, K. Hron, and C. Reimann.
Interpretation of multivariate outliers for compositional data.
Submitted to Computers and Geosciences.
}
\author{Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://cstat.tuwien.ac.at/filz/}}
\seealso{ \code{\link{mvoutlier.CoDa}}, \code{\link{arw}}, \code{\link{map.plot}}, \code{\link{uni.plot}} }
\examples{
data(humus)
el=c("As","Cd","Co","Cu","Mg","Pb","Zn")
dsel <- humus[,el]
data(kola.background) # contains different information (coast, borders, etc.)
coo <- rbind(kola.background$coast,kola.background$boundary,kola.background$borders)
XY <- humus[,c("XCOO","YCOO")]
set.seed(123)
res <- mvoutlier.CoDa(dsel)

par(ask=TRUE)
### Parallel coordinate plot:
## show for all obvervations (transp is only useful when generating e.g. a pdf):
# plot(res,onlyout=FALSE,bw=TRUE,which="parallel",symb=FALSE,symbtxt=FALSE,transp=0.3)
## show only outliers with special colors and labels in the margins:
plot(res,onlyout=TRUE,bw=FALSE,which="parallel",symb=TRUE,symbtxt=TRUE,transp=0.3)

### Biplot:
## show all data points, outliers are in different color and have different symbol:
# plot(res,onlyout=FALSE,which="biplot",bw=FALSE,symb=FALSE,symbtxt=FALSE)
## show only the outliers with special symbols and colors:
plot(res,onlyout=TRUE,which="biplot",bw=FALSE,symb=TRUE,symbtxt=TRUE)

### Map:
## show all data points, outliers are in different color and have different symbol:
# plot(res,coord=XY,map=coo,onlyout=FALSE,which="map",bw=FALSE,symb=FALSE,symbtxt=FALSE)
## show only the outliers with special symbols and colors:
plot(res,coord=XY,map=coo,onlyout=TRUE,which="map",bw=FALSE,symb=TRUE,symbtxt=TRUE)

### Univariate scatterplot:
## show all data points, outliers are in different color and have different symbol:
# plot(res,onlyout=FALSE,which="uni",symb=FALSE,symbtxt=FALSE)
## show only the outliers with special symbols and colors:
plot(res,onlyout=TRUE,which="uni",symb=TRUE,symbtxt=TRUE)
}
\keyword{multivariate}
\keyword{robust}

