\name{plot.mvoutlierCoDa}
\alias{plot.mvoutlierCoDa}
\title{Plots for interpreting multivatiate outliers of CoDa}
\description{
Plots the computed information by \code{mvoutlier.CoDa} for supporting the interpretation
of multivariate outliers in case of compositional data.}
\usage{
\method{plot}{mvoutlierCoDa}(x, ..., which = c("biplot", "map", "uni", "parallel"), choice = 1:2, coord = NULL, map = NULL, onlyout = TRUE, bw = FALSE, symb = TRUE, symbtxt = FALSE, col = NULL, pch = NULL, obj.cex = NULL, transp = 1)
}
\arguments{
  \item{x}{resulting object from function \code{mvoutlier.CoDa} }
  \item{\dots}{further plotting arguments }
  \item{which}{type of plot that should be made }
  \item{choice}{select the pair of PCs used for the biplot }
  \item{coord}{coordinates for the presentation in a map }
  \item{map}{coordinates for the background map }
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
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\seealso{ code{\link{mvoutlier.CoDa}}, \code{\link{arw}}, \code{\link{map.plot}}, \code{\link{uni.plot}} }
\examples{
data(humus)
d <- humus[,c("As","Cd","Co","Cu","Mg","Pb","Zn")]
res <- mvoutlier.CoDa(d)
plot(res,onlyout=FALSE,bw=FALSE,which="uni",symb=TRUE,symbtxt=FALSE)
}
\keyword{multivariate}
\keyword{robust}

