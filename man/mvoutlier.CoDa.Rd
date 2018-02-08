\name{mvoutlier.CoDa}
\alias{mvoutlier.CoDa}
\title{Interpreting multivatiate outliers of CoDa}
\description{
Computes the basis information for plot functions supporting the interpretation
of multivariate outliers in case of compositional data. }
\usage{
mvoutlier.CoDa(x, quan = 0.75, alpha = 0.025, 
   col.quantile = c(0, 0.05, 0.1, 0.5, 0.9, 0.95, 1), 
   symb.pch = c(3, 3, 16, 1, 1), symb.cex = c(1.5, 1, 0.5, 1, 1.5), 
   adaptive = TRUE)
}
\arguments{
  \item{x}{data set (matrix or data frame) containing the raw untransformed compositional data}
  \item{quan}{quantity of data used for robust estimation; between 0.5 and 1}
  \item{alpha}{maximum threshold for adaptive outlier detection}
  \item{col.quantile}{quantiles of an average concentration defining the colors}
  \item{symb.pch}{plotting character for symbols}
  \item{symb.cex}{plotting size for symbols }
  \item{adaptive}{if TRUE then the adaptive method for the outlier threshold is used }
}
\details{
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
squared Mahalanobis distances.
}
\value{
  \item{ilrvariables }{the ilr transformed data matrix}
  \item{outliers }{TRUE/FALSE vector; TRUE refers to outlier}
  \item{pcaobj }{object from PCA}
  \item{colcol }{vector with the colors}
  \item{colbw }{vector with the grey scales}
  \item{pchvec }{vector with plotting symbols}
  \item{cexvec }{vector with sizes of plot symbols}
}
\references{
P. Filzmoser, K. Hron, and C. Reimann.
Interpretation of multivariate outliers for compositional data.
Submitted to Computers and Geosciences.
}
\author{Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://cstat.tuwien.ac.at/filz/}}
\seealso{ \code{\link{plot.mvoutlierCoDa}}, \code{\link{arw}}, \code{\link{map.plot}}, \code{\link{uni.plot}} }
\examples{
data(humus)
d <- humus[,c("As","Cd","Co","Cu","Mg","Pb","Zn")]
res <- mvoutlier.CoDa(d)
str(res)
}
\keyword{multivariate}
\keyword{robust}
