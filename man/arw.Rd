\name{arw}
\alias{arw}
\title{Adaptive reweighted estimator for multivariate location and scatter}
\description{
Adaptive reweighted estimator for multivariate location and scatter
with hard-rejection weights. 
The multivariate outliers are defined 
according to the supremum of the difference between the empirical
distribution function of the robust Mahalanobis distance and the 
theoretical distribution function.
}
\usage{
arw(x, m0, c0, alpha, pcrit)
}
\arguments{
  \item{x}{Dataset (n x p)}
  \item{m0}{Initial location estimator (1 x p)}
  \item{c0}{Initial scatter estimator (p x p)}
  \item{alpha}{Maximum thresholding proportion (optional scalar, default: alpha = 0.025)}
  \item{pcrit}{Critical value obtained by simulations (optional scalar, default value
obtained from simulations)}
}
\value{
  \item{m}{Adaptive location estimator (p x 1)}
  \item{c}{Adaptive scatter estimator (p x p)}
  \item{cn}{Adaptive threshold ("adjusted quantile")}
  \item{w}{Weight vector (n x 1)}
}
\details{At the basis of initial estimators of location and scatter, the function arw 
performs a reweighting step to adjust the threshold for outlier rejection. 
The critical value pcrit was obtained by simulations using the MCD estimator
as initial robust covariance estimator. If a different estimator is used,
pcrit should be changed and computed by simulations for the specific dimensions
of the data x.
}

\author{Moritz Gschwandtner <\email{e0125439@student.tuwien.ac.at}>
\cr Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\references{P. Filzmoser, R.G. Garrett, and C. Reimann (2005).
Multivariate outlier detection in exploration geochemistry.
\emph{Computers & Geosciences}. To appear.}

\examples{
x <- cbind(rnorm(100), rnorm(100))
arw(x, apply(x,2,mean), cov(x))
}
\keyword{dplot}
