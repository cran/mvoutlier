\name{aq.plot}
\alias{aq.plot}
\title{Adjusted Quantile Plot}
\description{
The function aq.plot plots the ordered squared robust Mahalanobis
distances of the observations against
the empirical distribution function of the $MD^2_i$. In addition the distribution function of 
$chisq_p$ is plotted as well as two vertical lines corresponding to the chisq-quantile 
specified in the argument list (default is 0.975) and the so-called adjusted quantile.
Three additional graphics are created (the first showing the data, the second
showing the outliers detected by the specified quantile of the $chisq_p$ distribution and the
third showing these detected outliers by the adjusted quantile).}
\usage{aq.plot(x, delta=qchisq(0.975, df=ncol(x)), quan=1/2, alpha=0.025)}
\arguments{
\item{x}{matrix or data.frame containing the data; has to be at least two-dimensional}
\item{delta}{quantile of the chi-squared distribution with ncol(x) degrees of freedom. This
quantile appears as cyan-colored vertical line in the plot.}
\item{quan}{proportion of observations which are used for mcd estimations; 
has to be between 0.5 and 1, default ist 0.5}
\item{alpha}{Maximum thresholding proportion (optional scalar, default: alpha = 0.025)}
}
\value{
\item{outliers}{boolean vector of outliers}
}
\details{The function aq.plot plots the ordered squared robust Mahalanobis
distances of the observations 
against the empirical distribution function of the $MD^2_i$. The distance calculations are 
based on the MCD estimator. 

For outlier detection two different methods are used. The first one marks observations as 
outliers if they exceed a certain quantile of the chi-squared distribution. The second 
is an adaptive procedure searching for outliers specifically in the tails of the 
distribution, beginning at a certain chisq-quantile (see Filzmoser et al., 2005). 

The function behaves differently depending on the dimension of the data. If the data is 
more than two-dimensional the data are projected on the first two robust
principal components.
}
\author{Moritz Gschwandtner <\email{e0125439@student.tuwien.ac.at}>
\cr Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\references{P. Filzmoser, R.G. Garrett, and C. Reimann.
Multivariate outlier detection in exploration geochemistry.
\emph{Computers & Geosciences}, 31:579-587, 2005.}

\examples{
# create data:
x <- cbind(rnorm(100), rnorm(100), rnorm(100))
y <- cbind(rnorm(10, 5, 1), rnorm(10, 5, 1), rnorm(10, 5, 1))
z <- rbind(x,y)
# execute:
aq.plot(z, alpha=0.1)
}
\keyword{dplot}
