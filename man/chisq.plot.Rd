\name{chisq.plot}
\alias{chisq.plot}
\title{Chi-Square Plot}
\description{
The function chisq.plot plots the ordered robust mahalanobis distances of the data against
the quantiles of the Chi-squared distribution. By user interaction this plotting 
is iterated each time leaving out the observation with the greatest distance.}
\usage{chisq.plot(x, quan=1/2, ask=TRUE, ...)}
\arguments{
\item{x}{matrix or data.frame containing the data}
\item{quan}{amount of observations which are used for mcd estimations. 
has to be between 0.5 and 1, default ist 0.5}
\item{ask}{logical. specifies whether user interacton is allowed or not. default is TRUE}
\item{...}{additional graphical parameters}
}
\value{
\item{outliers}{Each time an observation is dropped, it´s index is added to a vector
named outliers, which is accessible within the global environment after cancellation of 
the function.}
}
\details{The function chisq.plot plots the ordered robust mahalanobis distances of the data 
against the quantiles of the Chi-squared distribution. If the data is normal distributed 
these values should approximately correspond to each other, so outliers can be detected 
visually. By user interaction this procedure is repeated, each time leaving out the 
observation with the greatest distance (the number of the observation is printed on the 
console). This method can be seen as an iterativ deletion of outliers until a straight 
line appears. 
}
\author{Moritz Gschwandtner <\email{e0125439@student.tuwien.ac.at}>
\cr Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\references{R.G. Garrett (1989).
The chi-square plot: a tools for multivariate outlier recognition.
\emph{Journal of Geochemical Exploration}, 32 (1/3), 319-341.}
\examples{
data(humus)
chisq.plot(log(humus[,c("Co","Cu","Ni")]))
# The identified outliers are in object "outliers"
}

\keyword{dplot}
