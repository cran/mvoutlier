\name{pbb}
\alias{pbb}
\title{BSS background Plot}
\description{
Plots the BSS background map
}
\usage{
pbb(map = "bss.background", add.plot = FALSE, ...)
}
\arguments{
  \item{map}{List of coordinates. For the correct format see also help(kola.background)}
  \item{add.plot}{logical. If true background is added to an existing plot}
  \item{\dots}{additional plot parameters, see help(par)}
}
\details{
  The list of coordinates is plotted as a polygon line.
}
\value{
  The plot is produced on the graphical device.
}
\references{
Reimann C, Siewers U, Tarvainen T, Bityukova L, Eriksson J, Gilucis A,
Gregorauskiene V, Lukashev VK, Matinian NN, Pasieczna A.
Agricultural Soils in Northern Europe: A Geochemical Atlas.
Geologisches Jahrbuch, Sonderhefte, Reihe D, Heft SD 5,
Schweizerbart'sche Verlagsbuchhandlung, Stuttgart, 2003.
}
\author{Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\seealso{See also \code{\link{pkb}}}
\examples{
data(bss.background)
data(bsstop)
plot(bsstop$XCOO,bsstop$YCOO,col="red",pch=3)
pbb(add=TRUE)
}
\keyword{datasets}
