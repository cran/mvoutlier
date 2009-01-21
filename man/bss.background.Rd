\name{bss.background}
\alias{bss.background}
\docType{data}
\title{Background map for the BSS project}
\description{
Coordinates of the BSS data background map
}
\usage{data(bss.background)}
\format{
  A data frame with 6093 observations on the following 2 variables.
  \describe{
    \item{\code{V1}}{a numeric vector with the x-coordinates}
    \item{\code{V2}}{a numeric vector with the y-coordinates}
  }
}
\details{
Is used by pbb()
}
\source{
BSS project
}
\references{
Reimann C, Siewers U, Tarvainen T, Bityukova L, Eriksson J, Gilucis A,
Gregorauskiene V, Lukashev VK, Matinian NN, Pasieczna A.
Agricultural Soils in Northern Europe: A Geochemical Atlas.
Geologisches Jahrbuch, Sonderhefte, Reihe D, Heft SD 5,
Schweizerbart'sche Verlagsbuchhandlung, Stuttgart, 2003.
}
\examples{
data(bss.background)
pbb()
}
\keyword{datasets}
