\name{bssbot}
\alias{bssbot}
\docType{data}
\title{Bottom Layer of the BSS Data}
\description{
The BSS data were collected in agrigultural soils from Northern
Europe. from an area of about 1,800,000 km2.
769 samples on an iregular grid were taken in two different
layers, the top layer (0-20cm) and the bottom layer.
This dataset contains the bottom layer of the BSS data.
It has 46 variables, including x and y coordinates.
}
\usage{data(bssbot)}
\format{
  A data frame with 768 observations on the following 46 variables.
  \describe{
    \item{ID}{a numeric vector}
    \item{CNo}{a numeric vector}
    \item{XCOO}{x coordinates: a numeric vector}
    \item{YCOO}{y coordinates: a numeric vector}
    \item{SiO2\_B}{a numeric vector}
    \item{TiO2\_B}{a numeric vector}
    \item{Al2O3\_B}{a numeric vector}
    \item{Fe2O3\_B}{a numeric vector}
    \item{MnO\_B}{a numeric vector}
    \item{MgO\_B}{a numeric vector}
    \item{CaO\_B}{a numeric vector}
    \item{Na2O\_B}{a numeric vector}
    \item{K2O\_B}{a numeric vector}
    \item{P2O5\_B}{a numeric vector}
    \item{SO3\_B}{a numeric vector}
    \item{Cl\_B}{a numeric vector}
    \item{F\_B}{a numeric vector}
    \item{LOI\_B}{a numeric vector}
    \item{As\_B}{a numeric vector}
    \item{Ba\_B}{a numeric vector}
    \item{Bi\_B}{a numeric vector}
    \item{Ce\_B}{a numeric vector}
    \item{Co\_B}{a numeric vector}
    \item{Cr\_B}{a numeric vector}
    \item{Cs\_B}{a numeric vector}
    \item{Cu\_B}{a numeric vector}
    \item{Ga\_B}{a numeric vector}
    \item{Hf\_B}{a numeric vector}
    \item{La\_B}{a numeric vector}
    \item{Mo\_B}{a numeric vector}
    \item{Nb\_B}{a numeric vector}
    \item{Ni\_B}{a numeric vector}
    \item{Pb\_B}{a numeric vector}
    \item{Rb\_B}{a numeric vector}
    \item{Sb\_B}{a numeric vector}
    \item{Sc\_B}{a numeric vector}
    \item{Sn\_B}{a numeric vector}
    \item{Sr\_B}{a numeric vector}
    \item{Ta\_B}{a numeric vector}
    \item{Th\_B}{a numeric vector}
    \item{U\_B}{a numeric vector}
    \item{V\_B}{a numeric vector}
    \item{W\_B}{a numeric vector}
    \item{Y\_B}{a numeric vector}
    \item{Zn\_B}{a numeric vector}
    \item{Zr\_B}{a numeric vector}
  }
}
\source{
BSS Project in Northern Europe
}
\references{
Reimann C, Siewers U, Tarvainen T, Bityukova L, Eriksson J, Gilucis A,
Gregorauskiene V, Lukashev VK, Matinian NN, Pasieczna A.
Agricultural Soils in Northern Europe: A Geochemical Atlas.
Geologisches Jahrbuch, Sonderhefte, Reihe D, Heft SD 5,
Schweizerbart'sche Verlagsbuchhandlung, Stuttgart, 2003.
}
\examples{
data(bssbot)
# classical versus robust correlation
cor.plot(log(bssbot[, "Al2O3_B"]), log(bssbot[, "Na2O_B"]))
}

\keyword{datasets}
