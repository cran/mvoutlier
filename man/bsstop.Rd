\name{bsstop}
\alias{bsstop}
\docType{data}
\title{Top Layer of the BSS Data}
\description{
The BSS data were collected in agrigultural soils from Northern
Europe. from an area of about 1,800,000 km2. 
769 samples on an iregular grid were taken in two different
layers, the top layer (0-20cm) and the bottom layer.
This dataset contains the top layer of the BSS data. 
It has 46 variables, including x and y coordinates.
}
\usage{data(bsstop)}
\format{
  A data frame with 768 observations on the following 46 variables.
  \describe{
    \item{*ID}{a numeric vector}
    \item{*CNo}{a numeric vector}
    \item{XCOO}{x coordinates: a numeric vector}
    \item{YCOO}{y coordinates: a numeric vector}
    \item{SiO2\_T}{a numeric vector}
    \item{TiO2\_T}{a numeric vector}
    \item{Al2O3\_T}{a numeric vector}
    \item{Fe2O3\_T}{a numeric vector}
    \item{MnO\_T}{a numeric vector}
    \item{MgO\_T}{a numeric vector}
    \item{CaO\_T}{a numeric vector}
    \item{Na2O\_T}{a numeric vector}
    \item{K2O\_T}{a numeric vector}
    \item{P2O5\_T}{a numeric vector}
    \item{(SO3)\_T}{a numeric vector}
    \item{(Cl)\_T}{a numeric vector}
    \item{(F)\_T}{a numeric vector}
    \item{LOI\_T}{a numeric vector}
    \item{(As)\_T}{a numeric vector}
    \item{Ba\_T}{a numeric vector}
    \item{Bi\_T}{a numeric vector}
    \item{Ce\_T}{a numeric vector}
    \item{Co\_T}{a numeric vector}
    \item{Cr\_T}{a numeric vector}
    \item{Cs\_T}{a numeric vector}
    \item{Cu\_T}{a numeric vector}
    \item{Ga\_T}{a numeric vector}
    \item{Hf\_T}{a numeric vector}
    \item{La\_T}{a numeric vector}
    \item{Mo\_T}{a numeric vector}
    \item{Nb\_T}{a numeric vector}
    \item{Ni\_T}{a numeric vector}
    \item{Pb\_T}{a numeric vector}
    \item{Rb\_T}{a numeric vector}
    \item{Sb\_T}{a numeric vector}
    \item{Sc\_T}{a numeric vector}
    \item{Sn\_T}{a numeric vector}
    \item{Sr\_T}{a numeric vector}
    \item{Ta\_T}{a numeric vector}
    \item{Th\_T}{a numeric vector}
    \item{U\_T}{a numeric vector}
    \item{V\_T}{a numeric vector}
    \item{W\_T}{a numeric vector}
    \item{Y\_T}{a numeric vector}
    \item{Zn\_T}{a numeric vector}
    \item{Zr\_T}{a numeric vector}

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
data(bsstop)
# classical versus robust correlation
cor.plot(log(bsstop[, "Al2O3_T"]), log(bsstop[, "Na2O_T"]))
}
\keyword{datasets}
