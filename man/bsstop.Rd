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
    \item{ID}{a numeric vector}
    \item{CNo}{a numeric vector}
    \item{XCOO}{x coordinates: a numeric vector}
    \item{YCOO}{y coordinates: a numeric vector}
    \item{SiO2_T}{a numeric vector}
    \item{TiO2_T}{a numeric vector}
    \item{Al2O3_T}{a numeric vector}
    \item{Fe2O3_T}{a numeric vector}
    \item{MnO_T}{a numeric vector}
    \item{MgO_T}{a numeric vector}
    \item{CaO_T}{a numeric vector}
    \item{Na2O_T}{a numeric vector}
    \item{K2O_T}{a numeric vector}
    \item{P2O5_T}{a numeric vector}
    \item{SO3_T}{a numeric vector}
    \item{Cl_T}{a numeric vector}
    \item{F_T}{a numeric vector}
    \item{LOI_T}{a numeric vector}
    \item{As_T}{a numeric vector}
    \item{Ba_T}{a numeric vector}
    \item{Bi_T}{a numeric vector}
    \item{Ce_T}{a numeric vector}
    \item{Co_T}{a numeric vector}
    \item{Cr_T}{a numeric vector}
    \item{Cs_T}{a numeric vector}
    \item{Cu_T}{a numeric vector}
    \item{Ga_T}{a numeric vector}
    \item{Hf_T}{a numeric vector}
    \item{La_T}{a numeric vector}
    \item{Mo_T}{a numeric vector}
    \item{Nb_T}{a numeric vector}
    \item{Ni_T}{a numeric vector}
    \item{Pb_T}{a numeric vector}
    \item{Rb_T}{a numeric vector}
    \item{Sb_T}{a numeric vector}
    \item{Sc_T}{a numeric vector}
    \item{Sn_T}{a numeric vector}
    \item{Sr_T}{a numeric vector}
    \item{Ta_T}{a numeric vector}
    \item{Th_T}{a numeric vector}
    \item{U_T}{a numeric vector}
    \item{V_T}{a numeric vector}
    \item{W_T}{a numeric vector}
    \item{Y_T}{a numeric vector}
    \item{Zn_T}{a numeric vector}
    \item{Zr_T}{a numeric vector}
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
corr.plot(log(bsstop[, "Al2O3_T"]), log(bsstop[, "Na2O_T"]))
}
\keyword{datasets}
