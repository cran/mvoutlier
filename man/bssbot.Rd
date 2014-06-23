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
    \item{SiO2_B}{a numeric vector}
    \item{TiO2_B}{a numeric vector}
    \item{Al2O3_B}{a numeric vector}
    \item{Fe2O3_B}{a numeric vector}
    \item{MnO_B}{a numeric vector}
    \item{MgO_B}{a numeric vector}
    \item{CaO_B}{a numeric vector}
    \item{Na2O_B}{a numeric vector}
    \item{K2O_B}{a numeric vector}
    \item{P2O5_B}{a numeric vector}
    \item{SO3_B}{a numeric vector}
    \item{Cl_B}{a numeric vector}
    \item{F_B}{a numeric vector}
    \item{LOI_B}{a numeric vector}
    \item{As_B}{a numeric vector}
    \item{Ba_B}{a numeric vector}
    \item{Bi_B}{a numeric vector}
    \item{Ce_B}{a numeric vector}
    \item{Co_B}{a numeric vector}
    \item{Cr_B}{a numeric vector}
    \item{Cs_B}{a numeric vector}
    \item{Cu_B}{a numeric vector}
    \item{Ga_B}{a numeric vector}
    \item{Hf_B}{a numeric vector}
    \item{La_B}{a numeric vector}
    \item{Mo_B}{a numeric vector}
    \item{Nb_B}{a numeric vector}
    \item{Ni_B}{a numeric vector}
    \item{Pb_B}{a numeric vector}
    \item{Rb_B}{a numeric vector}
    \item{Sb_B}{a numeric vector}
    \item{Sc_B}{a numeric vector}
    \item{Sn_B}{a numeric vector}
    \item{Sr_B}{a numeric vector}
    \item{Ta_B}{a numeric vector}
    \item{Th_B}{a numeric vector}
    \item{U_B}{a numeric vector}
    \item{V_B}{a numeric vector}
    \item{W_B}{a numeric vector}
    \item{Y_B}{a numeric vector}
    \item{Zn_B}{a numeric vector}
    \item{Zr_B}{a numeric vector}
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
corr.plot(log(bssbot[, "Al2O3_B"]), log(bssbot[, "Na2O_B"]))
}

\keyword{datasets}
