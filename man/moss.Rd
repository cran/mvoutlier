\name{moss}
\alias{moss}
\docType{data}
\title{Moss Layer of the Kola Data}
\description{
The Kola Data were collected in the Kola Project (1993-1998, 
Geological Surveys of Finland (GTK) and Norway (NGU) and Central
Kola Expedition (CKE), Russia). More than 600 samples in five
different layers were analysed, this dataset contains the moss layer.
}
\usage{data(moss)}
\format{
  A data frame with 598 observations on the following 34 variables.
  \describe{
    \item{ID}{a numeric vector}
    \item{XCOO}{a numeric vector}
    \item{YCOO}{a numeric vector}
    \item{Ag}{a numeric vector}
    \item{Al}{a numeric vector}
    \item{As}{a numeric vector}
    \item{B}{a numeric vector}
    \item{Ba}{a numeric vector}
    \item{Bi}{a numeric vector}
    \item{Ca}{a numeric vector}
    \item{Cd}{a numeric vector}
    \item{Co}{a numeric vector}
    \item{Cr}{a numeric vector}
    \item{Cu}{a numeric vector}
    \item{Fe}{a numeric vector}
    \item{Hg}{a numeric vector}
    \item{K}{a numeric vector}
    \item{Mg}{a numeric vector}
    \item{Mn}{a numeric vector}
    \item{Mo}{a numeric vector}
    \item{Na}{a numeric vector}
    \item{Ni}{a numeric vector}
    \item{P}{a numeric vector}
    \item{Pb}{a numeric vector}
    \item{Rb}{a numeric vector}
    \item{S}{a numeric vector}
    \item{Sb}{a numeric vector}
    \item{Si}{a numeric vector}
    \item{Sr}{a numeric vector}
    \item{Th}{a numeric vector}
    \item{Tl}{a numeric vector}
    \item{U}{a numeric vector}
    \item{V}{a numeric vector}
    \item{Zn}{a numeric vector}
  }
}
\source{
Kola Project (1993-1998)
}
\references{
Reimann C, Äyräs M, Chekushin V, Bogatyrev I, Boyd R, Caritat P de, Dutter R,
Finne TE, Halleraker JH, Jæger Ø, Kashulina G, Lehto O, Niskavaara H, Pavlov V,
Räisänen ML, Strand T, Volden T. Environmental Geochemical Atlas of the Central
Barents Region. NGU-GTK-CKE Special Publication, Geological Survey of Norway,
Trondheim, Norway, 1998.
}
\examples{
data(moss)
# classical versus robust correlation:
cor.plot(log(moss[,"Al"]), log(moss[,"Na"]))
}
\keyword{datasets}
