\encoding{latin1}
\name{bhorizon}
\alias{bhorizon}
\docType{data}
\title{B-horizon of the Kola Data}
\description{
The Kola data were collected in the Kola Project (1993-1998, 
Geological Surveys of Finland (GTK) and Norway (NGU) and Central
Kola Expedition (CKE), Russia). More than 600 samples in five
different layers were analysed, this dataset contains the B-horizon.
}
\usage{data(bhorizon)}
\format{
  A data frame with 609 observations on the following 48 variables.
  \describe{
    \item{\code{ID}}{a numeric vector}
    \item{\code{XCOO}}{a numeric vector}
    \item{\code{YCOO}}{a numeric vector}
    \item{\code{Ag}}{a numeric vector}
    \item{\code{Al}}{a numeric vector}
    \item{\code{Al_XRF}}{a numeric vector}
    \item{\code{As}}{a numeric vector}
    \item{\code{Ba}}{a numeric vector}
    \item{\code{Be}}{a numeric vector}
    \item{\code{Bi}}{a numeric vector}
    \item{\code{Ca}}{a numeric vector}
    \item{\code{Ca_XRF}}{a numeric vector}
    \item{\code{Cd}}{a numeric vector}
    \item{\code{Co}}{a numeric vector}
    \item{\code{Cr}}{a numeric vector}
    \item{\code{Cu}}{a numeric vector}
    \item{\code{EC}}{a numeric vector}
    \item{\code{Fe}}{a numeric vector}
    \item{\code{Fe_XRF}}{a numeric vector}
    \item{\code{K}}{a numeric vector}
    \item{\code{K_XRF}}{a numeric vector}
    \item{\code{LOI}}{a numeric vector}
    \item{\code{La}}{a numeric vector}
    \item{\code{Li}}{a numeric vector}
    \item{\code{Mg}}{a numeric vector}
    \item{\code{Mg_XRF}}{a numeric vector}
    \item{\code{Mn}}{a numeric vector}
    \item{\code{Mn_XRF}}{a numeric vector}
    \item{\code{Mo}}{a numeric vector}
    \item{\code{Na}}{a numeric vector}
    \item{\code{Na_XRF}}{a numeric vector}
    \item{\code{Ni}}{a numeric vector}
    \item{\code{P}}{a numeric vector}
    \item{\code{P_XRF}}{a numeric vector}
    \item{\code{Pb}}{a numeric vector}
    \item{\code{S}}{a numeric vector}
    \item{\code{Sc}}{a numeric vector}
    \item{\code{Se}}{a numeric vector}
    \item{\code{Si}}{a numeric vector}
    \item{\code{Si_XRF}}{a numeric vector}
    \item{\code{Sr}}{a numeric vector}
    \item{\code{Te}}{a numeric vector}
    \item{\code{Th}}{a numeric vector}
    \item{\code{Ti}}{a numeric vector}
    \item{\code{Ti_XRF}}{a numeric vector}
    \item{\code{V}}{a numeric vector}
    \item{\code{Y}}{a numeric vector}
    \item{\code{Zn}}{a numeric vector}
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
data(bhorizon)
# classical versus robust correlation
corr.plot(log(bhorizon[,"Al"]), log(bhorizon[,"Na"]))
}
\keyword{datasets}
