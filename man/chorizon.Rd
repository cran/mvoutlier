\name{chorizon}
\alias{chorizon}
\docType{data}
\title{C-horizon of the Kola Data}
\description{
The Kola Data were collected in the Kola Project (1993-1998, 
Geological Surveys of Finland (GTK) and Norway (NGU) and Central
Kola Expedition (CKE), Russia). More than 600 samples in five
different layers were analysed, this dataset contains the C-horizon.
}
\usage{data(chorizon)}
\format{
  A data frame with 606 observations on the following 110 variables.
  \describe{
    \item{*ID}{a numeric vector}
    \item{XCOO}{a numeric vector}
    \item{YCOO}{a numeric vector}
    \item{Ag}{a numeric vector}
    \item{Ag\_INAA}{a numeric vector}
    \item{Al}{a numeric vector}
    \item{Al2O3}{a numeric vector}
    \item{As}{a numeric vector}
    \item{As\_INAA}{a numeric vector}
    \item{Au\_INAA}{a numeric vector}
    \item{B}{a numeric vector}
    \item{Ba}{a numeric vector}
    \item{Ba\_INAA}{a numeric vector}
    \item{Be}{a numeric vector}
    \item{Bi}{a numeric vector}
    \item{Br\_IC}{a numeric vector}
    \item{Br\_INAA}{a numeric vector}
    \item{Ca}{a numeric vector}
    \item{Ca\_INAA}{a numeric vector}
    \item{CaO}{a numeric vector}
    \item{Cd}{a numeric vector}
    \item{Ce\_INAA}{a numeric vector}
    \item{Cl\_IC}{a numeric vector}
    \item{Co}{a numeric vector}
    \item{Co\_INAA}{a numeric vector}
    \item{EC}{a numeric vector}
    \item{Cr}{a numeric vector}
    \item{Cr\_INAA}{a numeric vector}
    \item{Cs\_INAA}{a numeric vector}
    \item{Cu}{a numeric vector}
    \item{Eu\_INAA}{a numeric vector}
    \item{F\_IC}{a numeric vector}
    \item{Fe}{a numeric vector}
    \item{Fe\_INAA}{a numeric vector}
    \item{Fe2O3}{a numeric vector}
    \item{Hf\_INAA}{a numeric vector}
    \item{Hg}{a numeric vector}
    \item{Hg\_INAA}{a numeric vector}
    \item{Ir\_INAA}{a numeric vector}
    \item{K}{a numeric vector}
    \item{K2O}{a numeric vector}
    \item{La}{a numeric vector}
    \item{La\_INAA}{a numeric vector}
    \item{Li}{a numeric vector}
    \item{LOI}{a numeric vector}
    \item{Lu\_INAA}{a numeric vector}
    \item{wt\_INAA}{a numeric vector}
    \item{Mg}{a numeric vector}
    \item{MgO}{a numeric vector}
    \item{Mn}{a numeric vector}
    \item{MnO}{a numeric vector}
    \item{Mo}{a numeric vector}
    \item{Mo\_INAA}{a numeric vector}
    \item{Na}{a numeric vector}
    \item{Na\_INAA}{a numeric vector}
    \item{Na2O}{a numeric vector}
    \item{Nd\_INAA}{a numeric vector}
    \item{Ni}{a numeric vector}
    \item{Ni\_INAA}{a numeric vector}
    \item{NO3\_IC}{a numeric vector}
    \item{P}{a numeric vector}
    \item{P2O5}{a numeric vector}
    \item{Pb}{a numeric vector}
    \item{pH}{a numeric vector}
    \item{PO4\_IC}{a numeric vector}
    \item{Rb}{a numeric vector}
    \item{S}{a numeric vector}
    \item{Sb}{a numeric vector}
    \item{Sb\_INAA}{a numeric vector}
    \item{Sc}{a numeric vector}
    \item{Sc\_INAA}{a numeric vector}
    \item{Se}{a numeric vector}
    \item{Se\_INAA}{a numeric vector}
    \item{Si}{a numeric vector}
    \item{SiO2}{a numeric vector}
    \item{Sm\_INAA}{a numeric vector}
    \item{Sn\_INAA}{a numeric vector}
    \item{SO4\_IC}{a numeric vector}
    \item{Sr}{a numeric vector}
    \item{Sr\_INAA}{a numeric vector}
    \item{SUM\_XRF}{a numeric vector}
    \item{Ta\_INAA}{a numeric vector}
    \item{Tb\_INAA}{a numeric vector}
    \item{Te}{a numeric vector}
    \item{Th}{a numeric vector}
    \item{Th\_INAA}{a numeric vector}
    \item{Ti}{a numeric vector}
    \item{TiO2}{a numeric vector}
    \item{U\_INAA}{a numeric vector}
    \item{V}{a numeric vector}
    \item{W\_INAA}{a numeric vector}
    \item{Y}{a numeric vector}
    \item{Yb\_INAA}{a numeric vector}
    \item{Zn}{a numeric vector}
    \item{Zn\_INAA}{a numeric vector}
    \item{ELEV}{a numeric vector}
    \item{*COUN}{a numeric vector}
    \item{*ASP}{a numeric vector}
    \item{TOPC}{a numeric vector}
    \item{LITO}{a numeric vector}
    \item{Al\_XRF}{a numeric vector}
    \item{Ca\_XRF}{a numeric vector}
    \item{Fe\_XRF}{a numeric vector}
    \item{K\_XRF}{a numeric vector}
    \item{Mg\_XRF}{a numeric vector}
    \item{Mn\_XRF}{a numeric vector}
    \item{Na\_XRF}{a numeric vector}
    \item{P\_XRF}{a numeric vector}
    \item{Si\_XRF}{a numeric vector}
    \item{Ti\_XRF}{a numeric vector}
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
data(chorizon)
# classical versus robust correlation
cor.plot(log(chorizon[,"Al"]), log(chorizon[,"Na"]))
}
\keyword{datasets}
