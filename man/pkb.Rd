\encoding{latin1}
\name{pkb}
\alias{pkb}
\title{Kola background Plot}
\description{
Plots the Kola background map
}
\usage{
pkb(map = "kola.background", which.map = c(1, 2, 3, 4), map.col = c(5, 1, 3, 4), 
    map.lwd = c(2, 1, 2, 1), add.plot = FALSE, ...)
}
\arguments{
  \item{map}{List of coordinates. For the correct format see also help(kola.background)}
  \item{which.map}{which==1 ... plot project boundary \# which==2 ... plot coast line \# which==3 ... plot country borders \# which==4 ... plot lakes and rivers}
  \item{map.col}{Map colors to be used}
  \item{map.lwd}{Defines linestyle of the background}
  \item{add.plot}{logical. if true background is added to an existing plot}
  \item{\dots}{additional plot parameters, see help(par)}
}
\details{
Is used by map.plot()
}
\author{Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}>
\url{http://www.statistik.tuwien.ac.at/public/filz/}}
\references{
Reimann C, Äyräs M, Chekushin V, Bogatyrev I, Boyd R, Caritat P de, Dutter R,
Finne TE, Halleraker JH, Jæger Ø, Kashulina G, Lehto O, Niskavaara H, Pavlov V,
Räisänen ML, Strand T, Volden T. Environmental Geochemical Atlas of the Central
Barents Region. NGU-GTK-CKE Special Publication, Geological Survey of Norway,
Trondheim, Norway, 1998.
}
\examples{
example(map.plot)
}
\keyword{datasets}
