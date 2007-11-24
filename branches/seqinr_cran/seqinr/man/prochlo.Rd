\name{prochlo}
\alias{prochlo}
\docType{data}
\title{ Zscore on three strains of Prochlorococcus marinus }
\encoding{latin1}
\description{
This dataset contains the zscores computed with the codon model on all
CDS from 3 strains of Procholorococcus marinus (as retrieved from Genome
Reviews database on June 16, 2005)
}
\usage{
data(prochlo)
}
\format{
  List of three dataframes of the zscore of each of the 16 dinucleotides
  on each CDS retrieved from the specific strain.

\describe{
  \item{BX548174}{strain adapted to living at a depth of 5
  meters (high levels of UV exposure)
  \code{base} model on each intergenic sequence}
  \item{AE017126}{strain adapted to living at a depth of 120
  meters (low levels of UV exposure)}
  \item{BX548175}{strain adapted to living at a depth of 135
  meters (low levels of UV exposure)}
}
}
\references{
  Palmeira, L., Guéguen, L. and Lobry JR. (2006) UV-targeted dinucleotides
  are not depleted in light-exposed Prokaryotic genomes.
  \emph{Molecular Biology and Evolution},
  \bold{23}:2214-2219.\cr
  \url{http://mbe.oxfordjournals.org/cgi/reprint/23/11/2214}\cr\cr
  
  \code{citation("seqinr")}
}
\seealso{ \code{\link{zscore}} }
\examples{
data(prochlo)
oneplot <- function(x){
  plot(density(prochlo$BX548174[, x]),
    ylim = c(0,0.4), xlim = c(-4,4), lty=3,
    main = paste(substr(x,1,1), "p", substr(x,2,2), " bias", sep = ""),
    xlab="",ylab="",las=1, type = "n")
  rect(-10,-1,-1.96,10, col = "yellow", border = "yellow")
  rect(1.96,-1,10,10, col = "yellow", border = "yellow")
  lines(density(prochlo$BX548174[, x]),lty=3)
  lines(density(prochlo$AE017126[, x]),lty=2)
  lines(density(prochlo$BX548175[, x]),lty=1)
  abline(v=c(-1.96,1.96),lty=5)
  box()
}
par(mfrow=c(2,2),mar=c(2,3,2,0.5) + 0.1)
oneplot("CT")
oneplot("TC")
oneplot("CC")
oneplot("TT")
}
\keyword{ datasets }
