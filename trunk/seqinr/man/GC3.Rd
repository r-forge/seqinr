\name{GC3}
\alias{GC3}
\title{Calculates the fractional GC3 content of nucleic acid sequences}
\description{
This calculates the fraction of G+C  in third position of the codon bases of the input nucleic acid sequence(s)
\usage{
GC3(seq)
}
\arguments{
  \item{seq}{any nucleic acid sequence}
}
\value{
  It returns the fraction of GC3 as a numeric vector
}
\author{D. Charif }
\examples{
   s=s2c("agtctggggggccccttttaagtagatagatagctagtcgta")
   GC3(s)
}
\keyword{manip}
