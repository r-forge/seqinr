\name{GC2}
\alias{GC2}
\title{Calculates the fractional GC2 content of nucleic acid sequences}
\description{
This calculates the fraction of G+C in the second position of the codon bases of the input nucleic acid sequence(s)}
\usage{
GC2(seq)
}
\arguments{
  \item{seq}{any nucleic acid sequence}
}
\value{
  It returns the fraction of GC2 as a numeric vector
}
\author{D. Charif }
\examples{
   s=s2c("agtctggggggccccttttaagtagatagatagctagtcgta")
   GC2(s)
}
\keyword{manip}
