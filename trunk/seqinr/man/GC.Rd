\name{GC}
\alias{GC}
\title{content in G+C from DNA sequences}
\description{
  This function computes the percentage of G+C in a DNA sequence.
}
\usage{
GC(seq)
}
\arguments{
  \item{seq}{an object containig a DNA sequence}
}
\value{
  It returns a numeric vector
}
\author{D. Charif }
\examples{
   s=s2c("agtctggggggccccttttaagtagatagatagctagtcgta")
   GC(s)
}
\keyword{manip}