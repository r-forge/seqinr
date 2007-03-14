\name{G+C Content}
\alias{GC}
\alias{GC1}
\alias{GC2}
\alias{GC3}
\title{Calculates the fractional G+C content of nucleic acid sequences}
\description{
  Calculates the fraction of G+C bases of the input nucleic acid
  sequence(s). It reads in nucleic acid sequences, sums the number of
  'G' and 'C' bases and writes out the result as the fraction (in the
  interval 0.0 to 1.0) of the length of the whole sequence. Global G+C
  content \code{GC}, G+C in the first position of the codon bases
  \code{GC1}, G+C in the second position of the codon bases
  \code{GC2}, and G+C in the third position of the codon bases
  \code{GC3} can be computed.
}
\usage{
GC(seq)
GC1(seq)
GC2(seq)
GC3(seq)
}
\arguments{
  \item{seq}{any nucleic acid sequence as a vector of single characters in lower case}
}
\value{
  \code{GC} returns the fraction of G+C as a numeric vector,
  \code{GC1} returns the fraction of G+C at first codon position as a numeric vector,
  \code{GC2} returns the fraction of G+C at second codon position  as a numeric vector,
  \code{GC3} returns the fraction of G+C at third codon position  as a numeric vector.
}
\references{ 
  \code{citation("seqinr")}
}
\seealso{You can use \code{\link{s2c}} to convert a string into a vetor of single
character and \code{\link{tolower}} to convert upper-case characters into
lower-case characters.
}
\author{D. Charif and L. Palmeira and J.R. Lobry}
\examples{
   mysequence <- s2c("agtctggggggccccttttaagtagatagatagctagtcgta")
   GC(mysequence)   # 0.4761905
   GC1(mysequence) # 0.6428571
   GC2(mysequence) # 0.3571429
   GC3(mysequence) # 0.4285714
#
# With upper-case characters:
#
  myUCsequence <- s2c("GGGGGGGGGA")
  GC(tolower(myUCsequence)) # 0.9
}
\keyword{manip}
