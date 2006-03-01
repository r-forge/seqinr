\name{G+C Content}
\alias{GC}
\alias{GC2}
\alias{GC3}
\title{Calculates the fractional G+C content of nucleic acid sequences}
\description{
  Calculates the fraction of G+C bases of the input nucleic acid
  sequence(s). It reads in nucleic acid sequences, sums the number of
  'G' and 'C' bases and writes out the result as the fraction (in the
  interval 0.0 to 1.0) of the length of the whole sequence. Global G+C
  content \code{GC}, G+C in the second position of the codon bases
  \code{GC2}, and G+C in the third position of the codon bases
  \code{GC3} can be computed.
}
\usage{
GC(seq)
GC2(seq)
GC3(seq)
}
\arguments{
  \item{seq}{any nucleic acid sequence}
}
\value{
  \code{GC} returns the fraction of G+C as a numeric vector,
  \code{GC2} returns the fraction of GC2 as a numeric vector,
  \code{GC3} returns the fraction of GC3 as a numeric vector.
}
\references{ 
To have an overview of the seqinR's functionnality, please consult this vignette: 
Charif, D., Lobry, J.R. (2005) SeqinR: a contributed package to the R project for statistical
computing devoted to biological sequences retrieval and analysis. Springer Verlag, \emph{Biological and Medical Physics/Biomedical Series}, in preparation. \cr
}
\author{D. Charif }
\examples{
   s=s2c("agtctggggggccccttttaagtagatagatagctagtcgta")
   GC(s)
   GC2(s)
   GC3(s)
}
\keyword{manip}
