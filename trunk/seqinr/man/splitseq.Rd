\name{splitseq}
\alias{splitseq}
\title{ split a sequence in sub-sequences }
\description{
 Split a sequence into sub-sequences of 3 (the default size) with no overlap between the sub-sequences. 
}
\usage{
splitseq(seq, frame = 0, word = 3)
}
\arguments{
  \item{seq}{ an object of class \code{seq} }
  \item{frame}{ an integer giving the frame to split the sequence }
  \item{word}{ an integer giving the size of the sub-sequences }
}

\value{
  This function returns a vector which contains the sub-sequences.
}

\author{ D. CHARIF }

\seealso{ \code{\link{split}} }

\examples{

a=s2c("aacgttgcqqqggtcgctqgctqcgtqgctqctgttt")
## To obtain the codon sequence in frame 0:
splitseq(a)

}
\keyword{ manip }

