\name{read.mase}
\alias{read.mase}
\title{ Read File in mase format }
\description{
  Read a file containing aligned sequences of proteins in mase format. This format is used to store nucleotide or protein multiple alignments. The beginning of the file must contain a header containing at least one line (but the content of this header may be empty). The header lines must begin by \code{;;}. The body of the file has the following structure: First, each entry must begin by one (or more) commentary line. Commentary lines begin by the character \code{;}. Again, this commentary line may be empty. After the commentaries, the name of the sequence is written on a separate line. At last, the sequence itself is written on the following lines.
}
\usage{
read.mase(File = system.file("sequences/test.mase", package = "seqinr"))
}
\arguments{
  \item{File}{ The name of the file which the aligned sequences in mase format are to be read from. If it does not contain an absolute path, the file name is relative to the current working directory, getwd() }
}
\details{
  By default \code{read.mase}  
}
\value{
  It returns an object of class \code{mase} which is a list with the following components:
  \item{seq.number}{ the number of aligned sequences }
  \item{site.number}{ the length of the alignment}
  \item{sequence}{ a vector of characters containing tha aligned sequences} 
  \item{comment}{ a vector of chars containing the commentaries for each sequence}
  \item{name}{ a vector of chars containing the name of the aligned sequences }
}
\author{ D. CHARIF }
\examples{
s = read.mase(File=system.file("sequences/test.mase", package = "seqinr"))
}
}
\keyword{ manip }
