\name{count}
\alias{count}
\title{Composition of dimer/trimer/etc nucleotides}
\description{
  Counts the number of times dimer/trimer/etc nucleotides occurs in a sequence.
}
\usage{
count(seq,word,frame=0)
}
\arguments{
  \item{seq}{ an object of class \code{seq} }
  \item{word}{ an integer giving the size of word (n-mer) to count }
  \item{frame}{	an integer (0,1,2) giving the frame }
}
\details{
 \code{count} counts the occurence of all words by moving a window of length \code{word} up by one each time.
 You can count only the words that occur in a single frame of \code{seq} by setting \code{frame} to a number other than zero.
 If you set it to 1 it will only count the words in frame 1, 2 will only count the words in frame 2 and so on.
}
\value{
  This function returns a factor whose levels are all the possible dimer/trimer/etc nucleotides. The word with zero count will be display. 		
}
\author{ D. CHARIF }

\seealso{ \code{\link{table}} }

\examples{

a=s2c("acgggtqcggtcccqtcgqt")
##To count the occurences of dinucleotides of a:
count(a,2)
##To count the occurences of trinucleotides of a in frame 2:
count(a,3,2)

}

\keyword{ manip }

