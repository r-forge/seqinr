\name{SeqFrag}
\alias{SeqFrag}
\alias{is.SeqFrag}
\alias{as.SeqFrag}
\title{Fragment of Sequence}
\description{
  \code{as.SeqFrag} is called by all methods of \code{getFrag}, but not directly by the users. It creates an object sequence of class \code{SeqFrag}. 
}
\usage{
 as.SeqFrag(object, begin, end, compl = F, name = "frag")
 is.SeqFrag(object)
}
\arguments{
  \item{object}{ an object sequence of class \code{seqFastadna}, \code{seqFastaAA}, \code{seqAcnucLocal}, \code{seqAcnucWeb} or \code{seqFrag} }
  \item{begin}{ the first base of the fragment to get }
  \item{end}{ the last base of the fragment to get }
  \item{compl}{ if TRUE, you must give a name for the sequence}
  \item{name}{ the name of the fragment of the sequence }
}
\value{
  \code{as.SeqFrag} returns a biological sequence representing by a vector of chars with the following attributes:
  \item{seqMother}{ the name of the sequence from which the sequence comes}
  \item{begin}{ the position of the first base of the fragment on the mother sequence}
  \item{end}{ the position of the last base of the fragment on the mother sequence }
  \item{class}{ \code{SeqFrag} which is the newest class of the sequence } 
  
  \code{is.SeqFrag} returns TRUE if the object is of class Seqfrag. 
}
\author{ D. Charif }
\seealso{\code{\link{getFrag}} }
\examples{
dontrun{
  choosebank(bank="demo")
  query("ecoli","sp=escherichia coli et t=cds")
  getFrag(ecoli$req[[1]],3,100)	
}
}
\keyword{utilities}
