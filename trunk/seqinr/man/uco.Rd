\name{uco}
\alias{uco}
\title{ Create codon usage table }
\description{
  \code{uco} calculates a codon frequencies or count table.
}
\usage{
uco( seq, frame = 0, freq = FALSE, as.data.frame = FALSE)
}

\arguments{
  \item{seq}{ an object of class \code{seq} }
  \item{frame}{ an integer (0,1,2) giving the frame }
  \item{freq}{ a logical value indicating whether the absolute frequence of codon must be return. 
	       If \code{FALSE} you will just have the effectif of the codon.}
  \item{as.data.frame}{ logical. If \code{TRUE}: Effectifs or frequences of codon will be in a data frame.}
}

\value{
  
  If \code{as.data.frame} is TRUE \code{uco} returns a data frame with three columns:
  \item{ aa }{ a factor containing the name of amino-acid }
  \item{ codon }{ a factor containing the corresponding codon }
  \item{ eff/freq }{ a factor containing effectifs or frequencies of codons }
  If \code{as.data.frame} is FALSE, it returns a factor whose levels are codons.
	
}

\author{ D. CHARIF }

\examples{

a=s2c("aaaccgtgtacgtagtcgaa")
uco( a, as.data.frame = T)

}

\keyword{ manip }

