\name{uco}
\alias{uco}
\title{ Create codon usage table }
\description{
  \code{uco} calculates some codon usage index: the codon count \code{eff}, the absolute frequencies \code{freq} or the Relative Synonymous Codon Usage \code{rscu}.
}
\usage{
uco( seq, frame = 0, index = c("eff","freq","rscu"), as.data.frame = FALSE) 
}

\arguments{
  \item{seq}{ a vector of chars }
  \item{frame}{ an integer (0,1,2) giving the frame }
  \item{index}{ One of this value }
  tabular{lll}{
	eff  \tab To get the codon count\cr
	freq \tab To get the absolute frequence of codon\cr
	rscu  \tab To get the RSCU indice
	}
  \item{as.data.frame}{ logical. If \code{TRUE}: The eff, freq and rscu index will be returned in a data frame.}
}

\value{
  
  If \code{as.data.frame} is TRUE \code{uco} returns a data frame with three columns: 
  \item{ aa }{ a factor containing the name of amino-acid }
  \item{ codon }{ a factor containing the corresponding codon }
  \item{ eff/freq }{ a factor containing effectifs or frequencies of codons }
  If \code{as.data.frame} and \code{RSCU} are TRUE, \code{uco} returns a data frame with four columns:
  \item{ aa }{ a factor containing the name of amino-acid }
  \item{ codon }{ a factor containing the corresponding codon }
  \item{ eff/freq }{ a factor containing effectifs or frequencies of codons }
  \item{ RSCU}{ a numeric vector containing the relative synonymous codon usage}
  If \code{as.data.frame} is FALSE, it returns a factor whose levels are codons.
	
}

\author{ D. CHARIF }

\examples{

s = read.fasta(File=system.file("sequences/malM.fasta",package = "seqinr"))
uco( s[[1]], index="freq")
uco( s[[1]], index="eff")
uco( s[[1]], index="rscu")
uco( s[[1]], as.data.frame = TRUE)

}

\keyword{ manip }

