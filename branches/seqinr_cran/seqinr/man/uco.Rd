\name{uco}
\alias{uco}
\title{ Codon usage indices }
\description{
  \code{uco} calculates some codon usage indices: the codon counts \code{eff}, the relative frequencies \code{freq} or the Relative Synonymous Codon Usage \code{rscu}.
}
\usage{
uco(seq, frame = 0, index = c("eff", "freq", "rscu"), as.data.frame = FALSE) 
}
\arguments{
  \item{seq}{ a coding sequence as a vector of chars }
  \item{frame}{ an integer (0, 1, 2) giving the frame of the coding sequence }
  \item{index}{ codon usage index choice, partial matching is allowed. 
                \code{eff} for codon counts, 
                \code{freq} for codon relative frequencies, 
                and \code{rscu} the RSCU index}
  \item{as.data.frame}{ logical. If \code{TRUE}: all indices are returned into a data frame.}
}
\details{
  Codons with ambiguous bases are ignored.
}
\value{
    If \code{as.data.frame} is TRUE \code{uco} returns a data frame with five columns: 
  \item{ aa }{ a vector containing the name of amino-acid }
  \item{ codon }{ a vector containing the corresponding codon }
  \item{ eff }{ a numeric vector of codon counts }
  \item{ freq }{ a numeric vector of codon relative frequencies }
  \item{ rscu }{ a numeric vector of RSCU index }
  If \code{as.data.frame} is FALSE, \code{uco} returns one of this index: 
  \item{ eff }{ a table of codon counts }
  \item{ freq }{ a table of codon relative frequencies }	
  \item{ rscu }{ a vector of relative synonymous codon usage values}
}
\references{
\code{citation("seqinr")} \cr
}
\author{ D. Charif, J.R. Lobry }
\examples{
## Show all possible codons:
words()
## Make a coding sequence from this:
(cds <- s2c(paste(words(), collapse = "")))
## Get codon counts:
uco(cds, index = "eff")
## Get codon relative frequencies:
uco(cds, index = "freq")
## Get RSCU values:
uco(cds, index = "rscu")
## Show what's happen with ambiguous bases:
uco(s2c("aaannnttt"))
## Use a real coding sequence:
rcds <- read.fasta(File = system.file("sequences/malM.fasta", package = "seqinr"))[[1]]
uco( rcds, index = "freq")
uco( rcds, index = "eff")
uco( rcds, index = "rscu")
uco( rcds, as.data.frame = TRUE)
}
\keyword{ manip }
