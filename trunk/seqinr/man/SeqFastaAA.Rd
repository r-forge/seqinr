\name{SeqFastaAA}
\alias{SeqFastaAA}
\alias{is.SeqFastaAA}
\alias{as.SeqFastaAA}
\alias{summary.SeqFastaAA}
\title{AA sequence in Fasta Format}
\description{
  \code{as.SeqFastaAA} is called by many functions as \code{read.fasta} but not directly by the user. It creates an object of class \code{SeqFastaAA}. 
  \code{is.SeqFastaAA} returns TRUE if the object is of class \code{SeqFastaAA}. 
  \code{summary.SeqFastaAA} gives the AA composition of an object of class \code{SeqFastaAA}.
}
\usage{
as.SeqFastaAA(object)
is.SeqFastaAA(object)
summary.SeqFastaAA(object,...)
}
\arguments{
  \item{object}{ a vector of chars representing a biological sequence}
}
\value{
  \code{as.SeqFastaAA} returns an object sequence of class \code{SeqFastaAA}.
  \code{summary.SeqFastaAA} returns a list which the following components:
	\item{composition}{ the AA counting of the sequence} 
	\item{GC}{ the percentage of each group of amino acid in the sequence. By example, the groups are small, tiny, aliphatic, aromatic ... }
}
\author{ D. Charif }
\examples{
  \dontrun{s = read.fasta(File=system.file("sequences/seqAA.fasta",package="seqinr"),seqtype="AA")}
  \dontrun{is.SeqFastaAA(s[[1]])}
  \dontrun{summary(s[[1]])}
}
\keyword{ utilities }
