\name{SeqFastadna}
\alias{SeqFastadna}
\alias{is.SeqFastadna}
\alias{as.SeqFastadna}
\alias{summary.SeqFastadna}
\title{DNA sequence in Fasta Format}
\description{
  \code{as.SeqFastadna} is called by many functions as \code{read.fasta} but not directly by the user. It creates an object of class \code{SeqFastadna}. 
  \code{is.SeqFastadna} returns TRUE if the object is of class \code{SeqFastadna}. 
  \code{summary.SeqFastadna} gives the base composition of an object of class \code{SeqFastadna}.
}
\usage{
as.SeqFastadna(object)
is.SeqFastadna(object)
summary.SeqFastadna(object,...)
}
\arguments{
  \item{object}{ a vector of chars representing a biological sequence} 
}
\value{
  \code{as.SeqFastadna} returns an object sequence of class \code{SeqFastadna}.
  \code{summary.SeqFastadna} returns a list which the following components:
	\item{composition}{ the base counting of the sequence} 
	\item{GC}{ the percentage of G+C in the sequence }
}
\author{ D. Charif }
\examples{
  s = read.fasta("sequences/malM.fasta")
  is.SeqFastadna(s[[1]])
  summary(s[[1]])
}
\keyword{ utilities }
