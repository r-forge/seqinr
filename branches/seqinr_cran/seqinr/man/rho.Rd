\name{rho}
\alias{rho}
\title{Statistical composition in dinucleotides}
\description{
  Computes the rho statistic of all 16 dinucleotides on the given sequence
}
\usage{rho(sequence)}
\arguments{
  \item{sequence}{ A nucleic acids sequence }
}
\details{
  The rho statistic, as presented by Karlin S., Cardon LR.  
}
\value{
  a table containing the rho statistic for each dinucleotide
}
\references{
  Karlin S. and Cardon LR. (1994) Computational DNA sequence analysis.
\emph{Annu Rev Microbiol}, \bold{48}, 619--54.
  citation("seqinr")
}
\author{ Leonor Palmeira }
\seealso{\code{\link{zscore}}}
\examples{
sequence=sample(s2c('acgt'),6000,rep=TRUE)
rho(sequence)
}
\keyword{ utilities }
