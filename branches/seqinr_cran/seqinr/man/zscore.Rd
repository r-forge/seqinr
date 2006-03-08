\name{zscore}
\alias{zscore}
\title{Statistical over- and under- representation in dinucleotides}
\description{
  Computes the z-score of dinucleotides according to a specific model of
  random sequence generation
}
\usage{zscore(sequence, simulations = NULL, modele, ... )}
\arguments{
  \item{sequence}{ A nucleic acids sequence }
  \item{simulations}{ If \code{NULL}, analytical solution is computed
    when available. Otherwise, it should be the number of permutations
    for the z-score computation }
  \item{modele}{ A string of characters describing the model chosen for
    the random generation }
}
\details{
  
}
\value{
  a table containing the z-score statistic for each dinucleotide
}
\references{
  citation("seqinr")
}
\author{ Leonor Palmeira }
\seealso{ \code{\link{permutation}} }
\examples{
sequence=sample(s2c('acgt'),6000,rep=TRUE)
zscore(sequence,modele='base')
zscore(sequence,modele='codon')
zscore(sequence,1000,modele='syncodon')
}
\keyword{ utilities }
