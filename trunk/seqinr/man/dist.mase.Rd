\name{dist.mase}
\alias{dist.mase}
\title{ Pairwise Distances from Protein Sequences }
\description{
  These function compute a matrix of pairwise distances from aligned protein sequences using similarity (Fitch matrix) or identity matrix.
}
\usage{
dist.mase(mase, matrix = "similarity")
}
\arguments{
  \item{mase}{ object of class \code{SeqMase} }
  \item{matrix}{ the euclidean matrix distance to be used }
}
\value{
  The distance matrix, obect of class \code{dist}, computed by using the specified distance measure.
}
\references{ 
Fitch, W.M. (1966) Mutation values for the interconversion of amino acid pair. 
\emph{J. Mol. Biol.}, \bold{16}:9-16.\cr
}
\author{D. Charif}
\seealso{ \code{\link{read.mase} }}
\examples{
 mase=read.mase()
 dist.mase(mase)
}
\keyword{ manip }
