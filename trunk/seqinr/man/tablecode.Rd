\name{tablecode}
\alias{tablecode}

\title{ to plot genetic code as in textbooks }
\description{
This function plots a genetic code table as in textbooks, that is
following the order U > C > A > G so that synonymous codons
are almost always in the same boxes.
}
\usage{
tablecode(numcode = 1, urn.rna = c("U", "C", "A", "G"), dia = TRUE)
}

\arguments{
  \item{numcode}{ The genetic code number as in \code{translate} }
  \item{urn.rna}{ The letters to display codons, use c("T","C",
"A","G") if you want the code in terms of DNA sequence}
  \item{dia}{ to produce a yellow/blue plot for slides }
}
\author{J.R. Lobry }

\seealso{ \code{\link{translate}} }
\examples{
 tablecode()
}
\keyword{ utilities }
