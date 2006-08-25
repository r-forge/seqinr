\name{tablecode}
\alias{tablecode}
\title{ to plot genetic code as in textbooks }
\description{
This function plots a genetic code table as in textbooks, that is
following the order \code{T > C > A > G} so that synonymous codons
are almost always in the same boxes.
}
\usage{
tablecode(numcode = 1, urn.rna = s2c("TCAG"), dia = FALSE, latexfile = NULL,
label = latexfile, size = "normalsize", caption = NULL)
}
\arguments{
  \item{numcode}{The genetic code number as in \code{translate} }
  \item{urn.rna}{The letters to display codons, use s2c("UCAG")
  if you want the code in terms of RNA sequence}
  \item{latexfile}{The name of a LaTex file if you want to redirect the output}
  \item{label}{The label for the LaTeX table}
  \item{size}{The LaTex size of characters for the LaTeX table}
  \item{caption}{The caption of the LaTeX table}
  \item{dia}{to produce a yellow/blue plot for slides }
}
\references{
  \code{citation("seqinr")}
}
\author{ J.R. Lobry }
\seealso{ \code{\link{translate}}, \code{\link{syncodons}} }
\examples{
#
# Show me the standard genetic code:
#

 tablecode()
}
\keyword{ utilities }
