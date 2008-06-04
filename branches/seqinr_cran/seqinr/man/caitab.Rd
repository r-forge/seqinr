\name{caitab}
\alias{caitab}
\docType{data}
\title{Codon Adaptation Index (CAI) tables}
\description{
Information about a preferred set of codons for highly expressed genes
in three species.
}
\usage{data(caitab)}
\format{
  A data frame with 64 rows for the codons and the following 3 columns:
  \describe{
    \item{ec}{\emph{Escherichia coli}}
    \item{bs}{\emph{Bacillus subtilis}}
    \item{sc}{\emph{Saccharomyces cerevisiae}}  
  }
}
\details{
}
\source{
The C program codonW version 1.4.4 writen by John Peden
available at \url{http://codonw.sourceforge.net/}. The data
are from the file codonW.h, see CAI STRUCT. 
}
\references{
For the \emph{Saccharomyces cerevisiae} data: 
Sharp and Cowe (1991) Yeast 7:657-678.

\code{citation("seqinr")}
}
\seealso{
  \code{\link{cai}} for an example using this dataset to compute CAI values.
}
\examples{
data(caitab)
}
\keyword{datasets}
