\name{kaks}
\alias{kaks}
\title{ to have an estimation of ka and ks }
\description{
  KS and KA distance are respectively the number of substitutions per synonymous site and per nonsynonymous site between two protein-coding genes. This function compute .
}
\usage{
kaks(x)
}
\arguments{
  \item{x}{ an object of class \code{alignment} }
}
\details{
 
}
\value{
  \item{ ks }{ matrix of ks }
  \item{ ka }{ matrix of ka }
  \item{ vks }{ variance matrix of ks }
  \item{ vka }{ variance matrix of ka }	
}
\references{ Li WH. (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. 
\emph{J Mol Evol.}, \bold{Jan};36(1):96-9.\cr 

    Li WH. Unbiased estimation of the rates of synonymous and nonsynonymous substitution. }
\author{ M. Gouy et D. Charif }
\note{}
\seealso{\code{\link{read.alignment}} }
\examples{
 s = read.alignment(File=system.file("sequences/test.phylip", package = "seqinr"),format="phylip")
 kaks(s)	
}
\keyword{ manip }

