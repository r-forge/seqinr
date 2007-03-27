\name{draw.oriloc}
\alias{draw.oriloc}
\title{Graphical representation for nucleotide skews in
  prokaryotic chromosomes.}
\description{
  Graphical representation for nucleotide skews in
  prokaryotic chromosomes.}

\usage{draw.rearranged.oriloc(ori)}

\arguments{

\item{ori}{A data frame obtained with the \code{oriloc}
  function. }



}

\details{
  
}

\references{
 
}
\author{Jean R. Lobry }
\seealso{ \code{\link{rearranged.oriloc}},
  \code{\link{extract.breakpoints}}  }

\examples{

### Example for Chlamydia trachomatis ####


\dontrun{ori <- oriloc(seq.fasta = system.file("sequences/ct.fasta",package = "seqinr"),
    g2.coord = system.file("sequences/ct.coord",package = "seqinr"))}


\dontrun{draw.oriloc(ori) }





\keyword{utilities}