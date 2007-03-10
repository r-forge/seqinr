\name{extract.breakpoints}
\alias{extract.breakpoints}
\title{Extraction of breakpoint positions on the
  rearranged nucleotide skews. }
\description{
  Extraction of breakpoint positions on the
  rearranged nucleotide skews. 
  }

\usage{extract.breakpoints(rearr.ori,type=c("atfw","atrev","gcfw","gcrev"),nbreaks,gridsize=100,it.max=500)}

\arguments{
\item{rearr.ori}{ }
\item{type}{}
\item{nbreaks}{}
\item{gridsize}{}
\item{it.max}{}
}
\details{

}

\value{

}

\references{
  Necsulea, A. and Lobry, J.R. (in prep) A novel method for assessing
  the effect of replication on DNA base composition asymmetry.
}
\author{A. Necsulea}
\seealso{ \code{\link{oriloc}}, \code{\link{draw.rearranged.oriloc}},
  \code{\link{rearranged.oriloc}}  }

\examples{


### Example for Chlamydia trachomatis ####

### Rearrange the chromosome and compute the nucleotide skews ###

r.ori=rearranged.oriloc(seq.fasta = system.file("sequences/ct.fasta",package = "seqinr"),
    g2.coord = system.file("sequences/ct.coord",package = "seqinr"))

### Extract the breakpoints for the rearranged nucleotide skews ###

breaks=extract.breakpoints(r.ori,type=c("atfw","atrev","gcfw","gcrev"),nbreaks=c(2,2,2,2),gridsize=100,it.max=500)

summary(breaks)

### Draw the rearranged nucleotide skews and place the position of the
breakpoints on the graphics ###

draw.rearranged.oriloc(r.ori,breaks.gcfw=breaks$gcfw$breaks,breaks.gcrev=breaks$gcrev$breaks,breaks.atfw=breaks$atfw$breaks,breaks.atrev=breaks$atrev$breaks)

}
\keyword{utilities}
