\name{get.ncbi}
\alias{get.ncbi}
\alias{ncbi.fna.url}

\title{ Bacterial complete genome data from ncbi ftp site }
\description{
Try to connect to ncbi ftp site to get a list of complete bacterial genomes.
}
\usage{
get.ncbi()
ncbi.fna.url(get.ncbi.out = get.ncbi() )
}
\arguments{
  \item{get.ncbi.out}{ the output from get.ncbi() }
}
\value{
Returns an n by 2 dataframe with species names in first column and
genbank accession numbers in second colum. Use ncbi.fna.url to
convert this into urls pointing to the corresponding fasta file.
}

\author{J.R. Lobry }

\section{WARNING }{
This function is highly dependant on ncbi ftp site conventions
for which we have no control.
The ftp connection apparently does not work when there is a proxy,
this problem is circumvented here in a rather crude way.
}

\examples{
\dontrun{bacteria <- get.ncbi()}
\dontrun{bacteria}
\dontrun{urls <- ncbi.fna.url( bacteria )}
\dontrun{urls}
\dontrun{seq <- read.fasta( urls[1] )}
}
\keyword{utilities}
