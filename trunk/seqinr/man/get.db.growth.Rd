\name{get.db.growth}
\alias{get.db.growth}
\alias{dia.db.growth}
\title{ Get the exponential growth of nucleic acid database content }
\description{
  Connects to the embl database to read the last release note about
the number of nucleotides in the DDBJ/EMBL/Genbank database content.
A log-linear fit is represented by dia.bd.gowth() with an estimate of
the doubling time in months.
}
\usage{
get.db.growth()
dia.db.growth( get.db.growth.out = get.db.growth() )
}
\arguments{
  \item{get.db.growth.out}{ the output from get.db.growth() }
}
\value{
  A dataframe with the statistics from the embl site.
}
\references{\url{ftp://ftp.ebi.ac.uk/pub/databases/embl/release/relnotes.txt} }
\author{ J.R. Lobry }

\examples{
\dontrun{data <- get.db.growth()}
\dontrun{dia.db.growth(data)}
}
\keyword{utilities}
