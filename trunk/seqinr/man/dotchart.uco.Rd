\name{dotchart.uco}
\alias{dotchart.uco}

\title{Cleveland plot for codon usage tables }
\description{
Draw a Cleveland dot plot for codon usage tables
}
\usage{
dotchart.uco(x, numcode = 1, aa3 = TRUE, cex = 0.7, alphabet = s2c("tcag"), pch = 21, gpch = 20, bg = par("bg"), color = par("fg"), gcolor = par("fg"), lcolor = "gray", xlim, main = NULL, xlab = NULL, ylab = NULL, ...)
}

\arguments{
  \item{x}{table of codon usage as computed by \code{uco}. }
  \item{numcode}{the number of the code to be used by \code{translate}.}
  \item{aa3}{logical. If TRUE use the three-letter code for amino-
acids. If FALSE use the one-letter code for amino-acids. }
  \item{cex}{the character size to be used. }
  \item{alphabet}{character for codons labels}
  \item{pch}{the plotting character or symbol to be used.}
  \item{gpch}{the plotting character or symbol to be used for group values. }
  \item{bg}{the background color to be used. }
  \item{color}{the color(s) to be used for points an labels. }
  \item{gcolor}{the single color to be used for group labels and values.}
  \item{lcolor}{the color(s) to be used for the horizontal lines.}
  \item{xlim}{horizontal range for the plot}
  \item{main}{overall title for the plot}
  \item{xlab}{X axis label}
  \item{ylab}{Y axis label}
  \item{\dots}{graphical parameters can also be specified as arguments}
}

\value{
An invisible list with components:
\item{x}{table of codon usage}
\item{labels}{codon names}
\item{groups}{amino acid factor}
\item{gdata}{sums by amino acid}
\item{ypg}{the y-axis coordinates for amino acids}
\item{ypi}{the y-axis coordinates for codons}
}

\references{ Cleveland, W. S. (1985) The Elements of Graphing Data. 
Monterey, CA: Wadsworth. }
\author{ J.R. Lobry }

\seealso{ \code{\link{dotchart}}, \code{\link{uco}}, \code{\link{aaa}},
\code{\link{translate}} }

\examples{
\dontrun{# Load dataset:}
\dontrun{data(ec999)}
\dontrun{# Compute codon usage for all coding sequences:}
\dontrun{ec999.uco <- lapply(ec999, uco) }
\dontrun{# Put it in a dataframe:}
\dontrun{df <- as.data.frame(lapply(ec999.uco, as.vector)) }
\dontrun{# Add codon names:}
\dontrun{row.names(df) <- names(ec999.uco[[1]])}
\dontrun{# Compute global codon usage:}
\dontrun{global <- rowSums(df)}
\dontrun{# Choose a title for the graph:}
\dontrun{title <- "Codon usage in 999 E. coli coding sequences"}
\dontrun{# Plot data:}
\dontrun{dotchart.uco(global, main = title) }

}
\keyword{hplot}
