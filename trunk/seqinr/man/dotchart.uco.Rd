\name{dotchart.uco}
\alias{dotchart.uco}

\title{Draw a Cleveland dot plot for codon usage tables }

\usage{
dotchart.uco(x, numcode = 1, aa3 = TRUE, cex = 0.7, alphabet = s2c("TCAG"), pch = 21, gpch = 20, bg = par("bg"), color = par("fg"), gcolor = par("fg"), lcolor = "gray", xlim = c(0, max(x)), main = NULL, xlab = NULL, ylab = NULL, ...)
}

\arguments{
  \item{x}{table of codon usage as computed by \code{uco}. }
  \item{numcode}{the number of the code to be used bt \code{translate}.}
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

\references{ Cleveland, W. S. (1985) The Elements of Graphing Data. 
Monterey, CA: Wadsworth. }
\author{ J.R. Lobry }

\seealso{ \code{\link{dotchart}}, \code{\link{uco}}, \code{\link{aaa}},
\code{\link{translate}} }
\examples{
\dontrun{data(ec999)    # load dataset }
\dontrun{ec999.uco <- lapply(ec999,uco) # compute codon usage for all CDS}
\dontrun{df <- as.data.frame(lapply(ec999.uco,as.vector)) # put it in a 
dataframe}
\dontrun{row.names(df) <- names(ec999.uco[[1]]) # add codon names}
\dontrun{global <- rowSums(df)}
\dontrun{dotchart.uco(global, main="Codon usage in 999 E. coli coding 
sequences") }

}
\keyword{ graphics }
