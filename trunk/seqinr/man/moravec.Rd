\name{moravec}
\alias{moravec}
\non_function{}
\title{ Evolution of Computer Power/Cost }
\usage{data(moravec)}
\description{
The \code{moravec} data frame has 134 rows and 7 columns.
Year, Performance and cost for 134 computers.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Machine}{a character vector for the computer name}
    \item{Year}{a numeric vector}
    \item{Cost.Kdol}{a numeric vector for the cost in Kdol}
    \item{Kdol.1997}{a numeric vector for the cost in 1997 Kdol}
    \item{MBytes}{a numeric vector for RAM capacity}
    \item{MIPS}{a numeric vector for CPU performance}
    \item{MIPS.per.Kdol.1997}{a numeric vector for the Power/Cost ratio}
  }
}

\source{
Data of appenxix A from Moravec, H. (1998) When will computer harware
match the humain brain? Journal of Evolution and Technology, Vol 1.
\url{http://www.jetpress.org/volume1/moravec.htm}
}
\examples{
data(moravec)
plot(moravec$Year, log(moravec$MIPS.per))
}
\keyword{datasets}
