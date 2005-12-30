\name{aacost}
\alias{aacost}
\docType{data}
\title{Aerobic cost of amino-acids in Escherichia coli}
\description{
The metabolic cost of amino-acid biosynthesis in \emph{E. coli} under aerobic conditions
from table 1 in Akashi and Gojobori (2002).
}
\usage{data(aacost)}
\format{
  A data frame with 20 rows for the amino-acids and the following 6 columns:
  \describe{
    \item{aaa}{amino-acid (three-letters code).}
    \item{a}{amino-acid (one-letter code).}
    \item{prec}{precursor metabolites (see details).}
    \item{p}{number of high-energy phosphate bonds contained in ATP and GTP molecules.}
    \item{h}{number of available hydrogen atoms carried in NADH, NADPH,
    and FADH2 molcules.}
    \item{tot}{total metabolic cost assuming 2 high-energy phosphate bonds
    per hydrogen atom.}
  }
}
\details{
precursor metabolites: penP, ribose 5-phosphate; PRPP,
    5-phosphoribosyl pyrophosphate; eryP, erythrose 4-phosphate; 3pg,
    3-phosphoglycerate; pep, phospho\emph{enol}pyruvate; pyr, pyruvate;
    acCoA, acetyl-CoA; akg, alpha-ketoglutarate; oaa, oxaloacetate. Negative
    signs on precursor metabolites indicate chemicals \emph{gained} through
    biosynthetic pathways. Costs of precursors reflect averages for growth on
    glucose, acetate, and malate (see Table 6 in the supporting information
    from Akashi and Gojobori 2002).
}
\source{
Akashi, H, Gojobori, T. (2002) Metabolic efficiency and amino acid composition
in the proteomes of \emph{Escherichia coli} and \emph{Bacillus subtilis}.
\emph{Proceedings of the National Academy of Sciences of the United States of America}, 
\bold{99}:3695-3700.
}

\references{
\code{citation("seqinr")}.
}

\examples{
data(aacost)
}
\keyword{datasets}
