#
# Conversion of the numeric encoding of a DNA sequence into
# a vector of chars
#
n2s <- function(nseq, levels = c("a", "c", "g", "t"))
{
  levels[nseq]
}
