#
# simple numerical encoding of a DNA sequence that by default
# is independent of locale.
#
s2n <- function(seq, levels = c("a", "c", "g", "t"), ... )
{
  codes(factor(seq, levels = levels , ...) )
}
