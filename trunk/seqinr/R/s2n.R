#
# simple numerical encoding of a DNA sequence that by default
# is independent of locale.
#
s2n <- function(seq, levels = c("a", "c", "g", "t"), base4 = TRUE, ... )
{
  if( base4 )
    codes(factor(seq, levels = levels , ...) ) - 1
  else
    codes(factor(seq, levels = levels , ...) )
}
