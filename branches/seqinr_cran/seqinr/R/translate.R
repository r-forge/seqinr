translate <- function(seq, frame = 0, sens = "F", numcode = 1, NAstring = "X")
{
  #
  # Take the reverse complementary strand when required:
  #

  if(sens == "R") seq <- comp(rev(seq))

  #
  # Transform the sequence in its numerical encoding equivalent
  # with textbook order, that is t = 0, c = 1, a = 2, g = 3
  #

  seq <- s2n(seq, levels = s2c("tcag"))

  #
  # Compute the length of the sequence when its length in codons
  # is an integer:
  #
  
  l <- 3*((length(seq) - frame) %/% 3)

  #
  # Compute the indices for the first codon positions:
  #
  
  c1 <- seq(from = frame + 1, to = frame + l, by = 3)
  
  #
  # Compute the indices of codons in the translation table:
  #

  tra <-  16*seq[c1] + 4*seq[c1 + 1] + seq[c1 + 2] + 1

  #
  # Get the translation table:
  #

  code <- s2c(SEQINR.UTIL$CODES.NCBI$CODES[numcode])

  #
  # Translate the sequence:
  #

  result <- code[tra]

  #
  # Replace missing values by the string for missing amino-acids:
  #

  result[is.na(result)] <- NAstring

  return( result )
}

