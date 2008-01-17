count <- function(seq, word, start = 0, freq = FALSE, alphabet = s2c("acgt"), frame = start){
  if(!missing(frame)) warning("argument frame is deprecated, use start instead")
#
# l is the number of elements in the sequence when starting at "start" position
# and ending so that the last word is documented:
#
  l <- length(seq) - word - start + 1
#
# s contains all words starting characters:
#
  s <- seq[(start + 1):(start + l)]
#
# After this statment seq contains all characters useful to build words:
#
  seq <- seq[(start + 1):length(seq)]
  if (word == 1){
    counts <- table(factor(seq, levels = levels(as.factor(words(1, alphabet = alphabet)))))
  } else {
    for(i in 2:word){
      #
      # Here we build all words by pasting the following characters:
      #
      s <- paste(s, seq[i:(i + l - 1)], sep = "")
    }
    counts <- table(factor(s, levels = levels(as.factor(words(word, alphabet = alphabet)))))
  }
  
  if (freq == FALSE){
    return(counts)
  }
  else{
    return((counts/(length(seq) - word + 1)))
  }
}
