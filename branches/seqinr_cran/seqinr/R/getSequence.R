#
# To get sequence data
#

getSequence <- function(object, as.string = FALSE, ...) UseMethod("getSequence")

getSequence.default <- function(object, as.string = FALSE, ...)
    stop(paste("no getSequence method for objects of class:", class(object)))

getSequence.list <- function(object, as.string = FALSE, ...)
  sapply(seq_len(length(object)), function(i) getSequence(object[[i]], as.string = as.string, ...))

getSequence.character <- function(object, as.string = FALSE, ...){
  is.single.string <- function(x) length(x) == 1 && nchar(x) > 1
  if(is.single.string(object)){
    if(as.string) return(as.list(object)) else return(s2c(object))
  } else {
    if(as.string) return(as.list(c2s(object))) else return(object)
  }
}

getSequence.SeqFastadna <- function(object, as.string = FALSE, ...){
  attributes(object) <- NULL # not needed here
  getSequence.character(object, as.string, ...)
}
getSequence.SeqFrag <- getSequence.SeqFastaAA <- getSequence.SeqFastadna

getSequence.SeqAcnucWeb <- function(object, as.string = FALSE, ..., socket = autosocket())
  getSequenceSocket(socket, object, start = 1, length = attr(object, "length"), as.string = as.string)

getSequence.qaw <- function(object, as.string = FALSE, ...) getSequence(object$req, ...)

getSequence.logical <- function (object, as.string = FALSE, ...)
  object # so that NA is returned for virtual lists
