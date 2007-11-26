#
# To get annotations corresponding to sequences
#
getAnnot <- function(object, ...) UseMethod("getAnnot")

getAnnot.default <- function(object, ...)
  stop(paste("no getAnnot method for objects of class:", class(object)))

getAnnot.list <- function(object, ...)
  sapply(seq_len(length(object)), function(i) getAnnot(object[[i]], ...))

getAnnot.SeqFastadna <- function(object, ...) attr(object, "Annot")
getAnnot.SeqFastaAA <- getAnnot.SeqFastadna

getAnnot.SeqAcnucWeb <- function(object, ..., nbl = 100, socket = autosocket())		
  readAnnots.socket(socket, name = object, nl = nbl)

getAnnot.qaw <- function(object, ...) getAnnot(object$req, ...)

getAnnot.logical <- function (object, ...)
   object # so that NA is returned for virtual lists

