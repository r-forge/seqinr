#
# To get Keywords associated with a sequence.
#

getKeyword <- function(object, ...) UseMethod("getKeyword")

getKeyword.default <- function(object, ...)
  stop(paste("no getKeyword method for objects of class:", class(object)))

getKeyword.list <- function(object, ...)
  lapply(seq_len(length(object)), function(i) getKeyword(object[[i]], ...))

getKeyword.SeqAcnucWeb <- function(object, ..., socket = autosocket())
  unlist(getKeywordsocket(socket, name = object))

getKeyword.qaw <- function(object, ...) getKeyword(object$req, ...)

getKeyword.logical <- function (object, ...)
  object # so that NA is returned for virtual lists
