#
# To get the location of subsequences from an ACNUC server
#

getLocation <- function(object, ...) UseMethod("getLocation")

getLocation.list <- function(object, ...)
  lapply(seq_len(length(object)), function(i) getLocation(object[[i]], ...))

getLocation.default <- function(object, ...)
  stop(paste("no getLocation method for objects of class:", class(object)))

getLocation.SeqAcnucWeb <- function(object, ..., socket = autosocket())
  unlist(getLocationSocket(socket, name = object))

getLocation.qaw <- function(object, ...) getLocation(object$req, ...)

getLocation.logical <- function (object, ...)
  object # so that NA is returned for virtual lists
