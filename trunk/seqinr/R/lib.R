.First.lib <- function(lib, pkg)
{
  library.dynam("seqinr", pkg, lib)
  file1 <- file.path(lib,pkg,"data","SEQINR.UTIL.RData")
  load (file1,.GlobalEnv)
  load (file1,envir=as.environment(match("package:seqinr", search())))
  if(exists("bankname",envir=globalenv())) rm(bankname,envir=globalenv())
  #if (!require(ade4))
     # stop("seqinr requires ade4, but ade4 couldn't be loaded")
}

.Last.lib <- function(libpath=.path.package(package = "seqinr"))
{
  if(exists("SEQINR.UTIL",envir=globalenv())) rm(SEQINR.UTIL,envir=globalenv())
}
