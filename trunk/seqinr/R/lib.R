.First.lib <- function(lib, pkg)
{
  library.dynam("seqinr", pkg, lib)
  file1 <- file.path(lib,pkg,"data","SEQINR.UTIL.RData")
  load (file1,.GlobalEnv)
  load (file1,envir=as.environment(match("package:seqinr", search())))
  #if (!require(ade4))
     # stop("seqinr requires ade4, but ade4 couldn't be loaded")
}

.Last.lib <- function(libpath=.path.package(package = "seqinr"))
{
  if(exists(bankname)) rm(bankname)
}