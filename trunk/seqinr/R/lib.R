.First.lib <- function(lib, pkg)
{
  library.dynam("seqinr", pkg, lib)
   file1 <- file.path(lib,pkg,"data",".SEQINR.UTIL.RData")
   load (file1,.GlobalEnv)
   require("methods")		
}

