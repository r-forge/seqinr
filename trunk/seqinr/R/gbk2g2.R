########################################################################
#
#                                 gbk2g2
#
# Conversion of a genbank format into a run-glimmer2 like format
#
#
########################################################################

gbk2g2 <- function(
  gbkfile = "ftp://pbil.univ-lyon1.fr/pub/logiciel/oriloc/ct.gbk",
  g2.coord = "g2.coord")
{
  #
  # First of all, check that this computer is not off the net:
  #
  if( ! capabilities("http/ftp") )
    stop("capabilities(\"http/ftp\") is not TRUE") 

#  
# BEGIN Proxy problem
#
# I have a problem here: the ftp connection apparently does
# not work when there is a proxy. I have fixed the bug
# this way, but this is a rather crude and unsatisfactory
# solution.
#
  ftp.proxy.bck <- Sys.getenv("ftp_proxy")

  if( ftp.proxy.bck != "" ) # there is a proxy
  {
    warning("I'am trying to neutralize proxies")
    Sys.putenv("no_proxy" = "") 
  }
#
# END Proxy problem
#  

  input <- readLines(gbkfile)
  
  outfile = file( description = g2.coord, open ="w")
  #
  # Keep lines with CDS flag:
  #
  input <- input[ substring(input,1,8) == "     CDS"]
  
  #
  # Extract boudaries strings
  #
  get.boundaries <- function( line )
  {
    line <- unlist(strsplit(line, split = " "))
    line <- line[nchar(line) > 1]
    return(line[2])
  }
  input <- sapply(input, get.boundaries)
  names(input) <- NULL
  
  #
  # Look for 5' partial genes:
  #
  idx5p <- grep(">", input)
  if( length( idx5p != 0 ) )
    warning("5' partial genes encountered (no output):", idx5p)
    
  #
  # Look for 3' partial genes:
  #
  idx3p <- grep("<", input)
  if( length( idx5p != 0 ) )
    warning("3' partial genes encountered (no output):", idx3p)
    
  #
  # Look for join in features:
  #
  idxjoin <- grep("join", input)
  if( length( idx5p != 0 ) )
    warning("join encountered (no output):", idxjoin)  
    
  #
  # Define partials and join:
  #
  censored <-  c(idx5p, idx3p, idxjoin)
  
  #
  # Extract boundaries:
  #
  for( i in seq(from = 1, to = length(input), by = 1 ) )
  {
    if( i %in% censored )
      next
    
    tmp <- unlist( strsplit(input[i], split="\\.\\.") )
    if( length( grep("complement", input[i]) ) == 1 )
    {
      end <- as.integer( substring(tmp[1], first = 12 ) ) + as.integer(3)
      start <- as.integer( substring(tmp[2], first = 1, last = nchar(tmp[2]) - 1 ))
      line <- sprintf(fmt = "%d\t%d\t%d", as.integer(i), start, end) 
    }
    else #direct strand
    {
      start <- as.integer(tmp[1])
      end   <- as.integer(tmp[2]) + as.integer(-3)
      line <- sprintf(fmt = "%d\t%d\t%d", as.integer(i), start, end) 
    }
    writeLines( line, con = outfile )
  }
  close(outfile)
  invisible(input)  
}

