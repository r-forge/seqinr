#
# Constant for the location of the complete bacterial genome sequences
# at ncbi ftp repository site
#

.BACTNCBI <- "ftp://ftp.ncbi.nih.gov/genomes/Bacteria/"

########################################################################
#                        get.ncbi
#
# Try to connect to ncbi to get a list of complete bacterial genomes.
# Returns an n by 4 dataframe with 
#   species names
#   accesion number
#   size in bp
#   type (chromosome or plasmid)
#   Last update time
#
#########################################################################
get.ncbi <- function(repository = .BACTNCBI)
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

  #
  # Try to get list of folder in ncbi repository. Note that R build-in ftp
  # does not allow to do this directly, so we rely on a system call that
  # will run only under Unix systems
  #
  cmd <- sprintf("echo \"ls\" | ftp %s", repository)
  brut <- readLines(pipe(cmd))
  
  
  #
  # Keep only lines corresponding to folders:
  #
  brut <- brut[grep("dr-xr-xr-x", brut)]


  #
  # Now there should be a vector of chr in "brut", each line looking like:
  #
  # "dr-xr-xr-x   2 ftp      anonymous     4096 Jul 15 15:41 Aeropyrum_pernix"
  #

  brut <- sapply( brut, strsplit, split=" ")
  names(brut) <- NULL

  #
  # Now each element in "brut" should be splited as in:
  #
  # [1] "dr-xr-xr-x"       ""                 ""                 "2"               
  # [5] "ftp"              ""                 ""                 ""                
  # [9] ""                 ""                 "anonymous"        ""                
  # [13] ""                 ""                 ""                 "4096"            
  # [17] "Jul"              "15"               "15:41"            "Aeropyrum_pernix"
  #

  get.last <- function( vector )
  {
    return( vector[length(vector)] )
  }
  brut <- sapply( brut, get.last)

  #
  # Now "brut" should contains folders names as in:
  # > brut[1:5]
  # [1] "Aeropyrum_pernix"                     
  # [2] "Agrobacterium_tumefaciens_C58_Cereon" 
  # [3] "Agrobacterium_tumefaciens_C58_UWash"  
  # [4] "Aquifex_aeolicus"                     
  # [5] "Archaeoglobus_fulgidus"               
  #

  #
  # Set vector types for results:
  #
  
  species <- character(0)
  accession <- character(0)
  size.bp <- integer(0)
  type <- character(0)
  lastupdate <- character(0)
  
  #
  # Main loop on folders to see what's inside
  #
  for( folder in brut )
  {
    where <- paste(repository, folder, "/", sep="", collapse="")
    cmd <- sprintf("echo \"ls\" | ftp %s", where)
    whatsin <- readLines(pipe(cmd))
    closeAllConnections()
    whatsin <- whatsin[ grep("\\.gbk", whatsin)] # Keep only files with ".gbk" extension
    
    for( i in seq(from=1, to=length(whatsin), by=1 )) # Loop on sequences data
    {
      #
      # Try to get the accession number of this entry:
      #
      accname <- unlist(strsplit(whatsin[i], split=" "))
      accname <- accname[length(accname)]
      accname <- unlist(strsplit(accname, split="\\."))
      accname <- accname[1] # The accession number should be in this variable
      #
      # Try to get the last update date of this entry
      #
      last <- unlist(strsplit(whatsin[i], split=" "))
      last <- last[nchar(last) > 0]
      last <- last[(length(last) - 3):(length(last) - 1)]
      last <- paste(last, collapse = " ")
      #
      # Try to get the size of this entry:
      #
      entry <- paste(repository, folder, "/", accname, ".gbk", sep="", collapse="")
      header <- readLines(entry, n=2)
      closeAllConnections()
      bp <- unlist(strsplit(header[1], split=" "))
      bp <- bp[nchar(bp) > 0]
      bp <- bp[3] # size in bp should be there
      species <- c(species, folder)
      accession <- c(accession, accname)
      size.bp <- c(size.bp, as.integer(bp))
      lastupdate <- c(lastupdate, last)
      #
      # Try to get the type (chromosome versus plasmid) of this entry:
      #
      if( length(grep("plasmid", tolower(header[2]))) != 0 )
        def <- "plasmid"
      else if(length(grep("chromosome", tolower(header[2]))) != 0)
        def <- "chromosome"
      else if(length(grep("genome", tolower(header[2]))) != 0)
        def <- "chromosome"
      else
        def <- NA
      type <- c(type, def)
      cat("\n",folder,accname,bp,def,last,"\n")
    }
  }

# shouldn't ftp_proxy be restored there ?
  return(data.frame(I(species), I(accession), size.bp, type, I(lastupdate)))
}




########################################################################
#             ncbi.fna.url 
#
#  Try to build urls to access complete genome sequences data
#  in fasta format from get.ncbi() output
#
########################################################################

ncbi.fna.url <- function( get.ncbi.out = get.ncbi() )
{
  build.url <- function( x )
  {
    ficname <- unlist(strsplit(x[2],"\\.")) # split prefix and suffix
    ficname <- ficname[1] # keep prefix
    ficname <- paste( ficname, ".fna", collapse="", sep="")
    urlname = paste("ftp://ftp.ncbi.nih.gov/genomes/Bacteria/", x[1],
                 "/",ficname, collapse="", sep="")
  }
  apply( get.ncbi.out, 1, build.url )  
}

########################################################################
#
#  Try to build urls to access complete genome sequences data
#  in genbank format from get.ncbi() output
#
########################################################################

ncbi.gbk.url <- function( get.ncbi.out = get.ncbi() )
{
  build.url <- function( x )
  {
    urlname = paste("ftp://ftp.ncbi.nih.gov/genomes/Bacteria/", x[1],
                 "/",x[2], collapse="", sep="")
  }
  apply( get.ncbi.out, 1, build.url )  
}
########################################################################
#  Try to build urls to access complete genome sequences data
#  file *.ptt from get.ncbi() output
#
########################################################################

ncbi.ptt.url <- function( get.ncbi.out = get.ncbi() )
{
  build.url <- function( x )
  {
    ficname <- unlist(strsplit(x[2],"\\.")) # split prefix and suffix
    ficname <- ficname[1] # keep prefix
    ficname <- paste( ficname, ".ptt", collapse="", sep="")
    urlname = paste("ftp://ftp.ncbi.nih.gov/genomes/Bacteria/", x[1],
                 "/",ficname, collapse="", sep="")
  }
  apply( get.ncbi.out, 1, build.url )  
}

########################################################################
#
# Try to get the number of cds and genome size
#
########################################################################

ncbi.stats <- function( get.ncbi.out = get.ncbi() )
{
  gbkurls <- ncbi.gbk.url( get.ncbi.out )
  ptturls <- ncbi.ptt.url( get.ncbi.out )
  get.genome.size <- function( url )
  {
    header <- readLines( url, n = 1 )
    tmp <- unlist(strsplit(header, split=" "))
    tmp <- tmp[nchar(tmp)>0]
    tmp <- tmp[3]
    as.integer(tmp)
  }
  get.n.prot <- function( url )
  {
    lines <- readLines( url, n = 3 )
    lines <- lines[3]
    lines <- unlist(strsplit(lines, split=" "))
    print(lines[1])
    as.integer(lines[1])
  }

  gsize <- sapply( gbkurls, get.genome.size)
  nprot <- sapply( ptturls, get.n.prot )
  data <- data.frame( cbind(get.ncbi.out[,1], gsize, nprot) )
  return( data )
}
