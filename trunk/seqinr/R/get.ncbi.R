#
#                        get.ncbi
#
# Try to connect to ncbi to get a list of complete bacterial genomes.
# Returns an n by 2 dataframe with species names in first column and
# genbank accession numbers in second colum.
#

get.ncbi <- function()
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

  brut <- readLines("ftp://ftp.ncbi.nih.gov/genomes/Bacteria/gbk")

  #
  # Now there should be a vector of chr in "brut", one element per
  # line in the "gbk" file, each line looking like:
  #
  # "Aeropyrum_pernix   NC_000854.gbk"
  #

  brut <- sapply( brut, strsplit, split=" ")

  #
  # Now each element in "brut" should be splited as in:
  #
  # "Aeropyrum_pernix" ""                 ""                 "NC_000854.gbk"
  #

  brut <- sapply( brut, unlist)

  #
  # Now "brut" should be an 4 by n matrix as in:
  # > brut[ ,1]
  # [1] "Aeropyrum_pernix" ""                 ""                 "NC_000854.gbk"
  #

  brut <- brut[c(1,4),]

  #
  # Now "brut" should be an 2 by n matrix as in:
  # > brut[ ,1]
  # [1] "Aeropyrum_pernix" "NC_000854.gbk" 
  #

  brut <- t(brut)
  names(brut) <- NULL
  row.names(brut) <- NULL
  brut <- as.data.frame(brut)
  names(brut) <- c("species","gbk.acc.nbr")

  #
  # Now "brut" should be an n by 2 data.frame as in:
  # > brut[1, ] 
  #           species   gbk.acc.nbr
  # 1 Aeropyrum_pernix NC_000854.gbk
  #

# shouldn't ftp_proxy be restored there ?

  return(brut)
}

#
#             ncbi.fna.url 
#
#  Try to build urls to access complete genome sequences data
#  in fasta format from get.ncbi() output
#

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

#
#  Try to build urls to access complete genome sequences data
#  in genbank format from get.ncbi() output
#

ncbi.gbk.url <- function( get.ncbi.out = get.ncbi() )
{
  build.url <- function( x )
  {
    urlname = paste("ftp://ftp.ncbi.nih.gov/genomes/Bacteria/", x[1],
                 "/",x[2], collapse="", sep="")
  }
  apply( get.ncbi.out, 1, build.url )  
}
#
#  Try to build urls to access complete genome sequences data
#  file *.ptt from get.ncbi() output
#

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

#
# Try to get the number of cds and genome size
#

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
