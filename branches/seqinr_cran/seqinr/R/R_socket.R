###################################################################################################
#                                                                                                 #
# Functions to communicate with a remote ACNUC database through sockets.                          #
#                                                                                                 #
###################################################################################################




###################################################################################################
#                                                                                                 #
#                                         parser.socket                                           #
#                                                                                                 #
#                      Utility function to parse answers from ACNUC server.                       #
#                                                                                                 #
###################################################################################################

parser.socket <- function(onelinefromserver)
{
  if(is.null(onelinefromserver)){
    return(NULL)
  }
  #
  # Answers from server looks like : "code=0&lrank=2&count=150513&type=SQ&locus=F"
  # 
  loc <- gregexpr("=[^=&]*", onelinefromserver)[[1]]
  substring(onelinefromserver, loc + 1, loc + attr(loc, "match.length") - 1)
}


################################################################################
#                                                                               
#                                         getSequenceSocket                    #
#                                                                               
################################################################################

getSequenceSocket <- function(socket, name, start, length, as.string = FALSE){
  request <- paste("gfrag&name=", name, "&start=", formatC(start, format = "d"),
                   "&length=", formatC(length, format = "d"), sep = "")
  writeLines(request, socket, sep = "\n")
  answerFromServer <- readLines(socket, n = 1)

  #
  # Check that there is an answer from server:
  #
  if(length(answerFromServer) == 0){
    warning(paste("Empty answer from server with sequence name:", name))
    return(NA)
  } else {
    #
    # Check that no error code is returned by server:
    #
    if(substr(x = answerFromServer, start = 1, stop = 5) == "code="){
      warning(paste("Server returned error code:", answerFromServer, "with sequence name:", name))
      return(NA)
    }
    #
    # Extract sequence from server answer:
    #
    sequence <- unlist(strsplit(answerFromServer, split = "&"))[2]
    #
    # Returns the sequence either as a string or as a vector of single chars:
    #
    if( as.string ){
      return(sequence)
    } else {
      return(s2c(sequence))
    }
  }
}
  
###################################################################################################
#                                                                                                 #
#                                         getAttributsocket                                       #
#                                                                                                 #
# To get sequence attributes from server.                                                         #
#                                                                                                 #
###################################################################################################

getAttributsocket <- function( socket, name){
  request <- paste("isenum&name=", name, sep = "")
  writeLines( request, socket, sep = "\n")
  res <- readLines(socket, n = 1)
  p <- parser.socket(res)
  return( list(length = as.numeric(p[2]), frame = as.numeric(p[3]), gencode = as.numeric(p[4])) )
}

###################################################################################################
#                                                                                                 #
#                                         readAnnots.socket                                       #
#                                                                                                 #
###################################################################################################

readAnnots.socket <- function(socket, name, nl){
#
# Check arguments:
#
  if(nl <= 0){
    warning("Negative or zero value for argument nl, forced to 1.")
    nl <- 1
  }
#
# Build request:
#
  request <- paste("read_annots&name=", name, "&nl=", nl, sep = "")
  writeLines(request , socket, sep="\n")
#
# Read result from server:
#  
  res <- readLines(socket , n = nl)
#
# Remove the "nl=xx&" answer from server on first line:
#
  newfirstline <- unlist(strsplit(res[1], "&"))[2]
  res[1] <- newfirstline
  return(res)
}

###################################################################################################
#                                                                                                 #
#                                         getNumber.socket                                        #
#                                                                                                 #
###################################################################################################

getNumber.socket <- function( socket, name){
  request <- paste("isenum&name=", name, sep = "")
  writeLines(request, socket, sep = "\n")
  s <- readLines(socket, n = 1)
  return(parser.socket(s)[1])
}



###################################################################################################
#                                                                                                 #
#                                         print.qaw                                               #
#                                                                                                 #
###################################################################################################

print.qaw <- function(x, ...)
{
  cat("\n")
  cat("\n$socket: ")
  print(x$socket)
  cat("\n$banque: ")
  #cat(get("bankName",env=.GlobalEnv)) # Ca pas bon
  cat("\n$call: ")
  print(x$call)
  cat("$name: ")
  print(x$name)
  cat("\n")
  sumry <- array("", c(1, 4), list(1, c("list", "length", "mode", "content")))
  sumry[1, ] <- c("$req",length(x$req),"character","sequences")
  class(sumry) <- "table"
  print(sumry)
  cat("\n")
}

###################################################################################################
#                                                                                                 #
#                                         getKeywordsocket                                        #
#                                                                                                 #
###################################################################################################

getKeywordsocket <- function(socket, name){
  #modif simon
  writeLines(paste("isenum&name=", name, sep = ""), socket, sep = "\n")
  res <- readLines(socket, n = 1)
  number <- parser.socket(res)[1] 

  writeLines(paste("readsub&num=", number, sep = ""), socket, sep = "\n")
  res2 <- readLines(socket, n = 1) 
  rr <- parser.socket(res2)

  writeLines(paste("readshrt&num=", rr[7], sep = ""), socket, sep = "\n")
  res3 <- readLines(socket, n = 1)
  #modif simon   

  # Get the nb of kw (not used here)
  # nbkws <- parser.socket(res3)[2]

  #recupere la liste de paires val, next 
  tmpl <- unlist(strsplit(res3, "&"))
  #transforme en liste
  tmpl <- unlist(strsplit(tmpl[3],","))
  kwl <- unlist(tmpl)[c(TRUE, FALSE)]

  lapply(kwl, function(x){
    writeLines(paste("readkey&num=", x, sep = ""), socket, sep = "\n")  
    res4 <- readLines(socket, n = 1)
    res <-parser.socket(res4)[2]
    substring(res[1], 2, nchar(res[1]) - 1)
  })

} 

###################################################################################################
#                                                                                                 #
#                                         getLocationSocket                                       #
#                                                                                                 #
###################################################################################################

getLocationSocket <- function( socket, name){

   writeLines(paste("isenum&name=",name,sep=""),socket,sep="\n")
         res = readLines( socket , n=1 )
         number = parser.socket(res)[1] 
   
         writeLines(paste("readsub&num=",number,sep=""),socket,sep="\n")
         res2 = readLines( socket , n=1 ) 
         rr = parser.socket(res2)
  
   # Test si subsequence           
     
   l=list() 
         if(as.numeric(rr[5]) != 0){
     warning("It's a parent sequence\n")
     return( NA )
    }
         else {
    i=1
    writeLines(paste("readext&num=",rr[6],sep=""),socket,sep="\n")    
    res3 = readLines( socket , n=1 )
    r = parser.socket(res3)
    l[[i]] = as.numeric(c(r[3],r[4]))
    n=r[5] 
  }
        while(as.numeric(n) != 0){
    i=i+1
    writeLines(paste("readext&num=",n,sep=""),socket,sep="\n")   
    res4 = readLines( socket , n=1 )
    rrr = parser.socket(res4)
    l[[i]] = as.numeric(c(rrr[3],rrr[4]))
    n=rrr[5]
      }
  return(l)
} 


###################################################################################################
#                                                                                                 #
#                              readfirstrec                                                       #
#                                                                                                 #
#                Returns the record count of the specified ACNUC index file.                      #                                                                                                 #
#                                                                                                 #
# ==>   readfirstrec&type=[AUT|BIB|ACC|SMJ|SUB|LOC|KEY|SPEC|SHRT|LNG|EXT|TXT]                     #
# <==  code=xx&count=xx                                                                           #
# Returns the record count of the specified ACNUC index file.                                     #
# Code != 0 indicates error.                                                                      #
#                                                                                                 #
###################################################################################################

readfirstrec <- function(socket = autosocket(), type)
{
  allowedtype <- c("AUT", "BIB", "ACC", "SMJ", "SUB", "LOC", "KEY", "SPEC", 
                   "SHRT", "LNG", "EXT", "TXT")
  if(missing(type)){
    return(allowedtype)
  }
  #
  # Build the request:
  #
  request <- paste("readfirstrec&type=", type, sep = "", collapse = "")
  
  #
  # Send request:
  #
  
  writeLines(request, socket, sep = "\n") 
  #
  # Read answer from server:
  #
  
  s <- readLines(socket, n = 1)
  rep <- parser.socket(s)
  
  #
  # Check answer from server:
  #
  if(rep[1] != "0"){
    warning("Server returns an error")
    return(NA)
  } else {
    return(as.numeric(rep[2]))
  }
}



