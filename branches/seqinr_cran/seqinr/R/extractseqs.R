# Recuperation des sequences
# Simon Octobre 2006




###################################################################################################
#                                                                                                 #
#                                         extractseqs.socket                                       #
#                                                                                                 #
###################################################################################################

extractseqs <- function( lrankseqnum,socket = "auto", format="fasta",operation="simple",zlib=FALSE, verbose = FALSE, npaquets = -1, nzlines=1000){

debug<-0
if (verbose)
	debug <- 1

  if(verbose) cat("I'm checking the arguments...\n")

  if (socket == "auto"){
    if(verbose) cat("No socket were specified, using default.\n")
    socket <- banknameSocket$socket
    }
    
  if( !inherits(socket, "sockconn") ) stop(paste("argument socket = ", socket, "is not a socket connection."))
  if( !is.character(lrankseqnum) ) stop(paste("argument lrankseqnum = ", lrankseqnum, "is not a character string."))
  if(verbose) cat("... and everything is OK up to now.\n")
  
    
# Check arguments:
# Check  if  lrankseqnum is a list or a sequence 
# Check if  format is acnuc", "fasta", or "flat"
# Check if operation "is simple", "translate", "fragment", "feature" or "region"
#
# Build request:
#
# if  lrankseqnum is a list

  lrank <-glr(lrankseqnum)
  if(verbose) cat("The rank of the list ",lrankseqnum, "is ",lrank,".\n")
  if (is.na(lrank)) {
  	stop(paste("Problem in rank list!\n"))
  	}
  
  if (zlib == FALSE) {
  	if( !inherits(socket, "sockconn") ) stop(paste("argument socket = ", socket, "is not a socket connection."))
  	if( !is.character(lrankseqnum) ) stop(paste("argument lrankseqnum = ", lrankseqnum, "is not a character string."))
  	if(verbose) cat("... and everything is OK up to now.\n")
  	request <- paste("extractseqs&lrank=", lrank, "&format=", format, "&operation=", operation,"&zlib=F", sep = "")
  	writeLines(request , socket, sep="\n") 
   
	# Read result from server: 
  	res <- readLines(socket , n = npaquets)
  
	# Check if server answer is correct 
  	if(verbose) cat("I'm trying to analyse answer from server...\n")
  	p <- parser.socket(res[1])
  	if (p != 0) stop(paste("Problem in server answer!\n"))
  	if(verbose) cat("... and everything is OK up to now.\n")

	# Get results   
  	testend <-res[length(res)];
  	while (testend != "extractseqs END.") {
        	if(verbose) cat("extract new packet....\n")
  		newres <-  readLines(socket , n = npaquets)
  		res=c(res,newres)
  		testend <-res[length(res)];
  		}
  
  	ii=1;
  	jj=0;
  	lastres=c()
  	while (ii < length(res)) {
  		toto<-unlist(strsplit(res[ii], "count"))[1]
  		if (toto !="\033"){
			lastres[jj]=res[ii]
  			jj = jj +1
			}
  		ii=ii+1
  		}
  
  
  	} #  if (zlib == FALSE)
  else {
  
   	request <- paste("extractseqs&lrank=", lrank, "&format=", format, "&operation=", operation,"&zlib=T", sep = "")
   	if(verbose) cat("request : ",request,"\n")
   	writeLines(request , socket, sep="\n")
	
	# Read result from server: 
   	lastres <-  .Call("getzlibsock", socket, nzlines, debug,PACKAGE = "seqinr")
  	}
  return(lastres);
}
exseq <- extractseqs


zlibexseq <- function(soc,nlz,ver)  .Call("getzlibsock", soc, nlz, ver,PACKAGE = "seqinr")
