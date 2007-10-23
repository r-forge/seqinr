###########################################################################
# 
#                                choosebank
#
# To select an ACNUC database or to get the list of available databases
# from an ACNUC server.
# 
###########################################################################

choosebank <- function(bank = NA,
                       host = "pbil.univ-lyon1.fr",
                       port = 5558,
                       verbose = FALSE,
                       timeout = 5,
                       infobank = FALSE,
                       tagbank = NA){
  #
  # Print parameter values if verbose mode is on:
  #
  if(verbose){ 
    cat("Verbose mode is on, parameter values are:\n")
    cat(paste("  bank = ", deparse(substitute(bank)), "\n"))
    cat(paste("  host = ", deparse(substitute(host)), "\n"))
    cat(paste("  port = ", deparse(substitute(port)), "\n"))
    cat(paste("  timeout = ", deparse(substitute(timeout)), "seconds \n"))
    cat(paste("  infobank = ", deparse(substitute(infobank)), "\n"))
    cat(paste("  tagbank = ", deparse(substitute(tagbank)), "\n"))
  }
  
  #
  # Check parameter values (to be completed):
  #
  if( !is.na(tagbank) ){
    if(verbose) cat("I'm checking the tagbank parameter value...\n")
    if( !(tagbank %in% c("TEST", "TP", "DEV")) ){
      if(verbose) cat("... and I was able to detect an error.\n")
      stop("non allowed value for tagbank parameter.\n")
    } else {
      if(verbose) cat("... and everything is OK up to now.\n")
    }  
  }
  
  #
  # Check that sockets are available:
  #
  if(verbose) cat("I'm ckecking that sockets are available on this build of R...\n")
  if( !capabilities("sockets") ){
    stop("Sockets are not available on this build of R.")
   } else {
    if(verbose) cat("... yes, sockets are available on this build of R.\n")
  }
  
  # 
  # Try to open socket connection:
  #
  if(verbose) cat("I'm trying to open the socket connection...\n")
  oldtimeout <- getOption("timeout")
  options(timeout = timeout)
  socket <- try( socketConnection( host = host, port = port, server = FALSE, blocking = TRUE))
  options(timeout = oldtimeout)
  if(inherits(socket, "try-error")) {
    errmess <- paste("I wasn't able to open the socket connection:\n",
                     "  o Check that your are connected to the internet.\n",
                     "  o Check that port", port, "is not closed by a firewall.\n",
                     "  o Try to increase timeout value (current is", timeout, "seconds).\n")
    stop(errmess)
  } else {
    if(verbose) cat("... yes, I was able to open the socket connection.\n")
  }

  #
  # Read the answer from server:
  #
  if(verbose) cat("I'm trying to read answer from server...\n")
  rep1 <- readLines(socket, n = 1)
  if(verbose) cat(paste("... answer from server is:", rep1, "\n"))
  
  #
  # Client ID definition : seqinr + package version number
  #  (internal note: log file is: /mnt/users/ADE-User/racnuc/log)
  #
  clientID <- paste("seqinr_", packageDescription("seqinr")$Version, sep = "")
  if(verbose) cat(paste("I'm trying to identify myself as", clientID, "to the server...\n"))
  request <- paste("clientid&id=", clientID, sep = "")
  writeLines( request, socket, sep = "\n")
  rep <- readLines(socket, n = 1)  
  if(verbose) cat(paste("... answer from server is:", rep, "\n"))
  res <- parser.socket(rep)  
  if( res[1] == "0") {
    if(verbose) cat("... and everything is OK up to now.\n")  
  } else {
    stop("I don't know what this error code means for clientid, please contact package maintener.\n")
  }
           
  ###############################################################################
  #
  # If no bank name is given, return the list of available banks from server:
  #
  ###############################################################################
  resdf <- kdb(tag = tagbank, socket = socket)
  nbank <- nrow(resdf)
  
  if( is.na(bank) ){  
    if(verbose) cat("No bank argument was given...\n")
    if( !infobank ){
      if(verbose) cat("infobank parameter is FALSE, I'm just returning bank names\n")
      return(resdf$bank)
    } else {
      if(verbose) cat("infobank parameter is TRUE, I'm returning all bank infos\n")
      return(resdf) 
      }
  } else {

    ###############################################################################
    #
    # If a bank name is given, try to open it from server:
    #
    ###############################################################################
    
    # 
    # Try to open bank from server:
    #
    if(verbose) cat("I'm trying to open the bank from server...\n")
    resacnucopen <- acnucopen(bank, socket)
    if(verbose) cat("... and everything is OK up to now.\n")

    #
    # Try to get informations from HELP file: 
    #
    if(verbose) cat("I'm trying to get information on the bank...\n")
    bankhelp <- ghelp(item = "CONT", file = "HELP", socket = socket, catresult = FALSE)
    bankrel <- bankhelp[1]
    
    #
    # Try to get status info:
    #
    status <- "unknown"
    for(i in seq_len(nbank)){
     if (resdf[i,1] == bank) status <- resdf[i,2]
    }
  
    #
    # Set up the ACNUC server for following queries:
    #
    if(verbose) cat("I'm trying to set up the server for following queries...\n")
    writeLines("prep_requete", socket, sep = "\n")
    rep3 <- readLines(socket, n = 1)
    
    #
    # Re-patch pas beau:
    #
    if(length(rep3) == 0){
   if(verbose) cat("... answer from server is empty!\n")
   while(length(rep3) == 0){
     if(verbose) cat("... reading again.\n")
     rep3 <- readLines(socket, n = 1)
   }
    }
    if(verbose) cat("... answer from server is: ", rep3, "\n")
    res3 <- parser.socket(rep3)
    if( res3[1] == "0") {
   if(verbose) {
     cat("... and everything is OK up to now.\n")
     cat(paste("... and there are", res3[2], "free lists available from server.\n"))
   }
    } else {
   if(verbose) cat("I was able to detect an error while seting up remote bank.\n")
   stop("There was an error while seting up remote bank.\n")
    }
    
    #
    # Build result and assign it in the global environment:
    #
    res <- list(socket = socket,
     bankname = bank,
     banktype = resacnucopen$type,
     totseqs = resacnucopen$totseqs,
     totspecs = resacnucopen$totspecs,
     totkeys = resacnucopen$totkeys,
     release = bankrel,
     status = status,
     details = bankhelp)
    assign("banknameSocket", res, .GlobalEnv)
    invisible(res)
  }
} 
