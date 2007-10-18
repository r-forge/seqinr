################################################################################
#
#                              plot.SeqAcnucWeb
#
################################################################################

plot.SeqAcnucWeb <- function(x, types = getType()$sname, socket = "auto", ...){
  verbose <- FALSE # a passer en argument si besoin est
  #
  # Check arguments:
  #
  if(!inherits(x, "SeqAcnucWeb")) stop("Sequence of class SeqAcnucWeb is needed")

  #
  # Use default bank if no socket is given:
  #
  if(socket == "auto"){
    socket <- get("banknameSocket", .GlobalEnv)$socket
  }
  #
  # Save graphical parameters:
  #
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))

  if(verbose) cat(paste("types:", types, sep = "\n"))
  

  #
  # Get the parent sequence:
  #
  
  GiveMeTheMotherSequence <- paste("me n=", x, sep = "") 
  query(listname = "me", query = GiveMeTheMotherSequence, socket = socket)
  MotherLength <- as.numeric(getLength(get("me", .GlobalEnv)$req[[1]]))
  MotherName <- get("me", .GlobalEnv)$req[[1]]
  if(verbose) cat("\nMotherLength = ", MotherLength)

  #
  # Plot organization:
  #
  par(mar = c(2.1, 0.1, 4.1, 0.1), lend = "square", ljoin = "mitre")
  cx <- c(0, MotherLength)
  cy <- c(0, 1)
  plot(cx, cy, ann = FALSE, type = "n", axes = FALSE)
  axis(1)
  title(main = paste("Physical position of subsequences on the parent sequence",
          MotherName, "(", MotherLength, "bp )", sep=" "))

  # Si x est une sous séquence alors dire le type

  if(get("me", .GlobalEnv)$req[[1]] != x){
    if(verbose) cat("x est une sous-sequence\n")
    writeLines(paste("isenum&name=",x,sep=""),socket,sep="\n")
    res <- readLines( socket , n=1 )
    writeLines(paste("readsub&num=",as.numeric(parser.socket(res)[1]),sep=""),socket,sep="\n")
    res <- readLines( socket , n=1 )
    writeLines(paste("readsmj&num=",as.numeric(parser.socket(res)[4]),"&nl=1",sep=""),socket,sep="\n")
    res <- readLines( socket , n=2 )
    ty <- substring(noquote(parser.socket(res[2]))[2],4,nchar(noquote(parser.socket(res[2]))[2])-1)
    p <- getLocationSocket(socket,x)
    lapply(p,function(x){rect(x[1],0,x[2],1,col="red", border="red", lend = "butt" )})
    plot(c(0, MotherLength),c(0,10),type="n",axes=FALSE,ann=FALSE)
    title("Legend",font.main=4)
    legend(9,legend=ty,fill="red",bg="cornsilk",ncol = 1)
    names(p) <- x
    return(p)
    
  } else{
    if(verbose) cat("x n'est pas une sous-sequence\n")
    q <- paste("fi n=", x, sep = "")
    query(listname = "filles", query = q, socket = socket)
    
    if( get("filles", .GlobalEnv)$nelem == 1 && get("filles", .GlobalEnv)$req == x){
      if(verbose) cat("No subsequences case\n")
      rect(0, 1, getLength(x), 2, col = "red", border = "red", lend = "butt")
      legend(9, legend = x, fill = "red", bg = "cornsilk", ncol = 1)
      invisible(list(x)) 
    }

    n <- length(types) # number of potential subsequences types
    ispresent <- rep(FALSE, n) # will be TRUE if subsequence if found
    nb <- numeric(n) # count of subsequences for available types
    posi <- vector(mode = "list", length = n) # position of subsequences
    
    for(i in seq_len(n)){
      q <- paste("filles et t=", types[i], sep = "")
      if(verbose) cat("query = ", q, "\n")

      result <- try(query(socket = socket, listname = "tmp", query = q))
      if( inherits(result, "try-error")) next
      if(get("tmp", .GlobalEnv)$nelem == 0) next
      if(is.na(get("tmp", .GlobalEnv)$req[[1]])) next
      if(get("tmp", .GlobalEnv)$req[[1]] == x ) next

      ispresent[i] <- TRUE
      
      u <- lapply(get("tmp", .GlobalEnv)$req, getLocation)
      names(u) <- get("tmp", .GlobalEnv)$req
      nb[i] <- length(u)
      posi[[i]] <- u
    }

    #
    # Draw subsequences:
    #
    posi <- posi[ispresent]
    nb <- nb[ispresent]
    types <- types[ispresent]
    n <- length(types)
    for(i in seq_len(n)){
      for(j in seq_len(length(posi[[i]]))){
        xleft <- posi[[i]][[j]][1]
        ybottom <- (i - 1)/(n + 1)
        xright <- posi[[i]][[j]][2]
        ytop <- i/(n + 1)
        rect(xleft, ybottom, xright, ytop, col = i, border = "black", lend = "square", ljoin = "mitre" )
      }
    }
    
    #
    # Draw legend:
    #
    legend("top", legend = paste(types, "(", nb, ")", sep = ""), fill = seq_len(n), horiz = TRUE, bty = "n")

    resu <- lapply(posi,function(x){lapply(x,unlist)})
    names(resu) <- types
  }
  #
  #  workspace cleanup
  #
  rm("me", pos = .GlobalEnv)
  rm("filles", pos = .GlobalEnv)
  rm("tmp", pos = .GlobalEnv)
  #
  # Return invisibly the result:
  #
  invisible(resu)
}
