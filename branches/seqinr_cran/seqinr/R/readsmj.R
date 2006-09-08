###################################################################################################
#                                                                                                 #
#                              readsmj                                                            #
#                                                                                                 #
#                                                                                                 #
# ==>   readsmj&num=xx&nl=xx                                                                      #
# <==  code=xx&nl=xx                                                                              #
#     recnum=xx&name="xx"&plong=xx{&libel="xx"}                                                   #
#     ... a series of nl lines like that ...                                                      #
# Returns data from nl consecutive records starting from rank num of file SMJYT including label   #
# if not empty.                                                                                   #
# Code != 0 indicates error.                                                                      #
#                                                                                                 #
###################################################################################################

readsmj <- function(socket = "auto", num = 2, nl = 10, recnum.add = FALSE, nature.add = TRUE,
plong.add = FALSE, libel.add = FALSE, all.add = FALSE)
{ 
  #
  # Use default bank if no socket is given:
  #
  if (socket == "auto"){
    socket <- banknameSocket$socket
  }
  
  #
  # Turn all flags to TRUE when requested:
  #
  if(all.add) libel.add <- plong.add <- nature.add <- recnum.add <- TRUE

  
  #
  # Build the request:
  #
  request <- paste("readsmj&num=", num, "&nl=", nl, sep = "", collapse = "")
  
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
  # check answer from server:
  #
  if(rep[1] != "0") stop("Error from server")

  #
  # read answer from server:
  #
  n <- as.numeric(rep[2])
  ans <- readLines(socket, n = n)
  
  #
  # Put answer into a data.frame:
  #
  name <- character(n) # the only mandatory one
  if(recnum.add) recnum <- numeric(n)
  if(nature.add) nature <- factor(character(n), levels = c("00","01","02","03","04","05","06","07"))
  if(plong.add) plong <- numeric(n)
  if(libel.add) libel <- character(n)
  for(i in 1:n){
    tmp <- parser.socket(ans[[i]])
    name[i] <- substr(tmp[2], 2, nchar(tmp[2]) - 1)
    if(recnum.add) recnum[i] <- tmp[1]
    if(nature.add) nature[i] <- substr(name[i], 1, 2)
    if(plong.add) plong[i] <- tmp[3]
    if(libel.add){
      if(length(tmp) == 4){
        libel[i] <- substr(tmp[4], 2, nchar(tmp[4]) - 1)
      } else {
        libel[i] <- NA
      }
    }
  }
  df <- data.frame(list(name = I(name)))
  if(recnum.add) df$recnum <- recnum
  if(nature.add){
    levels(nature) <- c("status", "molecule", "journal", "year", "type", "organelle", "division", "dbstrucinfo")
    df$nature <- nature
  }
  if(plong.add) df$plong <- plong
  if(libel.add) df$libel <- libel
  
  return(df)
}
