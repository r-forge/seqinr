kaks <- function(x, debug = FALSE){
    #
    # Check argument class:
    #
    if(attr(x,"class") != "alignment") error("object x must be of class alignment")
    if(debug){
      cat("<--- Argument x storage is --->\n")
      print(str(x))
      cat("<--- Argument x storage is --->\n")
    }
    #
    # Check that there are at least two sequences in the alignment:
    #
    if(x$nb <= 1) return(0)
    #
    # Call internal C function:
    #
    l <- .Call("kaks", x$seq, x$nb, debug, PACKAGE = "seqinr")
    if(debug){
      cat("<--- Result l storage is --->\n")
      print(str(l))
      cat("<--- Result l storage is --->\n")
    }
    #
    # This is to compute the list of results:
    #
    mkresult <- function(k){
      if(! is.null(x$nam)){
        tmp <- matrix( k, x$nb, x$nb, byrow = TRUE, dimnames = list(x$nam,x$nam))
      } else {
        n <- paste("seq", c(1:x$nb), sep = "")
        tmp <- matrix( k, x$nb, x$nb, byrow = TRUE, dimnames = list(n,n))
      }
      as.dist(t(tmp))
    }
    m <- lapply(l[1:4], mkresult)
    names(m) <- c("ka", "ks", "vka", "vks")
    return(m)
}

