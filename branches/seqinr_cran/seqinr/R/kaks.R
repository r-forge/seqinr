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
    # Make a local copy of the alignment (why? should be explained):
    #
    tmp <- x
    
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
    # l[[5]] contains the sequences, so this is to copy the sequences into
    # the seq component of tmp. But the sequences were already there since
    # tmp is a local copy of x. I don't understand this, the alignment is
    # not supposed to be modified by kaks !!! This should be explained.
    #
    tmp$seq <- l[[5]]
    #
    #
    #
    if(debug) print(as.character(as.list(match.call())$x))
    #
    # Oh! oh! This is to modify the objet x in the user environment! Why
    # why why ?? Extremly dangerous and undocumented. This may explain
    # why two succesive calls to kaks do not yield the same result.
    #
    assign(as.character(as.list(match.call())$x), tmp, envir = globalenv())
    if(debug){
      cat("<--- tmp storage is --->\n")
      print(str(tmp))
      cat("<--- tmp storage is --->\n")
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
    
    #
    # There is no valid reason to round the results, this should be deleted:
    #
    # m = lapply(m,round,digits=6)
    #
    names(m) <- c("ka", "ks", "vka", "vks")
    return(m)
}

