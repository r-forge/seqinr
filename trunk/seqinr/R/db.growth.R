get.db.growth <- function()
{
  if (!capabilities("http/ftp")) 
    stop("capabilities(\"http/ftp\") is not TRUE")
  ftp.proxy.bck <- Sys.getenv("ftp_proxy")
  if (ftp.proxy.bck != "") {
      warning("I'am trying to neutralize proxies")
      Sys.putenv("no_proxy" = "")
  }

  embl <- "ftp://ftp.ebi.ac.uk/pub/databases/embl/release/relnotes.txt"
  tmp <- readLines( embl )
  idx <- grep("Release Month", tmp)
  tmp <- tmp[ (idx + 2):length(tmp) ]
  tmp <- strsplit( tmp, split = " " )
  not.empty <- function(x)
  {
    x <- x[nchar(x) > 0 ]
  }
  tmp <- sapply( tmp, not.empty )
  tmp <- data.frame( t(tmp) )
  names(tmp) <- c("Release", "Month", "Entries", "Nucleotides")

  tmp[,1] <- as.double( as.character(tmp[,1]))
  tmp[,3] <- as.double( as.character(tmp[,3]))
  tmp[,4] <- as.double( as.character(tmp[,4]))

  date  <- strsplit(tmp[,2], split="/")
  date.to.num <- function(x)
  {
    x <- as.double( x )
    return( (x[1]-1)/12 + x[2] )
  }
  date <- sapply(date, date.to.num)
  tmp <- data.frame( cbind(tmp, date) )
  return(tmp)
}

dia.db.growth <- function( get.db.growth.out = get.db.growth() )
{
  op <- par(no.readonly = TRUE)
  par( bg = "blue" )
  par( fg = "yellow" )
  par( col = "yellow" )
  par( col.axis = "yellow" )
  par( col.lab = "yellow" )
  par( col.main = "yellow" )
  par( col.sub = "yellow" )
  attach( get.db.growth.out )
  plot( date, log10(Nucleotides) , pch = 20,
    main = "The exponential growth of the DDBJ/EMBL/Genbank content",
    xlab = "Year", ylab = "Log10 number of nucleotides" )
  abline(lm(log10(Nucleotides)~date),col="yellow")
  lm1 <- lm(log(Nucleotides)~date)
  mu <- lm1$coef[2] # slope
  dbt <- log(2)/mu # doubling time
  dbt <- 12*dbt # in months
  text( x = 1995, y = 7, lab=paste("Doubling time:", round(dbt,1),"months"))
  par( op )
} 
