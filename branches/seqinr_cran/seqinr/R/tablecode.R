#
# Genetic code table as in Text Books
#

tablecode <- function(numcode = 1, urn.rna = s2c("TCAG"), dia = FALSE,
latexfile = NULL, label = latexfile, size = "normalsize", caption = NULL)
{
  aa1 <- a()
  aa3 <- aaa()
  codename <- SEQINR.UTIL$CODES.NCBI[numcode, "ORGANISMES"]
  urn <- s2c("tcag") # internal
#
# Make default caption for LaTeX table:
#
  if( is.null(caption) ){
    caption <- paste("\\caption{Genetic code number ", 
                     numcode, ": ", codename, ".}", sep = "")
  }
#
# As a LaTex table:
#  
  if( ! is.null(latexfile) ) {
    Tfile <- file(latexfile, open = "w")
    writeLines("\\begin{table}", Tfile)
    writeLines("\\begin{center}", Tfile)
#
# Character size:
#
    writeLines(paste("{\\", size, sep = ""), Tfile)
    
    writeLines("\\begin{tabular}{*{13}{l}}", Tfile)
    writeLines("\\hline", Tfile)
    writeLines("\\\\", Tfile)
    
    ncodon <- 1 # codon rank as in words(alphabet = s2c("tcag"))
    for( i in 0:3 )
    {
      for( j in 0:3 )
      {
        for( k in 0:3 )
        {
          codon <- c(urn[i+1], urn[k+1], urn[j+1])
          codon.urn <- paste(urn.rna[i+1], urn.rna[k+1], urn.rna[j+1], sep = "", collapse = "")
          aminoacid <- aa3[which(aa1 == translate(codon, numcode = numcode))]

          writeLines(paste(codon.urn, aminoacid, " &", sep = " & "), Tfile)
          ncodon <- ncodon + 1
        }
        writeLines("\\\\", Tfile)
      }
      writeLines("\\\\", Tfile)
    }
    writeLines("\\hline", Tfile)
    writeLines("\\end{tabular}", Tfile)
#
# Caption:
#
    writeLines(caption, Tfile)
#
# LaTeX label:
#
    writeLines(paste("\\label{", label, "}", sep = ""), Tfile)
#
# End character size:
#
    writeLines("}", Tfile)

    writeLines("\\end{center}", Tfile)
    writeLines("\\end{table}", Tfile)
    close(Tfile)
    return(invisible(NULL))
  }
#
# END LATEX
#
  if( dia )
  {  
    op <- par(no.readonly = TRUE)
    par(bg = "blue")
    par(fg = "yellow")
    par(col = "yellow")
    par(col.axis = "yellow")
    par(col.lab = "yellow")
    par(col.main = "yellow")
    par(col.sub = "yellow")
  }

  plot.new()
  plot.window(xlim=c(0,100),ylim=c(0,100))

  segments( 0, 102, 100, 102, lwd = 2)
  segments( 0, 0, 100, 0, lwd = 2)
  segments( 0, 97, 100, 97)


  text(x=0, y = 98.5, font = 2, adj = c(0, 0),
    lab = paste("Genetic code", numcode,":",codename))

  urn <- c("t","c","a","g") # internal
  for( i in 0:3 )
  {
    for( j in 0:3 )
    {
      for( k in 0:3 )
      {
        codon <- c(urn[i+1], urn[j+1], urn[k+1])

        text( x = 100*j/4, y = 95 - 100*i/4 -5*k, adj = c(-0.5,1.5),
        lab = urn.rna[i+1] )

        text( x = 100*j/4 + 3, y = 95 - 100*i/4 -5*k, adj = c(-0.5,1.5),
        lab = urn.rna[j+1] )

        text( x = 100*j/4 + 6, y = 95 - 100*i/4 -5*k, adj = c(-0.5,1.5),
        lab = urn.rna[k+1] )

        aminoacid <- aa3[which(aa1 == translate(codon, numcode = numcode))]
        text( x = 100*j/4 + 12, y = 95 - 100*i/4 -5*k, adj = c(-0.5,1.5),
        lab =  aminoacid )
      }
    }
  }
  if(dia)
    par(op)
}
