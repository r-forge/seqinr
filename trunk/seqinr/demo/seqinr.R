par(ask = TRUE)


s = read.fasta("sequences/malM.fasta")

getFrag(s,1,50)
getName(s)
getAnnot(s)

tablecode(1)

Translate(s)

a = read.fasta("sequences/seqAA.fasta",seqtype="AA")

summary(a)

sock = choosebank.socket("genbank")

query.socket(sock$socket,"felis", "sp=felis catus et t=cds")
felis$req[1:10]
getAnnot(felis$req[[4]])
getFrag(felis$req[[300]],45,78)


length = lapply(felis$req[[1:10]],getLength)
seq = lapply(felis$req[[1:10]],getSequence)
res = lapply(seq,count,2)
plot(res)



#Load dataset:
data(ec999)
# Compute codon usage for all coding sequences:
ec999.uco <- lapply(ec999, uco) 
# Put it in a dataframe:
df <- as.data.frame(lapply(ec999.uco, as.vector))
# Add codon names:
row.names(df) <- names(ec999.uco[[1]])
# Compute global codon usage:
global <- rowSums(df)
# Choose a title for the graph:
title <- "Codon usage in 999 E. coli coding sequences"
# Plot data:
dotchart.uco(global, main = title) 


par(ask = FALSE)