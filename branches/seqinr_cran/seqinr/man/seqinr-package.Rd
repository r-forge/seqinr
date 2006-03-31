\name{seqinr-package}
\alias{seqinr-package}
\alias{seqinr}
\docType{package}
\title{
Biological Sequences Retrieval and Analysis
}
\description{
Exploratory data analysis and data visualization for
biological sequence (DNA and protein) data. Include also
utilities for sequence data manadgment under the ACNUC system.
}
\details{
\tabular{ll}{
Package: \tab seqinr\cr
Version: \tab 1.0-4\cr
Date: \tab 2006-03-31\cr
License: \tab GPL version 2 or newer\cr
URL: \tab http://pbil.univ-lyon1.fr/software/SeqinR/seqinr_home.php, Mailing
list:
http://pbil.univ-lyon1.fr/software/SeqinR/seqinr3_angl.php\cr
ZipData: \tab no\cr
Packaged: \tab Tue Mar 28 16:28:59 2006; simon\cr
Built: \tab R 2.2.1; i486-pc-linux-gnu; 2006-03-31 13:45:40; unix\cr
}

Index:
\preformatted{
AAstat                  To Get Some Protein Statistics
AnoukResult             Expected numeric results for Ka and Ks
                        computation
EXP                     Vectors of coefficients to compute linear
                        forms.
GC                      Calculates the fractional G+C content of
                        nucleic acid sequences
SEQINR.UTIL             utility data for seqinr
SeqAcnucWeb             Sequence coming from an ACNUC data bases
                        located on the web
SeqFastaAA              AA sequence in Fasta Format
SeqFastadna             Class for DNA sequence in Fasta Format
SeqFrag                 Class for sub-sequences
a                       Converts amino-acid three-letter code into the
                        one-letter one
aaa                     Converts amino-acid one-letter code into the
                        three-letter one
aacost                  Aerobic cost of amino-acids in Escherichia coli
                        and G+C classes
c2s                     conversion of a vector of chars into a string
chargaff                Base composition in ssDNA for 7 bacterial DNA
choosebank              To select a database structured under ACNUC and
                        located on the web
closebank               To close a remote ACNUC database
comp                    complement a nucleic sequence
computePI               To Compute the Theoretical Isoelectric Point
count                   Composition of dimer/trimer/etc nucleotides
dia.bactgensize         Distribution of bacterial genome size from GOLD
dinucl                  Mean zscore on 242 complete bacterial
                        chromosomes
dist.alignment          Pairwise Distances from Aligned Protein or
                        DNA/RNA Sequences
dotchart.uco            Cleveland plot for codon usage tables
ec999                   999 coding sequences from E. coli
gb2fasta                conversion of GenBank file into fasta file
gbk2g2                  Conversion of a GenBank format file into a
                        glimmer-like one
get.db.growth           Get the exponential growth of nucleic acid
                        database content
get.ncbi                Bacterial complete genome data from ncbi ftp
                        site
getAnnot                Generic Function to obtain annotation for a
                        sequence
getFrag                 Generic Function to obtain a fragment of a
                        sequence
getKeyword              Generic Function to get keyword(s) associated
                        to a sequence
getLength               Generic Function To obtain the length of a
                        sequence
getLocation             Generic Function To obtain the Position of a
                        Sequence on the Parent Sequence
getName                 Generic Function To obtain the name of a
                        sequence
getSequence             Generic Function To Get a sequence
getTrans                Generic Function to translate a sequence
getType                 To Get Sequence Type
invers                  invers
kaks                    to Get an Estimation of Ka and Ks
lseqinr                 To see what's inside the package seqinr
n2s                     function to convert the numeric encoding of a
                        DNA sequence into a vector of characters
oriloc                  Prediction of origin and terminus of
                        replication in bacteria
permutation             Sequence permutation according to several
                        different models
plot.SeqAcnucWeb        To Plot Subsequences on the Parent Sequence
prochlo                 Zscore on three strains of Prochlorococcus
                        marinus
query                   To get a list of sequence names from an ACNUC
                        data base located on the web
read.alignment          Read aligned sequence files in mase, clustal,
                        phylip, fasta or msf format
read.fasta              Read file in fasta format
rho                     Statistical over- and under- representation of
                        dinucleotides in a sequence
s2c                     conversion of a string into a vector of chars
s2n                     simple numerical encoding of a DNA sequence.
splitseq                split a sequence into sub-sequences
syncodons               Synonymous codons
synsequence             Random synonymous coding sequence generation
tablecode               to plot genetic code as in textbooks
toyaa                   A toy example of amino-acid counts in three
                        proteins
toycodon                A toy example of codon counts in three coding
                        sequences
translate               Translate nucleic acid sequences
uco                     Codon usage indices
ucoweight               Weight of each synonymous codon
words                   To get all words from an alphabet.
words.pos               Positions of possibly degenerated motifs within
                        sequences
}
}
\author{
Delphine Charif <charif@biomserv.univ-lyon1.fr> and Jean Lobry
<lobry@biomserv.univ-lyon1.fr> and Leonor Palmeira
<palmeira@biomserv.univ-lyon1.fr>

Maintainer: Simon Penel <penel@biomserv.univ-lyon1.fr>
}
\keyword{ package }
