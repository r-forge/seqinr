#!/bin/csh

# 
# Remove all old symbolic links:
# 

awk ' { printf("\\rm -f %s\n", $0) } ' ".#acnucfile" > todelete.csh
chmod +x todelete.csh
todelete.csh

\rm todelete.csh

#
# Add symbolic link to relevant ACNUC sources:
#

awk ' { printf("ln -fs /bge/banques/csrc/%s  .\n", $0) } ' ".#acnucfile" > tolink.csh
chmod +x tolink.csh
./tolink.csh

\rm tolink.csh
