#!/bin/csh

setenv ACNUC_CSRC "/bge/banques/csrc"

# 
# Remove all old symbolic links:
# 

ls -F | grep "@" | tr -d "@" > todelete
awk ' { printf("\\rm %s\n", $0) } ' todelete > todelete.csh
chmod +x todelete.csh
./todelete.csh

\rm todelete
\rm todelete.csh

#
# Add symbolic link to relevant ACNUC sources:
#

ln -fs $ACNUC_CSRC/initf.c .
ln -fs $ACNUC_CSRC/gbemgener.c .
