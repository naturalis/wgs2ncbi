#!/bin/bash

# IN SOME CASES, YOUR SUBMISSION MAY NEED TO HAVE AN AGP FILE. TO THIS END IS PROVIDED
# A SEPARATE SCRIPT THAT GENERATES SUCH A FILE. BELOW ARE THE VARIABLES THAT THIS SCRIPT
# WILL NEED TO KNOW ABOUT:

# CENTER=NaturalisBC
# EVIDENCE=paired-ends
# GAP_TYPE=scaffold
# LINKAGE='yes'
# ORGANISM='Ophiophagus hannah'
# TAXID=8665
# NAME=PRJNA201683
# AGP=outfile.agp

# THIS WILL CREATE THE AGP FILE:

# perl agpmaker.pl -i $FASTA -e $EVIDENCE -g $GAP_TYPE -o $ORGANISM -t $TAXID -n $NAME -c $CENTER > $AGP

# explode the annotations into files 
wgs2ncbi prepare -conf share share/wgs2ncbi.ini

# create the table and fasta files
wgs2ncbi process -conf share/wgs2ncbi.ini

# run tbl2asn
wgs2ncbi convert -conf share/wgs2ncbi.ini

# compress for upload
wgs2ncbi compress -conf share/wgs2ncbi.ini
