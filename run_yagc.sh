#!/bin/bash

# AGP settings
# CENTER=NaturalisBC
# EVIDENCE=paired-ends
# GAP_TYPE=scaffold
# LINKAGE='yes'
# ORGANISM='Ophiophagus hannah'
# TAXID=8665
# NAME=PRJNA201683
# AGP=outfile.agp

# create the AGP, this step is not necessary in our case, according to NCBI
# perl agpmaker.pl -i $FASTA -e $EVIDENCE -g $GAP_TYPE -o $ORGANISM -t $TAXID -n $NAME \
# 	-c $CENTER > $AGP

# explode the annotations into files 
# if [ ! -d $GFF3DIR ]; then	
# 	mkdir -p $GFF3DIR
# 	perl -Ilib script/explode_gff3.pl -gff3 $GFF3 -source $SOURCE -f CDS -f gene \
# 	-f five_prime_UTR -f three_prime_UTR -d $GFF3DIR
# fi

# create the table and fasta files
# perl -Ilib -MBio::WGS2NCBI -e process - -conf share/wgs2ncbi.ini

perl -Ilib -MBio::WGS2NCBI -e convert - -conf share/wgs2ncbi.ini

# the interface should be:
# wgs2ncbi prepare -conf share/wgs2ncbi.ini
# wgs2ncbi process -conf share/wgs2ncbi.ini
# wgs2ncbi convert -conf share/wgs2ncbi.ini
# wgs2ncbi compress -conf share/wgs2ncbi.ini

# create the genbank and validation files
# tbl2asn -p $TBLDIR -t $TEMPLATE -M n -a r10k -l paired-ends -r $ASN1DIR -Z $DISCREP -V b