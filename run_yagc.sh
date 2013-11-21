#!/bin/bash

# verbosity, this writes extra logging to STDERR, increment as '-v' or '-v -v'
VERBOSE=

# gff3 directory, must NOT exist if the intention is to re-explode the GFF3
GFF3DIR=share/gff3

# directory for tbl and fasta files
TBLDIR=share/tblfasta-new
if [ ! -d $TBLDIR ]; then
	mkdir -p $TBLDIR
fi

# validation directory, contains val reports, ASN.1 output (and genbank files)
ASN1DIR=share/asn1val-new
if [ ! -d $ASN1DIR ]; then
	mkdir -p $ASN1DIR
fi

# the source of the annotations we include in the output, optional, this
# corresponds with column 2 in the GFF3
SOURCE=maker

# prefix below which to generate locus tags, should be registered as bioproject
PREFIX=L345_

# input genome
FASTA=share/king_cobra_scaffolds_spring_2011-cleaned.fasta

# input annotation file
GFF3=share/cobra.functional-cleaned.gff

# number of contigs to chunk together in one file (<=10k)
CHUNKSIZE=5000

# extra metadata to include in the fasta header lines
INFO=share/info.ini

# AGP settings
CENTER=NaturalisBC
EVIDENCE=paired-ends
GAP_TYPE=scaffold
LINKAGE='yes'
ORGANISM='Ophiophagus hannah'
TAXID=8665
NAME=PRJNA201683
AGP=outfile.agp

# create the AGP, this step is not necessary in our case, according to NCBI
# perl agpmaker.pl -i $FASTA -e $EVIDENCE -g $GAP_TYPE -o $ORGANISM -t $TAXID -n $NAME \
# 	-c $CENTER > $AGP

# namespace below which to scope protein and transcript IDs
AUTHORITY="gnl|${CENTER}|"

# limit number of written scaffolds, for testing
LIMIT=0

# minimum length of scaffolds to retain
MINLENGTH=200

# minimum length of introns to believe
MININTRON=10

# discrepancy report file, compare with https://www.ncbi.nlm.nih.gov/genbank/asndisc
DISCREP=share/discrep-new.txt

# the submission template as produced by http://www.ncbi.nlm.nih.gov/WebSub/template.cgi
TEMPLATE=share/template.sbt

# explode the annotations into files 
if [ ! -d $GFF3DIR ]; then	
	mkdir -p $GFF3DIR
	perl -Ilib script/explode_gff3.pl -gff3 $GFF3 -source $SOURCE -f CDS -f gene \
	-f five_prime_UTR -f three_prime_UTR -d $GFF3DIR
fi

# create the table and fasta files
perl -Ilib -MBio::WGS2NCBI -e run - -d $TBLDIR -s $SOURCE -p $PREFIX -f $FASTA \
	-g $GFF3DIR -info $INFO -a $AUTHORITY -c $CHUNKSIZE -l $LIMIT -minlength $MINLENGTH \
 	-minintron $MININTRON $VERBOSE 

# create the genbank and validation files
tbl2asn -p $TBLDIR -t $TEMPLATE -M n -a r10k -l paired-ends -r $ASN1DIR -Z $DISCREP -V b