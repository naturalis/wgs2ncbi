#!/bin/bash
# usage: wgs2ncbi_trim.sh $DATA $OUTFILE

# make temporary directory, assign location to variable $tmpdir
tmpdir=`mktemp -d 2>/dev/null || mktemp -d -t 'tmpdir'`

# extract FASTA files and feature tables from input archive, write to $tmpdir
unzip $1 *.fsa *.tbl -d $tmpdir

# run the trip operation
wgs2ncbi trim -datadir $tmpdir

# remove backup files from $tmpdir
rm $tmpdir/*.bak

# compress $tmpdir to provided output name (needs .zip extension)
zip -jrq $2.zip $tmpdir

# rename to final output name
mv $2.zip $2

# remove temporary directory
rm -rf $tmpdir
