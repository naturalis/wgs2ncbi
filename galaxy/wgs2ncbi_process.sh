#!/bin/bash
# usage: wgs2ncbi_process.sh $FASTA $CONF_INI $OUTFILE

tempdir=`mktemp -d 2>/dev/null || mktemp -d -t 'tempdir'`
datadir=$tempdir/datadir
gff3dir=$tempdir/gff3dir
mkdir $datadir
mkdir $gff3dir
tar -xzf $3 -C $gff3dir --strip-components=1
wgs2ncbi process -conf $2 -datadir $datadir -datafile $1 -gff3dir $gff3dir -verbose
tar -C $datadir -cvzf $4 . 2>/dev/null
rm -rf $tempdir