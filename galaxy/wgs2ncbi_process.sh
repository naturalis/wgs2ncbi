#!/bin/bash
#                             $1          $2      $3
# usage: wgs2ncbi_process.sh '$argstring' '$gff3' '$out_file1'

tempdir=`mktemp -d 2>/dev/null || mktemp -d -t 'tempdir'`
datadir=$tempdir/datadir
gff3dir=$tempdir/gff3dir
mkdir $datadir
mkdir $gff3dir
unzip -qq -d $gff3dir $2
wgs2ncbi process $1 -datadir $datadir -gff3dir $gff3dir
zip -jrq $3.zip $datadir
mv $3.zip $3
rm -rf $tempdir
