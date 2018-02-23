#!/bin/bash
#                             $1       $2      $3      $4       $5          $6      $7                             
# usage: wgs2ncbi_process.sh '$fasta' '$conf' '$info' '$masks' '$products' '$gff3' '$out_file1'

tempdir=`mktemp -d 2>/dev/null || mktemp -d -t 'tempdir'`
datadir=$tempdir/datadir
gff3dir=$tempdir/gff3dir
mkdir $datadir
mkdir $gff3dir
tar -xzf $6 -C $gff3dir --strip-components=1
wgs2ncbi process -conf $2 -datadir $datadir -datafile $1 -info $3 -masks $4 -products $5 -gff3dir $gff3dir -verbose
tar -C $datadir -cvzf $7 . 2>/dev/null
rm -rf $tempdir