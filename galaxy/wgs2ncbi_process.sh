#!/bin/bash
#                             $1          $2      $3
# usage: wgs2ncbi_process.sh '$argstring' '$gff3' '$out_file1'

tempdir=`mktemp -d 2>/dev/null || mktemp -d -t 'tempdir'`
datadir=$tempdir/datadir
gff3dir=$tempdir/gff3dir
mkdir $datadir
mkdir $gff3dir
tar -xzf $2 -C $gff3dir --strip-components=1
wgs2ncbi process $1 -datadir $datadir -gff3dir $gff3dir
tar -C $datadir -cvzf $3 . > /dev/null 2>&1
rm -rf $tempdir > /dev/null 2>&1
