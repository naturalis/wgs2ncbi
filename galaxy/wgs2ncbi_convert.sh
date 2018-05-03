#!/bin/bash
#                             $1      $2          $3
# usage: wgs2ncbi_convert.sh '$args' '$out_file1' '$data'

tempdir=`mktemp -d 2>/dev/null || mktemp -d -t 'tempdir'`
outdir=$tempdir/outdir
datadir=$tempdir/datadir
mkdir $outdir
mkdir $datadir
#tar -xzf $3 -C $datadir --strip-components=1
unzip -qq -d $datadir $3
wgs2ncbi convert $1 -datadir $datadir -outdir $outdir
#tar -C $outdir -cvzf $2 . > /dev/null 2>&1
zip -jrq $2 $outdir
rm -rf $tempdir > /dev/null 2>&1
