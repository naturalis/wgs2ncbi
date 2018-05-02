#!/bin/bash
#                             $1      $2          $3
# usage: wgs2ncbi_convert.sh '$args' '$out_file1' '$data'

tempdir=`mktemp -d 2>/dev/null || mktemp -d -t 'tempdir'`
outdir=$tempdir/outdir
datadir=$tempdir/datadir
mkdir $outdir
mkdir $datadir
tar -xzf $2 -C $datadir --strip-components=1
wgs2ncbi convert $1 -datadir $datadir -outdir $outdir
tar -C $datadir -cvzf $2 . > /dev/null 2>&1
rm -rf $tempdir > /dev/null 2>&1
