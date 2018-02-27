#!/bin/bash
#                             $1      $2      $3         $4
# usage: wgs2ncbi_convert.sh '$conf' '$data' '$template' '$outfile'

tempdir=`mktemp -d 2>/dev/null || mktemp -d -t 'tempdir'`
outdir=$tempdir/outdir
datadir=$tempdir/datadir
discrep=$tempdir/discrep.txt
mkdir $outdir
mkdir $datadir
tar -xzf $2 -C $datadir --strip-components=1
wgs2ncbi convert -conf $1 -datadir $datadir -template $3 -outdir $outdir -discrep $discrep -verbose
tar -C $datadir -cvzf $4 . 2>/dev/null
rm -rf $tempdir