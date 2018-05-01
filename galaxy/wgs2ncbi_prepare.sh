#!/bin/bash
# usage: wgs2ncbi_prepare.sh $GFF3 $CONF_INI $OUTFILE

gff3dir=`mktemp -d 2>/dev/null || mktemp -d -t 'gff3dir'`
wgs2ncbi prepare -conf $2 -gff3dir $gff3dir -gff3file $1
tar -C $gff3dir -cvzf $3 . > /dev/null 2>&1
rm -rf $gff3dir > /dev/null 2>&1
