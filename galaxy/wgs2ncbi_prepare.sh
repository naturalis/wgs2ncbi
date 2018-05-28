#!/bin/bash
# usage: wgs2ncbi_prepare.sh $GFF3 $CONF_INI $OUTFILE

gff3dir=`mktemp -d 2>/dev/null || mktemp -d -t 'gff3dir'`
wgs2ncbi prepare -conf $2 -gff3dir $gff3dir -gff3file $1
zip -jrq $3.zip $gff3dir
mv $3.zip $3
rm -rf $gff3dir
