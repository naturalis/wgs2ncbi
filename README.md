WGS2NCBI - pipeline for preparing genomes for submission to NCBI
================================================================

The process of going from an annotated genome to a valid NCBI submission is somewhat 
cumbersome. "Boutique" genome projects might produce a scaffolded assembly in FASTA format
and predicted genes in GFF3 tabular format (e.g. produced by the "maker" pipeline) but no 
convenient tools appear to exist to turn these results in the format that NCBI requires.
This project remedies this by providing some Perl scripts (with no dependencies) to do the 
re-formatting. Included is also a shell script that chains the Perl scripts together and
runs NCBI's tbl2asn on the result. This shell script is intended as an example and should
be edited or copied to provide the right values.

Intro: What needs to be submitted
---------------------------------

NCBI requires that "whole genome shotgunning" (WGS) genomes are submitted as .sqn files
(under the new rules of Spring 2013, one file for each scaffold). A sqn file is a file 
in ASN.1 syntax that contains both the sequence, its features, and the metadata about the 
submission (i.e. the authors, the publication title, the organism, etc.). sqn files are 
normally produced by the SeqIn program, which has a graphical user interface and which is 
therefore not practical for the potentially many thousands of files that comprise an 
entire genome. It is therefore preferred to use the program tbl2asn (command line), which 
takes a directory with FASTA files (.fsa) and corresponding files with the gene 
features in tabular format (.tbl), and a submission template (template.sbt) to produce
the sqn files. We therefore need to do some data processing to prepare the inputs for 
tbl2asn.

Since our starting material is one giant FASTA file that contains the scaffolds and a GFF 
file with the features, we need to do the following:

1. use the GenBank web form to create the submission template
2. explode the GFF3 file into smaller ones, one for each scaffold
3. split the FASTA file into scaffolds and feature tables
4. run tbl2asn on the folder with the intermediate files
5. verify output from tbl2asn
6. upload .sqn files

Here now follow more details about each of these steps:

Creating the submission template
--------------------------------

GenBank provides a web form that produces the sbt file. This form needs to be filled out
with the correct metadata, i.e. all the authors of the publication, the publication title,
the organism, etc. The included template.sbt file contains and example. [The form to create
such files is here](http://www.ncbi.nlm.nih.gov/WebSub/template.cgi)

Splitting the GFF3 file
-----------------------

The genome annotation file (GFF3 format) may have the following issues that may prevent
quick lookups of features for a given scaffold:

* very big, so scanning it takes a long time
* include annotation sources we don't trust (e.g. multiple annotation pipelines)
* include features we don't care for (e.g. anything not gene/CDS/3' UTR/5' UTR)
* include FASTA sequence data

To remedy this we "explode" the GFF3 file into separate files, one for each scaffold. This
allows us to quickly find the annotations for a given scaffold (i.e. random access) and we
can filter out included things we don't want. To this end is provided the Perl script
"explode_gff3.pl", which is run as:

 `perl explode_gff3.pl -gff3 $GFF3 -source $SOURCE -f CDS -f gene -f five_prime_UTR \  
 -f three_prime_UTR -d $GFF3DIR`

Where the argument values need to be set to the following:

* GFF3    = the input GFF3 file
* SOURCE  = source of annotations to trust, i.e. the 2nd column in the GFF3, e.g. "maker"
* GFF3DIR = the output directory. This needs to exist already
* the -f <feature> is used multiple times and specifies features to include.

Splitting the FASTA file
------------------------

Included in this archive is a Perl script ("yagc.pl") that takes the big FASTA
file and chops it up into separate FASTA files and tbl files, which it writes into a 
folder called "submission" (or whatever was provided on the command line).

The script is run as:

 `perl yagc.pl -d $TBLDIR -p $PREFIX -f $FASTA -g $GFF3DIR -i $INFO \  
 -a $AUTHORITY -l $LIMIT $VERBOSE`

Where the argument values need to be set to the following:

* TBLDIR    = the output directory to write to. This needs to exist already.
* PREFIX    = a locus_tag prefix that becomes part of transcript/protein IDs, e.g. "OPHA_"
* FASTA     = the input genome in FASTA format
* GFF3DIR   = the directory with exploded annotations
* INFO      = a simple file with key=value pairs, e.g. see info.ini
* AUTHORITY = prefix that identifies a naming authority, e.g. "gnl|NaturalisBC|"
* LIMIT     = optional number of scaffolds to produce, for test runs, e.g. 20
* VERBOSE   = increments the number of progress messages, e.g. "-v" or "-v -v"

The script has no additional dependencies so it should be useable by people who aren't me.

A word of caution: this script produces in some cases tens of thousands of files, each of
which have a name that matches the first word in the FASTA definition line and the *.fsa 
extension. Generally speaking you want to avoid having to look inside the folder that 
contains these files because graphical interfaces (like the windows explorer or the mac 
finder) have a hard time dealing with this.

Running tbl2asn
---------------

Once the submission template, the FASTA files and the feature tables are produced, the
tbl2asn program provided by NCBI needs to be run on the folder that contains these files.
A typical invocation is something like:

 `tbl2asn -p $TBLDIR -t $TEMPLATE -M $MASTERFLAG -a $TYPE -l $LINKAGE -r $ASN1DIR \  
 -Z $DISCREP -V $VERIFY`

Where the argument values need to be set to the following:

* TBLDIR     = the directory that was populated by yagc.pl
* TEMPLATE   = submission template from http://www.ncbi.nlm.nih.gov/WebSub/template.cgi
* MASTERFLAG = the "master genome flag", e.g. "n" for normal
* TYPE       = the file type, e.g. r10k = Runs of 10+ Ns are gaps, 100 Ns are known length
* LINKAGE    = evidence to assert linkage across assembly gaps, e.g. "paired-ends"
* ASN1DIR    = output directory to write to
* DISCREP    = file name of discrepancy report
* VERIFY     = output for verification purposes, e.g. "b" for genbank flat files

When run, this will have produced the required .sqn files. In addition there will be a 
file that reports potential discrepancies ($DISCREP) which needs to be vetted against 
[NCBI's instructions](https://www.ncbi.nlm.nih.gov/genbank/asndisc). Lastly, if additional
$VERIFY arguments were provided (e.g. "b") there will be .gbf files inside $ASN1DIR which
will show (roughly) what the results will look like on the NCBI website.