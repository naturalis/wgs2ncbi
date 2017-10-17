WGS2NCBI - toolkit for preparing genomes for submission to NCBI
===============================================================

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

1. use the GenBank web form to 
   [create the submission template](#creating-the-submission-template)
2. [explode the GFF3 file](#splitting-the-gff3-file) into smaller ones, one for each 
   scaffold: `wgs2ncbi prepare -conf <config.ini>`
3. [split the FASTA file](#splitting-the-fasta-file) into scaffolds and feature tables: 
   `wgs2ncbi process -conf <config.ini>`
4. [run tbl2asn](#running-tbl2asn) on the folder with the intermediate files: 
   `wgs2ncbi convert -conf <config.ini>`
5. verify output from [`tbl2asn`](#running-tbl2asn), make updates to the 
   [protein names](share/products.ini). If need be, keep re-running step 4 and updating 
   the names file until there are no warnings in the descrepancy report.
6. compress .sqn files: `wgs2ncbi compress -conf <config.ini>`, and upload to NCBI. They
   will perform a contaminants screen. If there are suspicious sequences (e.g. untrimmed
   adaptors), mask these using the [adaptors](share/adaptors.ini) file.

In other words, the pipeline mostly consists of invocations of the [wgs2ncbi](script/wgs2ncbi)
script. Each invocation is followed by a verb (prepare, process, convert, compress), followed
by a set of arguments that point to a [configuration file](share/wgs2ncbi.ini), which in
turn points to other files. By perusing the examples of these configuration files you 
should get a pretty good idea how to prepare your own versions of these files. To be able to
run the script, you will need to install it locally. One way to do that is as follows:

1. download the [repository](https://github.com/naturalis/wgs2ncbi/archive/master.zip)
2. unzip it, open a terminal window, and move into the top-level folder
3. `perl Makefile.PL`
4. `sudo make install`

Here now follow more details about each of the steps of the pipeline:

Creating the submission template
--------------------------------

GenBank provides a web form that produces the sbt file. This form needs to be filled out
with the correct metadata, i.e. all the authors of the publication, the publication title,
the organism, etc. The included [template.sbt](template.sbt) file contains an example. 
[The form to create such files is here](http://www.ncbi.nlm.nih.gov/WebSub/template.cgi)

Splitting the GFF3 file
-----------------------

The genome annotation file (GFF3 format) may have the following issues that may prevent
quick lookups of features for a given scaffold:

* very big, so scanning it takes a long time
* include annotation [sources](share/wgs2ncbi.ini#L48) we don't trust (e.g. multiple 
  annotation pipelines)
* include [features](share/wgs2ncbi.ini#L51-L54) we don't care for (e.g. anything not 
  gene/CDS/3' UTR/5' UTR)
* include FASTA sequence data

To remedy this we "explode" the GFF3 file into separate files, one for each scaffold. This
allows us to quickly find the annotations for a given scaffold (i.e. random access) and we
can filter out included things we don't want. This is done using the following command:

    wgs2ncbi prepare -conf <config.ini>

Splitting the FASTA file
------------------------

Once the annotations are exploded, we then need to take the big FASTA file and chop it 
up into multiple FASTA files and [tbl](https://www.ncbi.nlm.nih.gov/projects/Sequin/table.html) 
files, which need to be written into an output folder. The naive behavior is to write 
each scaffold (and its features) to a separate FASTA file. This, however, may result in 
very many files. Therefore, you can provide a parameter to indicate that sequences and 
feature tables are lumped together with up to [`chunksize`](share/wgs2ncbi.ini#L57) 
sequences per file, where `chunksize` may not exceed 10000 according to NCBI guidelines. 
The default for this is 5000.

    wgs2ncbi process -conf <config.ini>

A word of caution: this script produces in some cases tens of thousands of files, each of
which have a name that matches the first word in the FASTA definition line (**so this should
be a unique identifier!**) and the `*.fsa` extension. Generally speaking you want to avoid 
having to look inside the folder that contains these files because graphical interfaces 
(like the windows explorer or the mac finder) have a hard time dealing with this. If you 
use the `chunksize` parameter (which is the default behaviour) the numberof files will be 
a lot lower, and each will have a name matching `combined_xxx-yyy.(fsa|tbl)`, where `xxx` 
and `yyy` are the start and end rank of the sequences in the file.

Running tbl2asn
---------------

Once the submission template, the FASTA files, and the feature tables are produced, the
tbl2asn program provided by NCBI needs to be run on the folder that contains these files.
A typical invocation using the wrapper goes like this:

    wg2ncbi convert -conf <config.ini>

In other words, this command will run the [`tbl2asn`](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) 
command for you with the right command line arguments (provided you have installed it on 
your system and made sure [it can be found](share/wgs2ncbi.ini#L71)). Pay attention to the 
output as this is running, and inspect the [discrepancy report](share/wgs2ncbi.ini#L39). 
This will tell you which protein names are discouraged by NCBI. For these you need to create 
a [mapping](share/products.ini). Note that this conversion step may consequently be an
iterative process: if the descrepancy report raises issues about names you will need to 
address and re-run the step.

Uploading to NCBI
-----------------

Tip: Note that NCBI **does** accept .tar.gz archives, which means you can prepare your
[upload](share/wgs2ncbi.ini#L42) as follows:

    wgs2ncbi compress -conf <config.ini>

Once you upload the archive, you will get a verdict from whoever is handling this submission
at NCBI. It is possible that there will be stretches of sequence in your submission that
NCBI will consider contaminants (based on a pipeline they run). One way to deal with those
is to blank them out of the data, using a [configuration file](share/adaptors.ini) that
specifies the coordinates of stretches to NNN. In addition, NCBI might have additional
issues with certain protein names, so you may have to update the 
[names mapping file](share/products.ini). Then rerun the convert step, rebuild the archive,
and do another upload.
