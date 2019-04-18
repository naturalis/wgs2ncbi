[![DOI](https://zenodo.org/badge/11924393.svg)](https://zenodo.org/badge/latestdoi/11924393)
[![status](http://joss.theoj.org/papers/72b726934646d999001e9997ad2fdb54/status.svg)](http://joss.theoj.org/papers/72b726934646d999001e9997ad2fdb54)

WGS2NCBI - preparing genomes for submission to NCBI
===================================================

The process of going from an annotated genome to a valid NCBI submission is somewhat 
cumbersome. "Boutique" genome projects typically produce a scaffolded assembly in FASTA 
format as produced by any of a variety of de-novo assemblers and predicted genes in GFF3 
tabular format, e.g. as produced by the 
[maker](http://www.yandell-lab.org/software/maker.html) pipeline, but no convenient tools 
appear to exist to turn these results in a format and to a standard that NCBI accepts.

NCBI requires that "whole genome shotgunning" (WGS) genomes are submitted as `.sqn` files. 
A sqn file is a file in ASN.1 syntax that contains both the sequence, its features, and 
the metadata about the submission, i.e. the authors, the publication title, the organism, 
etc.. `.sqn` files are normally produced by the 
[sequin](https://www.ncbi.nlm.nih.gov/Sequin/) program, which has a graphical user 
interface. Sequin works fine for a single gene or for a small genome (e.g. a mitochondrial 
genome) but for large genomes with thousands of genes spread out over potentially 
thousands of scaffolds the submission process done in this way would be unworkable.

The alternative is to use the [tbl2asn](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) 
command line program, which takes a directory with FASTA files (`.fsa`) and corresponding 
files with the gene features in tabular format (`.tbl`), and a submission template 
(template.sbt) to produce the `.sqn` files. The trick thus becomes to convert the assembly
FASTA file and the annotation GFF3 file into a collection of FASTA chunks with
corresponding feature tables. This is doable in principle - several toolkits provide
generic convertors -, but NCBI places quite a few restrictions on what are permissible 
things to have in the FASTA headers, what coordinate ranges are credible as gene features,
and what gene and gene product names are acceptable.

This project remedies these challenges by providing a command-line utility (with no 3rd 
party dependencies except [URI::Escape](http://search.cpan.org/dist/URI-Escape)) to do 
the required data re-formatting and cleaning. Included is also a shell script that chains 
the Perl scripts together and runs NCBI's tbl2asn on the result. This shell script is 
intended as an example and should be edited or copied to provide the right values.

Installation
============

The WGS2NCBI release is organized in a way that is standard for software releases written
in the Perl5 programming language. This means that it can be installed using a series
of commands that either you yourself, or your systems administrator, are likely already
familiar with. The first step is to install a required dependency using the Perl5 package
manager ([cpan](https://perldoc.perl.org/cpan.html)), as follows:

    $ sudo cpan -i URI::Escape
    
The next steps assume that you have downloaded the WGS2NCBI release - for example from the
[git repository](https://github.com/naturalis/wgs2ncbi/archive/master.zip) - have unzipped
it, and have moved into the root folder of the release in your terminal. The next steps
then are as follows:

    $ perl Makefile.PL
    $ make test
    $ sudo make install

The second command (`make test`) performs a number of basic tests of the software on your
system. These should all pass without problems. If you do encounter issues, it is best
_not_ to proceed to the following step for the actual installation, but rather to try to
resolve the outstanding problems, for example by submitting an
[issue report](https://github.com/naturalis/wgs2ncbi/issues), so that the authors can help
you out.

In addition to the preceding steps, you also need to install the `tbl2asn` program. The 
instructions for this are [here](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/).

Usage
=====

WGS2NCBI is used by following a number of steps, which are detailed below:

- [Before you start](#before-you-start) - set up all the input files, prepare a submission
  template
- [Subcommand `prepare`](#subcommand-prepare) - pre-process the annotation file for rapid
  access in the following steps
- [Subcommand `process`](#subcommand-process) - convert the genome file and annotations
  to FASTA chunks and feature tables
- [Subcommand `convert`](#subcommand-convert) - runs tbl2asn to convert the FASTA chunks
  and feature tables to SeqIn files
- [Subcommand `compress`](#subcommand-compress) - collates the SeqIn files into a single
  archive for upload to NCBI      

Before you start
----------------

Before issuing any commands, the following steps need to be taken:

1. The installation (see above) needs to be completed.
2. You need to have the genome assembly available as a FASTA file, and the annotations
   as a GFF3 file.
3. You will need to prepare a 
   [submission template](http://www.ncbi.nlm.nih.gov/WebSub/template.cgi). The file
   [template.sbt](share/template.sbt) is an example of what these files look like.
4. You need to have created a number of `.ini` files correctly. Using the linked files as 
   examples, the following need to be prepared:
   - [wgs2ncbi.ini](share/wgs2ncbi.ini) - the main configuration file, in which you 
     specify the locations of the input files and output directories. In addition, here
     you will specify the prefixes for the identifiers that will be inserted in the 
     feature tables and various parameters for what to filter on. The file is well 
     documented with comments.
   - [info.ini](share/info.ini) - a file with key/value pairs whose contents will be 
     inserted in the FASTA headers of the sequence files. These key/value pairs have to
     do with the organism that was sequenced, such as the taxon name, its sex, its
     developmental stages, what tissues were sampled, and so on.
   - [adaptors.ini](share/adaptors.ini) - this is a file that contains the coordinates 
     of sequence fragments that NCBI considers inadmissible. What will happen over the
     course of your submission is that NCBI will scan your sequence data for suspicious
     sequence fragments. These might be adaptor sequences of various sequencing platforms,
     and fragments that NCBI thinks might be contaminants. Hence, during your first pass
     it is more or less impossible to get the values right in this file: this part will
     be an iterative process where you blank out parts of your data that NCBI really will
     not accept. Start out with an empty file, and populate it based on the feedback you
     will get, making sure you follow the same syntax as the provided example file.
   - [products.ini](share/products.ini) - this is a file that contains mappings from 
     (parts of) the gene names that you assigned during the annotation process to names
     that NCBI will accept. Again, this is impossible to predict during the first pass:
     you will get feedback on which names NCBI doesn't like (for example because there are
     things in the names that look like database identifiers, organism names, molecular
     weights, etc.) and in this file you map these to allowed names.

Subcommand `prepare`
--------------------

Once the preparation is done, you will now run the `prepare` subcommand, as follows:

    $ wgs2ncbi prepare -conf <wgs2ncbi.ini>

The value of the `-conf` argument specifies the location of the 
[wgs2ncbi.ini](share/wgs2ncbi.ini) configuration file. In all following steps you will
also need to provide the location of this same file. 

What happens during this step is that the GFF3 file is pre-processed so that the following
steps will have quicker access to the relevant contents than they would have if they had
to scan through the entire file every time. To be precise, the following happens:

- only 'true' annotation data is retained. GFF3 files may also contain their own bits of
  fasta data, but these are filtered out.
- only the annotations produced by the specified annotation
  [source](https://github.com/naturalis/wgs2ncbi/blob/master/share/wgs2ncbi.ini#L48), e.g.
  the `maker` pipeline, are retained.
- only those features specified under 
  [feature](https://github.com/naturalis/wgs2ncbi/blob/master/share/wgs2ncbi.ini#L51-L54)
  are retained.
- the remaining data are written to the
  [gff3dir](https://github.com/naturalis/wgs2ncbi/blob/master/share/wgs2ncbi.ini#L32), one
  file for each contig.

As such, this step is to do initial filtering and pre-processing. Typically you will only
need to run this step once.

Subcommand `process`
--------------------

This subcommand is issued as follows:

    $ wgs2ncbi process -conf <wgs2ncbi.ini>

i.e. by providing the location of the [wgs2ncbi.ini](share/wgs2ncbi.ini) configuration 
file to the `-conf` argument.

The `process` subcommand contains most of the "intelligence" (such as it is). In this step
the following happens:

- the genome assembly, i.e. the large FASTA file, is chopped up into smaller FASTA files.
  All but the last of these output files will contain as many FASTA records as specified
  by [chunksize](https://github.com/naturalis/wgs2ncbi/blob/master/share/wgs2ncbi.ini#L57)
  with the last one containing the remainder. If all your contigs are longer than the
  [minlength](https://github.com/naturalis/wgs2ncbi/blob/master/share/wgs2ncbi.ini#L60)
  then the number of files thus produced will be the number of contigs, divided by 
  `chunksize`, rounded up to the nearest integer. However, contigs smaller than
  `minlength`, if you have them, will be omitted, as NCBI won't accept these.
- the FASTA data that will be written will have any stretches specified in
  [adaptors.ini](share/adaptors.ini) replaced with `NNNs`. These will be sequence 
  fragments that NCBI will specify as inadmissible because they might be sequence adaptors 
  (i.e. vendor-specific synthetic DNA) or contaminants.
- the FASTA files will have the `.fsa` file extension, as required by `tbl2asn`.  
- the annotations from the GFF3 file, pre-processed in the previous step, will be written
  out as feature tables (required extension: `.tbl`). There will be as many `.tbl` files 
  as there are `.fsa` files.
- any gene annotations that have introns that are shorter than
  [minintron](https://github.com/naturalis/wgs2ncbi/blob/master/share/wgs2ncbi.ini#L67)
  will be converted to pseudogenes, as NCBI does not believe these could be real.
- any gene product names that are unacceptable to NCBI, and for which you have provided 
  a mapping in [products.ini](share/products.ini), will be mapped to the names
  you have provided.  
- both the `.fsa` and the `.tbl` files will be written in the same directory, specified by
  [datadir](https://github.com/naturalis/wgs2ncbi/blob/master/share/wgs2ncbi.ini#L26)

Since the results of this step depend on the settings in the 
[products.ini](share/products.ini) and [adaptors.ini](share/adaptors.ini) files, and since
you will hear from NCBI what needs to go in these files, this step and the following ones
are something that you will probably run multiple times until NCBI is happy.

Subcommand `convert`
--------------------

This subcommand will run `tbl2asn`. As such, it is essential that this program is 
installed successfully according to NCBI's instructions, which are
[here](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/). If the program is installed such
that it is available on the [PATH](http://www.linfo.org/path_env_var.html) you can 
proceed with this step without making any changes. If you've had to install it in a 
location where it is not on the `PATH`, you can specify an alternative location
under [tbl2asn](https://github.com/naturalis/wgs2ncbi/blob/master/share/wgs2ncbi.ini#L71)
in the main configuration file. Check that the executable is ready to run, e.g. by
issuing `which tbl2asn` if it is on the `PATH` or by running it from its alternative 
location.

Once you are all set, issue the subcommand as follows:

    $ wgs2ncbi convert -conf <wgs2ncbi.ini>

i.e. by providing the location of the [wgs2ncbi.ini](share/wgs2ncbi.ini) configuration 
file to the `-conf` argument. The following will then happen:

- the Sequin files, with the `.sqn` extension, will be written to the directory
  [outdir](https://github.com/naturalis/wgs2ncbi/blob/master/share/wgs2ncbi.ini#L36)
- the discrepancy report, containing all the problems that `tbl2asn` diagnosed, will be
  written to 
  [discrep](https://github.com/naturalis/wgs2ncbi/blob/master/share/wgs2ncbi.ini#L39)

Note that this discrepancy report will give you the first suggestions for problematic
gene product names (which you can deal with in [products.ini](share/products.ini)), but
this will not be exhaustive: NCBI will likely point out additional problems, and any
checks for contaminations or spurious adaptors will only be performed by NCBI. In other
words, citing their website:

> The Discrepancy Report is an evaluation of a single or multiple ASN.1 files, looking for 
> suspicious annotation or annotation discrepancies that NCBI staff has noticed commonly 
> occur in genome submissions, both complete and incomplete (WGS). A few of the problems 
> that this function was written to find include inconsistent locus_tag prefixes, missing 
> protein_id's, missing gene features, and suspect product names. The function is 
> available in specially configured Sequin, as an argument for tbl2asn, or with the 
> command-line program asndisc.
>
> If you have questions about the Discrepancy Report, please contact us by email at 
> genomes@ncbi.nlm.nih.gov prior to sending us your submission.
Source: https://www.ncbi.nlm.nih.gov/genbank/asndisc/

Subcommand `compress`
---------------------

The final step simply takes the `.sqn` files from the previous step and combines them in
a single `.tar.gz` archive for upload to the NCBI submission portal. No data processing of
any kind takes place, this is purely for convenience and is executed as follows:

    $ wgs2ncbi compress -conf <config.ini>
    
i.e. by providing the location of the [wgs2ncbi.ini](share/wgs2ncbi.ini) configuration 
file to the `-conf` argument. The following will then happen:

- all .sqn files are combined in a single archive, whose location is specified by
  [archive](https://github.com/naturalis/wgs2ncbi/blob/master/share/wgs2ncbi.ini#L42)

You will then upload the produced archive to the submission portal. Once you upload the 
archive, you will get a verdict from whoever is handling this submission at NCBI. 
Depending on their feedback, you will likely have to update the configuration files a few
more times to correct for spurious sequence data and gene product names, after which you
will re-run the `process` subcommand (and onwards to `convert` and `compress`).

About this software
===================

WGS2NCBI is implemented as a Perl5 package. It is open source software made available
under the [BSD3 license](LICENSE).

If you experience any difficulties with this software, or you have suggestions, or want
to contribute directly, you have the following options:

- submit a bug report or feature request to the 
  [issue tracker](https://github.com/naturalis/wgs2ncbi/issues)
- contribute directly to the source code through the 
  [github](https://github.com/naturalis/wgs2ncbi) repository. 'Pull requests' are
  especially welcome.
