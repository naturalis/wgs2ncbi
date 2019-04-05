#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::WGS2NCBI::Logger;
use constant EXCLUDE => 1;
use constant TRIM    => 2;
use constant DEDUP   => 3;

=pod

[1] We ran the sequences through the contamination screen.  The screen found 110 sequences
 to exclude, 103 sequences with locations to mask/trim (see attached).

- fixed by reading in the report file and excluding/trimming the sequences it includes

[2] In addition, the contamination screen found 48 duplicated sequences (see attached).  
We report these back in case they are an artifact of the assembler and should be removed.  
If they belong in the assembly, then that's fine.

- fixed by reading in the report file and excluding/trimming the sequences it includes

[3] Please remove any contigs that are shorter than 200 nucleotides since they generally 
add little information. We can make an exception if a contig is important for the assembly 
because it is part of a large multi-component scaffold. However, we would not want short 
contigs that were only mate pairs or that were singletons.

- fixed by adding a $minlength option below which sequences are ignored

[4] There are several short sequences that are more than 5% N's in the discrepancy report 
(see report uploaded to submission portal).  If these are low quality sequences, please 
remove them.

- won't do this: NNN's are gaps because of paired-end assembly

[5] There are some short introns that appear to be frameshifts.  We prefer not to include 
translations for CDS features when the translation is known to be incorrect.  Please 
remove the CDS/mRNA features and annotate this with a single gene feature across the 
entire span.  Include a note on the gene with a brief description.  For example:

- won't do this: we don't think they're frameshifts

[6] There are a few exons that begin or end in a gap.  Are you certain you have the 
correct translation (for example, you aren't missing an exon in the middle of the 
translation)?  If the translations are not correct, split the features so you have two 
partial CDS features (and mRNAs) abutting the gap, with a single gene over the whole 
locus.  Alternatively, one of the partial CDS/mRNA features may be deleted if it is very 
short and there is little or no supporting evidence.  If you have a single gene and two 
partial CDS/mRNA features, you should: (1) add a note to each CDS referencing the other 
half of the gene, (2) add a note to the gene and CDS features stating, "gap found within 
coding sequence."

- FIXME

[7] 1444 of the CDS features have invalid translations (no stop codon). Every CDS feature 
should have a valid start and stop codon, and should not have any internal stops.  The 
only exception is if the CDS is partial at the end of a sequence or at an intron/exon 
boundary.  If the CDS is partial, it must have the appropriate partial symbols.

- the CDS features without stop codons are annotated as partial.

[8] Is this OK or should this gene be removed?:

WARNING: valid [SEQ_FEAT.MultipleGeneOverlap] Gene contains 5 other genes FEATURE: Gene: 
PER2 <34689> [lcl|scaffold1852.1-size148322:<33631-138407] [lcl
|scaffold1852.1-size148322: delta, dna len= 148322]

- FIXME

[9] If you have enough information to add a gene name, we would expect that you could have 
a better product name.

- fixed, we now use the gene name for the product, if it exists

[10] Some of the product names need to be improved (see the discrepancy report output 
uploaded to submission portal).  Remove any database identifiers, locus_tags or systematic 
names.  We only want to use biological product names.

- these are genuine gene names, not identifiers.

[11] It would be better if you removed the size information from the contig names.

- fixed

[12] When you submit your files, include as much source information as you can in the 
fasta header (eg strain, host, collection-date, isolation-source, country where the sample 
was collected).

- the fasta header contains the only applicable source information: the organism name

=cut

my $report    = 'share/contamination_dups';
my $fasta     = 'share/king_cobra_scaffolds_spring_2011.fasta';
my $gff3      = 'share/cobra.functional.gff';
my $outfasta  = 'share/king_cobra_scaffolds_spring_2011-cleaned.fasta';
my $outgff3   = 'share/cobra.functional-cleaned.gff';
my $minlength = 200;

GetOptions(
	'verbose+'    => \$Bio::WGS2NCBI::Logger::Verbosity,
	'report=s'    => \$report,
	'fasta=s'     => \$fasta,
	'gff3=s'      => \$gff3,
	'outfasta=s'  => \$outfasta,
	'outgff3=s'   => \$outgff3,
	'minlength=i' => \$minlength,
);

# read the NCBI report with sequence IDs to exclude or trim
my ( %exclude, %trim );
{
	my $section = 0;
	open my $fh, '<', $report or die $!;
	INFO "going to read $report";
	LINE: while(<$fh>) {
		chomp;
		if ( /(?:Exclude|Trim|DuplicatedContent):\s*$/ ) {
			$section++;
			next LINE;
		}
	
		# only need to record contig name
		if ( $section == EXCLUDE && /^(\S+)/ ) {
			my $contig = $1;
			$contig =~ s/-.+//;
			$exclude{$contig} = 1;
			INFO "going to exclude $contig";
		}
	
		# need to record coordinates to trim
		elsif ( $section == TRIM && /^(\S+)\t\d+\t(\d+)\.\.(\d+)/ ) {
			my ( $contig, $start, $stop ) = ( $1, $2, $3 );
			$contig =~ s/-.+//;
			$trim{$contig} = [ $start, $stop ];
			INFO "going to trim range $start .. $stop from $contig";
		}
	
		# need to record other contigs to dedup
		elsif ( $section == DEDUP && /\S/ ) {
			my @fields = split /\s+/;
			for my $contig ( @fields[1..$#fields - 2] ) {
				$contig =~ s/^lcl\|//;
				$contig =~ s/-.+//;
				$exclude{$contig} = 1;
				INFO "going to remove duplicate $contig";
			}
		}
	}
	INFO "done reading $report";
}

# read the fasta file
{
	my ( $id, $seq );
	INFO "going to clean sequences from $fasta";
	open my $fh, '<', $fasta or die $!;
	open my $out, '>', $outfasta or die $!;
	while(<$fh>) {
		chomp;
		if ( />(.+)/ ) {
			finalize_seq( $id, \$seq, $out ) if $id and $seq;
			$id = $1;
			$id =~ s/\|.+//;
			$seq = '';
		}
		else {
			$seq .= $_;
		}
	}
	finalize_seq( $id, \$seq, $out );
}

# exclude and trim sequences, fix contig names
sub finalize_seq {
	my ( $id, $seqref, $fh ) = @_;
	INFO "finalizing $id";
	my $ncount = $$seqref =~ tr/N/N/;
	if ( ($ncount/length($$seqref)) > 0.05 ) {
		WARN "$id has more than 5% gaps";
	}
	if ( $minlength > length($$seqref) ) {
		WARN "$id is too short, skipping";
		$exclude{$id} = 1; # this so that we also skip the annotation
		return;
	}
	if ( $trim{$id} ) {
		INFO "going to trim contig $id";
		my ( $start, $stop ) = @{ $trim{$id} };
		my $length = $stop - $start + 1;
		substr $$seqref, $start - 1, $length, 'N' x $length;
	}
	if ( not $exclude{$id} ) {
		INFO "writing $id";
	
		# print basic header
		print $fh '>', $id, "\n";
	
		# fold lines at 80 characters
		print $fh join "\n", unpack "(a80)*", $$seqref;
		
		# end the sequence record
		print $fh "\n";
	}
}

# read the annotation file, exclude and trim annotations, fix contig names
{
	my $id_i = 0;

	open my $fh, '<', $gff3 or die $!;
	INFO "going to adjust annotation coordinates for $gff3";
	open my $out, '>', $outgff3 or die $!;
	LINE: while(<$fh>) {
		next LINE if not /\t/;
		my @fields = split /\t/;
		my $id = $fields[$id_i];
		$id =~ s/-.+//;
		next LINE if $exclude{$id};
		
		# we can now fix the contig name
		$fields[$id_i] = $id;
		
		# print result contig names
		print $out join("\t", @fields), "\n";
	}
}
