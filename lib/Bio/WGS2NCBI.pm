package Bio::WGS2NCBI;
use strict;
use warnings;
use Pod::Usage;
use URI::Escape;
use Archive::Tar;
use Bio::WGS2NCBI::Seq;
use Bio::WGS2NCBI::Logger;
use Bio::WGS2NCBI::Config;
use Bio::WGS2NCBI::Feature;
use Bio::WGS2NCBI::FeatureSet;
use Bio::WGS2NCBI::TableReader;
use File::Path 'make_path';

# version number for the whole release is defined here
our $VERSION = 'v0.1.2';

# export the run function
require Exporter;
use base 'Exporter';
our @EXPORT = qw(process prepare compress convert prune);

# GFF3 column numbers
my $chr_idx    = 0;
my $source_idx = 1;
my $type_idx   = 2;
my $start_idx  = 3;
my $end_idx    = 4;
my $strand_idx = 6;
my $codon_idx  = 7;
my $meta_idx   = 8;

=pod

=head1 NAME

Bio::WGS2NCBI - module to assist in submitting whole genome sequencing projects to NCBI

=head1 DESCRIPTION

This module documents the four actions (prepare, process, convert and compress) that are
available to users of the L<wgs2ncbi> script. Each of these steps is configured by one or
more configuration files. In the documentation below, the relevant fields from these
configuration files are listed. To understand how the configuration system itself works,
consult the documentation at L<Bio::WGS2NCBI::Config>.

=head1 prepare

The C<prepare> action takes the annotation file (in GFF3 format) and extracts the relevant
information out of it, writing it to a potentially large set of files. This is done 
because GFF3 annotation files can become quite large, so that finding the annotations 
for any particular scaffold might take a long time if this is done by scanning through
the whole file. Instead, the set of annotations is reduced, by taking the following steps:

=over

=item B<Remove all sequence data from the file> - Embedded FASTA data is permissible 
according to the GFF3 standard, but this makes files needlessly bulky if the FASTA data 
are available separately as well.

=item B<Remove all annotations from unrecognized sources> - GFF3 files can contain 
annotations from sources you don't particularly trust and want to ignore in your 
submission. This is configurable.

=item B<Remove all irrelevant features> - any sequence features in GFF3 that are not 
recognized in NCBI feature tables are discarded. This is configurable.

=back

Subsequently, the annotations are split such that there is a separate file for each 
scaffold (or chromosome, if that is how your annotations are organized). This way, the
relevant information for any scaffold can be found much more quickly through the file 
system rather than by scanning through a large file.

In order for this action to succeed, the following configuration values need to be 
provided:

=over

=item C<gff3file>

The location of the input annotation file in GFF3 format.

=item C<gff3dir>

The location of the output directory for the split annotation files.

=item C<source>

Which annotation source to trust.

=item C<feature>

Which features to retain.

=back

=cut

sub prepare {
	my $config   = Bio::WGS2NCBI::Config->new;
	my %types    = map { $_ => 1 } $config->feature;
	my $gff3file = $config->gff3file;
	my $gff3dir  = $config->gff3dir;
	my $source   = $config->source; # XXX maybe comma separated list?
	my $chr_idx  = 0;
	my $src_idx  = 1;
	my $type_idx = 2;
	my @queue;
	my %handles;
	
	# check the configuration
	my $quit;
	if ( not $gff3file ) {
		ERROR "No GFF3 file provided.";
		$quit = 1;
	}
	elsif ( not -e $gff3file ) {
		ERROR "No GFF3 file at location '$gff3file'";
		$quit = 1;
	}
	if ( not $gff3dir ) {
		ERROR "No GFF3 directory provided.";
		$quit = 1;	
	}
	elsif ( not -d $gff3dir ) {
		ERROR "No GFF3 directory at location '$gff3dir'";
		$quit = 1;	
	}
	if ( $quit ) {
		ERROR "Quitting. Try 'wgs2ncbi help' for more info.";
		exit(1);
	}

	# open file handle
	INFO "Going to pre-process annotations in $gff3file";
	open my $fh, '<', $gff3file or die;
	
	# iterate over lines
	my $line = 0;
	LINE: while(<$fh>) {
		next LINE if /^#/;
		my @line = split /\t/, $_;
		
		# we will automatically skip over sequence lines as they won't have columns
		next LINE if not $line[$src_idx] or $line[$src_idx] ne $source;
		next LINE if not $line[$type_idx] or not $types{$line[$type_idx]};
		my $chr = $line[$chr_idx];
	
		# we keep up to 100 file handles (which is max on POSIX). once we
		# reach that number we close existing ones on a first in, first out basis
		if ( not $handles{$chr} ) {
			if ( scalar @queue == 100 ) {
				my $fifo = shift @queue;
				close $handles{$fifo};
				delete $handles{$fifo};
			}
			open $handles{$chr}, '>', "${gff3dir}/${chr}.gff3" or die $!;
			push @queue, $chr;
		}
		my $fh = $handles{$chr};
		print $fh $_;
		
		# log line #
		if ( not ++$line % 10000 ) {
			INFO "Processed line $line";
		}
	}
	INFO "Done. Files written to $gff3dir";
}

=head1 process

The C<process> action takes the sequence file (in FASTA format, with a record for each
scaffold or chromosome) and the pre-processed annotations and converts these into feature
table files and (masked) sequence files.

In order for this action to succeed, the following configuration values need to be 
provided:

=over

=item C<datafile>

The location of the input FASTA file, with a record for each scaffold or chromosome.

=item C<info>

An INI-style configuration file that contains key/value pairs that will be embedded in the
FASTA sequence headers of the produced output files. This is typically used for metadata
about the sampled organism, such as its sex, collection locality, collected cell type, 
etc.

=item C<masks>

An INI-style configuration file that contains the coordinates of sequence segments to 
mask. This may be needed because NCBI will do a strict screen to check for unclipped 
adaptor sequences or contaminants. In the report that is returned by NCBI it will state 
the sequence coordinates of segments that NCBI will not accept in a submission. By putting 
these coordinate in this file the offending segments will be replaced with NNNs.

=item C<products>

An INI-style configuration file that contains the corrected names for protein products. 
The rationale is that your genome annotation process may introduce protein names that NCBI
would like to deprecate, such as names that include molecular weights, database 
identifiers, references to 'homology', and so on. The discrepancy report that is produced
in the L<Bio::WGS2NCBI/convert> step will be a first guide in composing corrected names,
but the validation that NCBI will perform will likely point out additional errors.

=item C<gff3dir>

The location of the pre-processed GFF3 files as produced by L<Bio::WGS2NCBI/prepare>.

=item C<datadir>

The location of the output dir where the (potentially 'chunked', see below) sequence files
and feature tables will be written.

=item C<prefix>

A short character sequence that is prefixed to every sequence record identifier that is
generated. NCBI will provide submitters with this prefix when the submission is 
initialized.

=item C<authority>

This is a naming authority that will be applied to all sequence record identifiers. A
reasonable value for this could be the name of the lab or institution that leads the 
project resulting in the submission. NCBI intends this authority, in combination with the
C<prefix> as a way to ensure that sequences are globally uniquely identifiable.

=item C<minlength>

The minimum length of a scaffold to be retained in a submission. This should be 200 or 
above.

=item C<minintron>

The minimum length of an intron to be retained in a submission. Introns shorter than this
are interpreted (by NCBI) to be spurious and should therefore by discarded. As a 
consequence, the gene that contains such an intron will be annotated as a pseudogene. This
value must be 10 or above.

=item C<chunksize>

The output that is produced can be combined into chunks of more than one scaffold per 
file. To keep the number of files manageable it is convenient to set this to a large 
value, but less than or equal to 10,000.

=item C<limit>

This parameter allows you to run the process on only a limited set of scaffolds. This is
provided for testing, "dry run" purposes. For real usage this value must be set to 0.

=back

=cut

sub process {
    my $config = Bio::WGS2NCBI::Config->new;
    
    # read the info file, if any
    my %info = $config->read_ini( $config->info );
    my $desc = join ' ', map { sprintf('[%s=%s]',$_,$info{$_}) } keys %info;
    
    # read the masks, if any
    my %masks;
    if ( my $file = $config->masks ) {
	
        # read the masks file, make values uniformly array refs
        %masks = $config->read_ini( $file );
        for my $key ( keys %masks ) {
            $masks{$key} = [ $masks{$key} ] if not ref $masks{$key};
        }
        INFO "Read masks for ".scalar(keys(%masks))." sequences";
    }
    
    # open the fasta file handle
    open my $fastaFH, '<', $config->datafile or die $!;
    my $fastaPos = 0;   
    
    # counter to generate autoincrementing IDs for gene, mrna and protein
    my $counter = {};
    my $seq_counter;
    
    # we will re-use these handles
    my ( $tblFH, $scaffoldFH );
    
    # iterate over the genome
    SEQ: while( not eof $fastaFH ) {
    	++$seq_counter;
        last if $config->limit and $seq_counter == $config->limit;
    
        # advance to the next scaffold/chromosome 
        my ( $offset, $seq ) = ( 0 );
        ( $fastaPos, $seq ) = Bio::WGS2NCBI::Seq->read_fasta( $fastaFH, $fastaPos );
        $seq->desc($desc);
        my $length = $seq->length;
        my $chr    = $seq->id;
                
        # compute offset at beginning, if any. e.g. if NNCGTNN, $offset is 2
        if ( $offset = $seq->get_non_missing_index ) {
            INFO "Leading ${offset}bp gap in $chr, will strip this and apply offset";
        }
        
        # check what we have left, here if NNCGTNN, $lnmi == 4
        my $last_non_missing_index = $seq->get_non_missing_index('reverse');
		$length -= ( $length - 1 - $last_non_missing_index ) + $offset;
		if ( $length < $config->minlength ) {
			WARN "Remaining seq $chr is too short ($length bp), skipping";
			$seq_counter--;
			next SEQ;
		}
        
        # open handle to the features table
        if ( ( $seq_counter % $config->chunksize ) == 1 ) {
            
            # close if already open
            if ( $tblFH and $scaffoldFH ) {
                close $tblFH;
                close $scaffoldFH;
            }
        
            # generate a new name indicating the range
            my $upper = $seq_counter + $config->chunksize - 1;
            my $stem = "combined_${seq_counter}-${upper}";
            open $tblFH, '>', $config->datadir . "/${stem}.tbl" or die $!;   
            open $scaffoldFH, '>', $config->datadir . "/${stem}.fsa" or die $!;     
        }
        
        # get the features for that scaffold/chromosome, if we have them
        my $gff3 = $config->gff3dir . "/${chr}.gff3";
        if ( -e $gff3 ) {
            open my $gff3FH, '<', $gff3 or die $!;  
            my $features = Bio::WGS2NCBI::FeatureSet->read_features( 
            	$gff3FH,  # file handle of focal file
            	$counter, # counter for generating locus tags
            	$config,  # config object
            	$seq,     # sequence object
            	$offset,  # start offset
            	$last_non_missing_index, # last true seq character
            );  
            $features->offset($offset);     
            $features->write_features( $tblFH );        
        }
        else {
            print $tblFH '>Features ', $chr, "\n";
        }
        
		# truncate, mask, serialize
        my $trunc = $seq->trunc( $offset + 1, $offset + $length );
        $trunc->mask( @{$masks{$chr}} ) if $masks{$chr};        
        $trunc->write_fasta( $scaffoldFH );        
    }   
}

=head1 convert

The C<convert> action runs the C<tbl2asn> program provided by NCBI with the right 
settings. This requires the following configuration settings:

=over

=item C<datadir>

The location of the dir where the (potentially 'chunked', see below) sequence files
and feature tables were written by L<Bio::WGS2NCBI/process>.

=item C<template>

The location of the template file produced with the form at:
L<http://www.ncbi.nlm.nih.gov/WebSub/template.cgi>

=item C<outdir>

The location where to write the resulting ASN.1 files.

=item C<discrep>

The location where to write the discrepancy report.

=item C<tbl2asn>

The location where the tbl2asn executable is located.

=back

=cut

sub convert {
	my $config = Bio::WGS2NCBI::Config->new;
	my $INDIR   = $config->datadir;
	my $TMPL    = $config->template;
	my $OUTDIR  = $config->outdir;
	my $DISCREP = $config->discrep;
	my $TBL2ASN = $config->tbl2asn;	
	my $command = "$TBL2ASN -p $INDIR -t $TMPL -M n -a r10k -l paired-ends -r $OUTDIR -Z $DISCREP > /dev/null 2>&1";
	INFO "going to execute command '$command'";
	exec $command;
}

=head1 trim

The C<trim> action trims stretches of leading or trailing NNNs from sequence records, and
updates the coordinates in the associated feature tables accordingly. In cases where a 
feature falls within a trimmed region, the feature is removed entirely.

=over

=item C<datadir>

The location of the dir where the (potentially 'chunked', see below) sequence files
and feature tables were written by L<Bio::WGS2NCBI/process>.

=back

=cut

sub trim {
	my $config = Bio::WGS2NCBI::Config->new;
	my $INDIR  = $config->datadir;
	
	# iterate over files in folder, read FASTA files
	opendir my $dh, $INDIR or die $!;
	while( my $file = readdir $dh ) {
		
		# have a FASTA file
		if ( $file =~ /(.+)\.fsa$/ ) {
			my $stem = $1;
		
			# make backup of FASTA file
			rename "${INDIR}/${file}", "${INDIR}/${file}.bak";
			
			# read file, look op non-missing residue positions, write truncated
			open my $fh,  '<', "${INDIR}/${file}.bak" or die $!;
			open my $out, '>', "${INDIR}/${file}"     or die $!;
			my ( $pos, $seq, %coord );
			while( not eof($fh) ) {
				( $pos, $seq ) = Bio::WGS2NCBI::Seq->read_fasta( $fh, $pos );
				my $id = $seq->id;
				my $i1 = $seq->get_non_missing_index;
				my $i2 = $seq->get_non_missing_index(1);
				INFO "$id\t$i1 .. $i2";
				$coord{$id} = [ $i1, $i2  ];
				$seq->trunc( $i1 + 1, $i2 + 1 )->write_fasta($out);	
			}
			
			# make backup of TBL file, open handle for writing		
			rename "${INDIR}/${stem}.tbl", "${INDIR}/${stem}.tbl.bak";
			open my $outtbl, '>', "${INDIR}/${stem}.tbl" or die $!;
			
			# initialize variables
			my $tr = Bio::WGS2NCBI::TableReader->new( 
				'-file' => "${INDIR}/${stem}.tbl.bak",
				'-cb'   => sub {
					my $id = shift;
					print $outtbl '>Features ', $id, "\n";
				}
			);
						
			# iterate over features
			my ( $oldid, $drop, $id ) = ( '' );
			while( my $f = $tr->next_feature ) {
				$id = $tr->seq;
				if ( $f->isa('Bio::WGS2NCBI::GeneFeature') ) {
					if ( $f->lies_within( @{ $coord{$id} } ) ) {
						$drop = 0;
					}
					else {
						$drop = 1;
					}
				}				
				if ( not $drop ) {
					
					# shift features leftward		
					if ( my $diff = $coord{$id}->[0] ) {
						my @r = $f->range;
						for my $r ( @r ) {
							my @coord;
							for my $coord ( @$r ) {
								if ( $coord =~ /^([^0-9]*)(\d+)$/ ) {
									my $prefix = $1;
									my $number = $2;
									$number -= $diff;
									push @coord, $prefix . $number;
								}
							}
							$r->[0] = $coord[0];
							$r->[1] = $coord[1];
						}
					}
					print $outtbl $f->to_string;
				}
				$oldid = $id;
			}
		}
	}
}

=head1 prune

The C<prune> action reads a discrepancy file as supplied by NCBI, parses out errors that
have locations in them, which are then pruned from the table files in $config->datadir.

This requires the following configuration settings:

=over

=item C<datadir>

The location of the dir where the (potentially 'chunked', see below) sequence files
and feature tables were written by L<Bio::WGS2NCBI/process>.

=item C<validation>

The location where to read the validation report from NCBI.

=item C<prefix>

The ID prefix that was assigned to you by NCBI when you created your submission, something 
like 'CR513_'

=item C<authority>

The naming authority prefix that you chose for your identifiers, something like 
'gnl|aceprd|'

The

=back

=cut

sub prune {
	my $config  = Bio::WGS2NCBI::Config->new;
	my $datadir = $config->datadir;
	my $discrep = $config->validation;
	my $auth    = $config->authority; # gnl|aceprd|
	my $prefix  = $config->prefix;    # CR513_
	my %locations;
	my %objects;
	
	# read locations from coordinates in report
	open my $fh, '<', $discrep or die $!;
	while(<$fh>) {
	
		# [(lcl|contig_3124:c<739-1821, 2000-2254)]
		my $coord = qr/c?[<>]?\d+/;
		my $seq = qr/[^:]+/;
		if ( /lcl\|($seq):($coord-$coord)((?:, $coord-coord)*)/ ) {
			my ( $seq, $first, $remainder ) = ( $1, $2, $3 );
			$locations{$seq} = [] if not $locations{$seq};
			push @{ $locations{$seq} }, [ split /-/, $first ];
			INFO "Going to prune annotations on $seq:$first";
			
			for my $r ( split /, /, $remainder ) {
				if ( $r =~ /($coord)-($coord)/ ) {
					my ( $start, $stop ) = ( $1, $2 );
					push @{ $locations{$seq} }, [ $start => $stop ];
					INFO "Extending range to prune $start => $stop";
				}
			}
		}
		if ( / \Q$auth\E($prefix\d+)/ ) {
			my $id = $1;
			INFO "Going to prune annotations identified by $id";
			$objects{$id}++;
		}
	}
	
	# start reading feature tables
	opendir my $dh, $datadir or die $!;
	while( my $entry = readdir $dh ) {
		if ( $entry =~ /\.tbl$/ ) {
		
			# make pruned copy
			INFO "Going to start writing to $datadir/$entry.pruned";
			open my $fh, '>', "$datadir/$entry.pruned" or die $!;
		
			# instantiate reader and start scanning
			INFO "Going to start reading $datadir/$entry";
			my $read = Bio::WGS2NCBI::TableReader->new(  
				'-file' => "$datadir/$entry",
				'-cb'   => sub { print $fh ">Features @_\n" }
			);
			
			# iterate over features, write out with sequence separators
			my $seq = '';
			FEAT: while( my $feat = $read->next_feature ) {
				
				# scan all the ranges
				$seq = $read->seq;
				if ( $locations{$seq} ) {
					for my $loc ( @{ $locations{$seq} } ) {
						if ( $feat->lies_within(@$loc) ) {
							WARN "Pruning from location $seq:".$loc->[0].'-'.$loc->[1];
							next FEAT;
						}
					}				
				}
				
				# check object ids
				if ( $feat->can('locus_tag') and $objects{$feat->locus_tag} ) {
					WARN "Pruning $feat ".$feat->locus_tag;
					next FEAT;
				}
				if ( $feat->can('protein_id') and $objects{$auth . $feat->protein_id} ) {
					WARN "Pruning $feat ".$feat->protein_id;
					next FEAT;
				}
				
				# still here? then print the feature
				print $fh $feat->to_string;
			}
		}	
	}	
}

=head1 compress

The C<compress> action bundles the ASN.1 files produced by C<Bio::WGS2NCBI/convert> into
a .tar.gz archive that can be uploaded to NCBI. This requires the following configuration 
settings:

=over

=item C<outdir>

The location where the ASN.1 files were written.

=item C<archive>

The name and location of the archive to produce.

=back

=cut

sub compress {
	my $config  = Bio::WGS2NCBI::Config->new;
	my $tar     = Archive::Tar->new;
	my $sqndir  = $config->outdir;
	my $archive = $config->archive;
	
	# open directory handle
	INFO "going to add *.sqn files from $sqndir to $archive";
	opendir my $dh, $sqndir or die $!;
	
	# iterate over files
	my $counter;
	while( my $entry = readdir $dh ) {
		if ( $entry =~ /\.sqn$/ ) {
			$counter++;
			$tar->add_files( "${sqndir}/${entry}" );
		}	
	}
	INFO "added $counter file(s)";
	
	# compress results
	INFO "going to compress $archive";
	$tar->write( $archive, COMPRESS_GZIP );
}

=head1 help

Displays module documentation (which you are reading now).

=cut

sub help {
	pod2usage({
		'-verbose'   => 2,
		'-exitval'   => 0,
		'-noperldoc' => 1,
	});
}

1;

