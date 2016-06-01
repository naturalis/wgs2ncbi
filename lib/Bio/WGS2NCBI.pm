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
use File::Path 'make_path';

# version number for the whole release is defined here
our $VERSION = 'v0.1.0';

# export the run function
require Exporter;
use base 'Exporter';
our @EXPORT = qw(process prepare compress convert);

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
        ( $fastaPos, $seq ) = read_fasta( $fastaFH, $fastaPos );
        $seq->desc($desc);
        my $length = $seq->length;
        my $chr    = $seq->id;
                
        # compute offset at beginning, if any. e.g. if NNCGTNN, $offset is 2
        if ( $offset = get_non_missing_index($seq) ) {
            INFO "Leading ${offset}bp gap in $chr, will strip this and apply offset";
        }
        
        # check what we have left, here if NNCGTNN, $lnmi == 4
        my $last_non_missing_index = get_non_missing_index($seq,'reverse');
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
            my $features = read_features( 
            	$gff3FH,  # file handle of focal file
            	$counter, # counter for generating locus tags
            	$config,  # config object
            	$seq,     # sequence object
            	$offset,  # start offset
            	$last_non_missing_index, # last true seq character
            );  
            $features->offset($offset);     
            write_features( $features, $tblFH );        
        }
        else {
            print $tblFH '>Features ', $chr, "\n";
        }
        
		# truncate, mask, serialize
        my $trunc = $seq->trunc( $offset + 1, $offset + $length );
        mask_seq($trunc, @{$masks{$chr}}) if $masks{$chr};        
        write_fasta( $trunc, $scaffoldFH );        
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
	my $command = "$TBL2ASN -p $INDIR -t $TMPL -M n -a r10k -l paired-ends -r $OUTDIR -Z $DISCREP";
	INFO "going to execute command '$command'";
	exec $command;
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

sub help {
	pod2usage({-verbose => 2});
}

sub read_fasta {
    my ( $fh, $pos ) = @_;
    
    # advance to the current position
    seek $fh, $pos, 0;
    DEBUG "reading starts at $pos";
    
    # these will store what's on the def line and the subsequent seq data
    my ( $chr, $seq, $desc );
    LINE: while(<$fh>) {
        chomp;
        
        # we have a definition line, break if we've already processed a defline
        if ( /^>(.+)$/ ) {
            last LINE if $chr;
            my $defline = $1;
            ( $chr, $desc ) = split /\s+/, $defline, 2;
            DEBUG "going to concatenate sequence data for $chr";
        }
        else {
            $seq .= $_;
            
            # we store the current position here so that it is set to 
            # just before the defline when we encounter it
            $pos = tell $fh;
        }       
    }   
    DEBUG "reading ended at $pos";
    return $pos, Bio::WGS2NCBI::Seq->new( 
    	'-id'   => $chr, 
    	'-desc' => $desc,
    	'-seq'  => \$seq,     	
    );
}

sub write_fasta {
    my ( $seq, $fh ) = @_;
    
    # print basic header
    print $fh '>', $seq->id;
    
    # add extra info
    if ( $seq->desc ) {
        print $fh ' ', $seq->desc, "\n";
    }
    
    # fold lines at 80 characters
    for my $line ( unpack "(a80)*", $seq->seq ) {
        print $fh $line, "\n";
    }
    DEBUG "wrote sequence ".$seq->id;
}

sub read_features {
    my ( $fh, $counter, $config, $seq, $offset, $last_non_missing_index ) = @_;
    INFO "Reading features for ".$seq->id;
    
    # instantiate new set
    my $set = Bio::WGS2NCBI::FeatureSet->new( 'seqid' => $seq->id );
    
    # this is re-set for every gene
    my ( $gene, $cds, $mrna, $skipgene );
    
    LINE: while(<$fh>) {
        chomp;
        my @line = split /\t/, $_;
        
        # raw start and end coordinates when reading left-to-right. these need to be
        # reconciled with the strandedness, after reducing the left-to-right to the
        # last index for non-missing characters in the raw sequence. if we don't do
        # that we end up with features that run beyond the end of their sequence.
        my ( $ltrs, $ltre ) = @line[$start_idx,$end_idx];
        
        # last_non_missing_index is the 0-based index for the last character that is
        # not an N in the sequence. since feature tables use 1-based coordinates we
        # add one. No $ltre may exceed this adjusted end.
        my $adjusted_end = $last_non_missing_index + 1;
        $ltre = $adjusted_end if $ltre > $adjusted_end;
        my ( $start, $end ) = $line[$strand_idx] eq '+' ? ($ltrs,$ltre) : ($ltre,$ltrs);
        
        # initialize shared constructor args
        my %args = ( 'type'  => $line[$type_idx] );
        
        # found a new gene, start creating objects
        if ( $line[$type_idx] eq 'gene' ) {
        
        	# ignore genes that fall outside the range because we stripped
        	# leading gaps
        	if ( $line[$start_idx] - $offset <= 0 ) {
        		$skipgene = 1;
        		WARN "First gene on ".$seq->id()." outside of allowed range, skipping";
        		next LINE;
        	}
        	else {
        		$skipgene = 0;
        	}

			# finish the previous gene, if we have any
        	if ( $gene ) {
            	annotate_short_introns($gene,$cds,$mrna,$set,$config);
            	annotate_partials($gene,$cds,$mrna,$seq);
            }
            
            # clear the caches
            $gene = $cds = $mrna = undef;
            
            # add additional constructor args
            $args{'strand'}    = $line[$strand_idx];
            $args{'locus_tag'} = $config->prefix . sprintf "%05s", ++$counter->{'gene'};
            $args{'range'}     = [ [ $start, $end ] ];
            
            # read the line
            read_gene_line(\%args,@line);
            DEBUG "Gene: " . $args{'product'};
            
            # create the object
            $gene = Bio::WGS2NCBI::Feature->new(%args);
            $set->add($gene);
        }
        elsif ( $line[$type_idx] eq 'CDS' ) {
        	next LINE if $skipgene;
        
            # we only create one CDS and one mRNA object, which both span multiple ranges
            if ( not $cds and not $mrna ) {
                DEBUG "Instantiating new CDS and mRNA";
            
                # create the CDS                
                $args{'product'}       = $gene->product;
                $args{'protein_id'}    = $config->authority . $gene->locus_tag;
                $args{'transcript_id'} = $config->authority . $gene->locus_tag . '.mrna';
                $cds = Bio::WGS2NCBI::Feature->new(
                    %args,
                    'phase'   => $line[$codon_idx],
                    'note'    => $gene->note,
                    'db_xref' => $gene->db_xref,
                    'range'   => [ [ $start, $end ] ],
                );
                $set->add($cds);    
                
                # create the mRNA       
                $args{'type'} = 'mRNA';
                $mrna = Bio::WGS2NCBI::Feature->new( %args,'range' => [ [ $start, $end ] ]);
                $set->add($mrna);
            }
            else {
                DEBUG "growing range for CDS and mRNA ($start..$end)";
                $cds->range(  [ $start, $end ] );
                $mrna->range( [ $start, $end ] );
            }
        }
        elsif ( $line[$type_idx] eq 'five_prime_UTR' ) {
            $gene->five_prime_UTR( [ $line[$start_idx], $line[$end_idx] ] );
        }
        elsif ( $line[$type_idx] eq 'three_prime_UTR' ) {
            $gene->three_prime_UTR( [ $line[$start_idx], $line[$end_idx] ] );
        }
    }
    
	# finish the final gene, if we have any
	if ( $gene ) {
		annotate_short_introns($gene,$cds,$mrna,$set,$config);
		annotate_partials($gene,$cds,$mrna,$seq);
	}
            
    return $set;
}

sub annotate_short_introns {
    my ( $gene, $cds, $mrna, $set, $config ) = @_;
    DEBUG "Checking for short introns in ".$gene->product;     
    my @cds_ranges = $cds->range;
    if ( @cds_ranges > 1 ) {
    	for my $i ( 1 .. $#cds_ranges ) {
    		my $intron_length = abs($cds_ranges[$i]->[0] - $cds_ranges[$i-1]->[1]);
    		if ( $intron_length <= $config->minintron ) {
    			$set->remove($cds,$mrna);
    			$gene->note('nonfunctional due to frameshift');
    			WARN "Found $intron_length nt intron in '".$gene->product."', marked pseudo-gene";
    			return 0;
    		}
    	}
    }
    else {
    	DEBUG "Gene had only one exon, no intron lengths to validate in ".$gene->product;
    }
    $gene->note('');
	return 1;
}

sub annotate_partials {
    my ( $gene, $cds, $mrna, $seq ) = @_;
    DEBUG "Finalizing gene ".$gene->product; 
    
    # the 5' and 3' untranslated regions need to have been observed
    # in order to be sure about the gene (and, transitively, CDS and mRNA) coordinates
    my @five_prime_UTR  = $gene->five_prime_UTR;
    my @three_prime_UTR = $gene->three_prime_UTR;
    my ($generange) = $gene->range;
    my @mrna_ranges = $mrna->range;
    my @cds_ranges  = $cds->range;  
    my $strand      = $gene->strand;
        
    # if no 5' UTR was seen, all we know is that the gene, CDS and mRNA may start
    # before the coordinate we have seen. Prefix the start coordinate of the first
    # segment as 5' partial ('<'). If there is a UTR, the mRNA needs to be extended
    # to cover it.
    if ( not @five_prime_UTR ) {
		$generange->[0] = '<' . $generange->[0];
		$mrna_ranges[0]->[0] = '<' . $mrna_ranges[0]->[0];
		$cds_ranges[0]->[0]  = '<' .  $cds_ranges[0]->[0];
		DEBUG "No 5' UTR";
    }
    else {
		$mrna_ranges[0]->[0] = $five_prime_UTR[0]->[0];     
		DEBUG "Found 5' UTR";
    }
    
    # if no 3' UTR was seen we don't know where the gene and mRNA ended, but we
    # do know for the CDSs (because of the stop codons) - unless there isn't one
    if ( not @three_prime_UTR ) {
    	$generange->[1] = '>' . $generange->[1];
		$mrna_ranges[-1]->[1] = '>' . $mrna_ranges[-1]->[1];
    	
    	# parse the stop codon from the strand
    	my $stop_codon;
    	my $cds_end   = $cds_ranges[-1]->[1];    	
    	if ( $strand eq '+' ) {
    		$stop_codon = $seq->trunc( $cds_end - 2, $cds_end )->seq;
    	}
    	else {
    		$stop_codon = $seq->trunc( $cds_end, $cds_end + 2 )->revcom->seq;
    	}        
        if ( $stop_codon !~ /^(?:TAG|TAA|TGA)$/i ) {
			DEBUG "No stop codon ($stop_codon) on $strand in ".$gene->product;
			$cds_ranges[-1]->[1] = '>' . $cds_ranges[-1]->[1];
		}
		else {
			DEBUG "Stop codon $stop_codon on $strand in ".$gene->product;
		}
    }
    else {
		$mrna_ranges[-1]->[1] = $three_prime_UTR[-1]->[1];
		DEBUG "Found 3' UTR";
    }
}

sub read_gene_line {
    my ( $args, @line ) = @_;
    $args->{'product'} = 'hypothetical protein';
    
    # here we parse out the whole note and the gene symbol and product
    if ( $line[$meta_idx] =~ /Note=([^;]+)/ ) {
        $args->{'note'} = $1;
        
        # some GFF3 features are URL-encoded, this needs to be corrected
        if ( $args->{'note'} =~ /%/ ) {
            INFO "[URLENCODED] possible URL encoding in note: ".$args->{'note'};
            $args->{'note'} = uri_unescape($args->{'note'});
            INFO "[URLENCODED] corrected as: ".$args->{'note'};         
        }            
        if ( $args->{'note'} =~ /Similar to ([^:]+):/ ) {
            $args->{'gene'} = $1;
        }
        if ( $args->{'note'} =~ /: (.+?) \([^\(\)]+\)$/ ) {
            $args->{'product'} = $1;
            
            # product should not end with 'domain'
            if ( $args->{'product'} =~ /domain$/i ) {
                INFO "[DOMAIN] product ends with 'domain': ".$args->{'product'};
                $args->{'product'} .= ' protein';
                INFO "[DOMAIN] corrected as: ".$args->{'product'};
            }
            
            # Implies evolutionary relationship; change to -like protein
            if ( $args->{'product'} =~ /homolog/i ) {
                INFO "[HOMOLOG] product contains 'homolog': ".$args->{'product'};
                $args->{'product'} =~ s/ homolog/-like/i;
                INFO "[HOMOLOG] corrected as: ".$args->{'product'};
            }
            
            # features ends with 'like', 
            # Consider adding 'protein' to the end of the product name
            if ( $args->{'product'} =~ /like$/i ) {
                INFO "[LIKE] product ends with like: ".$args->{'product'};
                $args->{'product'} .= ' protein';
                INFO "[LIKE] corrected as: ".$args->{'product'};
            }
            
            # feature contains 'gene'
            if ( $args->{'product'} =~ /\bgene\b/i ) {
                INFO "[GENE] product contains 'gene': ".$args->{'product'};
            }
            
            # product may contain database identifier
            if ( $args->{'product'} =~ /\d{3}/ ) {
                INFO "[DBID] product may contain database identifier: ".$args->{'product'};
            }
            
            # product contains organelle
            if ( $args->{'product'} =~ /mitochondrial/i ) {
                INFO "[ORGANELLE] product contains organelle: ".$args->{'product'};
            }
            
            # product starts with 'Probable'
            if ( $args->{'product'} =~ /^Probable\b/i ) {
            	$args->{'product'} =~ s/^Probable/putative/i;
            	INFO "[PROBABLE] product starts with probable, changed to putative";
            }

            # product starts with 'Probable'
            if ( $args->{'product'} =~ /^Uncharacterized\b/i ) {
            	$args->{'product'} =~ s/^Uncharacterized/putative/i;
            	INFO "[UNCHARACTERIZED] product starts with Uncharacterized, changed to putative";
            }
        }
        if ( $args->{'gene'} and $args->{'product'} eq 'hypothetical protein' ) {
        	$args->{'product'} = $args->{'gene'};
        }
    }
    
    # there are potentially multiple database cross-references, e.g.
    # to Pfam and InterPro
    if ( $line[$meta_idx] =~ /Dbxref=([^;]+)/ ) {
        my $refs = $1;
        $refs =~ s/Pfam/PFAM/g;
        $args->{'db_xref'} = [ split /,/, $refs ];
    }
    
    # check if we need to replace from INI
	my $config = Bio::WGS2NCBI::Config->new;
	if ( my $file = $config->products ) {
		my %map = $config->read_ini( $file );
		if ( my $fixed = $map{ $args->{'product'} } ) {
			if ( ref $fixed and ref $fixed eq 'ARRAY' ) {
				$fixed = pop(@{$fixed});			
				WARN "Multiple mappings for the same product, will use last seen: '"
				. $args->{'product'}."' => '".$fixed."'";
			}
			INFO "Replacing '".$args->{'product'}."' with '$fixed' from '$file'";
			$args->{'product'} = $fixed; 
		}
	}
	
	# check to see if we should delete uninformative gene symbols
	if ( $args->{'gene'} and ( $args->{'gene'} eq 'hypothetical protein' ) ) {
		WARN "Gene symbol was set to 'hypothetical protein', deleting...";
		delete $args->{'gene'};
	}
	if ( $args->{'gene'} and ( $args->{'product'} eq 'hypothetical protein' ) ) {
		my $sym = $args->{'gene'};
		WARN "Product was 'hypothetical protein', moving gene name '$sym' to note";
		delete $args->{'gene'};	
		if ( not defined $args->{'note'} ) {	
			$args->{'note'} = $sym;
		}
		elsif ( not ref $args->{'note'} ) {
			$args->{'note'} = [ $args->{'note'}, $sym ];
		}
		else {
			push @{ $args->{'note'} }, $sym;
		}
	}
}

sub write_features {
    my ( $features, $fh ) = @_;
    print $fh $features->to_string;
}

sub get_non_missing_index {
	my ( $seq, $reverse ) = @_;
	
	# make a 1-based index array
	my @i = 1 .. $seq->length;
	
	# turn into decrementing if we read from the end
	@i = reverse @i if $reverse;
	
	# iterate over indices, return 0-based index of first non-missing residue
	for my $i ( @i ) {
		return $i - 1 if uc( $seq->subseq($i,$i) ) ne 'N';
	}
}

sub mask_seq {
	my ( $seq, @masks ) = @_;
	my $raw = $seq->seq;
	my $id  = $seq->id;
	for my $m ( @masks ) {
		if ( $m =~ /(\d+)\.\.(\d+)/ ) {
			my ( $start, $stop ) = ( $1, $2 );
			my $index = $start - 1;
			my $length = $stop - $index;
			substr $raw, $index, $length, ( 'N' x $length );
			WARN "Masked ${id}:${start}-${stop}";
		}		
	}
	$seq->seq( \$raw );
}

1;

