package Bio::WGS2NCBI;
use strict;
use warnings;
use URI::Escape;
use Bio::WGS2NCBI::Logger;
use Bio::WGS2NCBI::Config;
use Bio::WGS2NCBI::Feature;
use Bio::WGS2NCBI::FeatureSet;
use File::Path 'make_path';

# version number for the whole release is defined here
our $VERSION = 1.0;

# export the run function
require Exporter;
use base 'Exporter';
our @EXPORT = qw(run);

# GFF3 column numbers
my $chr_idx    = 0;
my $source_idx = 1;
my $type_idx   = 2;
my $start_idx  = 3;
my $end_idx    = 4;
my $strand_idx = 6;
my $codon_idx  = 7;
my $meta_idx   = 8;

sub read_fasta {
    my ( $fh, $pos ) = @_;
    
    # advance to the current position
    seek $fh, $pos, 0;
    DEBUG "reading starts at $pos";
    
    # these will store what's on the def line and the subsequent seq data
    my ( $chr, $seq );
    LINE: while(<$fh>) {
        chomp;
        
        # we have a definition line, break if we've already processed a defline
        if ( />(.+)$/ ) {
            last LINE if $chr;
            $chr = $1;
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
    return $pos, $chr, \$seq;
}

sub write_fasta {
    my ( $id, $seq, $fh, $offset, $length, %info ) = @_;
    
    # print basic header
    print $fh '>', $id;
    
    # add extra info
    if ( %info ) {
        print $fh ' ', join(' ', map { "[$_=$info{$_}]" } keys %info), "\n";
    }
    
    # fold lines at 80 characters
    for my $line ( unpack "(a80)*", substr $$seq, $offset, $length ) {
        print $fh $line, "\n";
    }
    INFO "wrote sequence $id";
}

sub read_features {
    my ( $fh, $chr, $counter, $config, $seq, $offset, $last_non_missing_index ) = @_;
    INFO "reading features for $chr";
    
    # instantiate new set
    my $set = Bio::WGS2NCBI::FeatureSet->new( 'seqid' => $chr );
    
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
        		WARN "first gene on $chr outside of allowed range, skipping";
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
            INFO "reading gene " . $args{'product'};
            
            # create the object
            $gene = Bio::WGS2NCBI::Feature->new(%args);
            $set->add($gene);
        }
        elsif ( $line[$type_idx] eq 'CDS' ) {
        	next LINE if $skipgene;
        
            # we only create one CDS and one mRNA object, which both span multiple ranges
            if ( not $cds and not $mrna ) {
                INFO "instantiating new CDS and mRNA";
            
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
                INFO "growing range for CDS and mRNA ($start..$end)";
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
    INFO "checking for short introns in ".$gene->product;     
    my @cds_ranges = $cds->range;
    if ( @cds_ranges > 1 ) {
    	for my $i ( 1 .. $#cds_ranges ) {
    		my $intron_length = abs($cds_ranges[$i]->[0] - $cds_ranges[$i-1]->[1]);
    		if ( $intron_length <= $config->minintron ) {
    			$set->remove($cds,$mrna);
    			$gene->note('nonfunctional due to frameshift');
    			WARN "found $intron_length nt intron in ".$gene->product;
    			return 0;
    		}
    	}
    }
    else {
    	INFO "gene had only one exon, no intron lengths to validate in ".$gene->product;
    }
    $gene->note('');
	return 1;
}

sub annotate_partials {
    my ( $gene, $cds, $mrna, $seq ) = @_;
    INFO "finalizing gene ".$gene->product; 
    
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
    }
    else {
		$mrna_ranges[0]->[0] = $five_prime_UTR[0]->[0];     
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
    		$stop_codon = substr $$seq, $cds_end - 3, 3;
    	}
    	else {
    		$stop_codon = substr $$seq, $cds_end - 1, 3;
    		$stop_codon = reverse uc $stop_codon;
    		$stop_codon =~ tr/ACGT/TGCA/;
    	}        
        if ( $stop_codon !~ /^(?:TAG|TAA|TGA)$/ ) {
			WARN "no stop codon ($stop_codon) on $strand in ".$gene->product;
			$cds_ranges[-1]->[1] = '>' . $cds_ranges[-1]->[1];
		}
		else {
			INFO "stop codon $stop_codon on $strand in ".$gene->product;
		}
    }
    else {
		$mrna_ranges[-1]->[1] = $three_prime_UTR[-1]->[1];
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
        	$args->{'product'}  = $args->{'gene'};
        }
    }
    
    # there are potentially multiple database cross-references, e.g.
    # to Pfam and InterPro
    if ( $line[$meta_idx] =~ /Dbxref=([^;]+)/ ) {
        my $refs = $1;
        $refs =~ s/Pfam/PFAM/g;
        $args->{'db_xref'} = [ split /,/, $refs ];
    }   
}

sub write_features {
    my ( $features, $fh ) = @_;
    print $fh $features->to_string;
}

sub get_non_missing_index {
	my ( $seq, $reverse ) = @_;
	if ( $reverse ) {
		for ( my $i = length($$seq) - 1; $i >= 0; $i-- ) {
			return $i if substr( $$seq, $i, 1 ) ne 'N';
		}
	}
	else {
		for my $i ( 0 .. length($$seq) ) {
			return $i if substr( $$seq, $i, 1 ) ne 'N';
		}	
	}
}

sub run {
    my $config = Bio::WGS2NCBI::Config->new;
    
    # read the info file, if any
    my %info = $config->read_ini( $config->info );
    
    # open the fasta file handle
    open my $fastaFH, '<', $config->fasta or die $!;
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
        my ( $offset, $chr, $seq ) = ( 0 );
        ( $fastaPos, $chr, $seq ) = read_fasta( $fastaFH, $fastaPos );
        my $length = length $$seq;
        
        # compute offset, if any
        if ( $offset = get_non_missing_index($seq) ) {
            INFO "leading ${offset}bp gap in $chr, will strip this and apply offset";
        }
        
        # check what we have left
        my $last_non_missing_index = get_non_missing_index($seq,'reverse');
		$length -= ( $length - 1 - $last_non_missing_index ) + $offset;
		if ( $length < $config->minlength ) {
			WARN "remaining seq $chr is too short ($length bp), skipping";
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
            open $tblFH, '>', $config->dir . "/${stem}.tbl" or die $!;   
            open $scaffoldFH, '>', $config->dir . "/${stem}.fsa" or die $!;     
        }
        
        # get the features for that scaffold/chromosome, if we have them
        my $gff3 = $config->gff3 . "/${chr}.gff3";
        if ( -e $gff3 ) {
            open my $gff3FH, '<', $gff3 or die $!;  
            my $features = read_features( 
            	$gff3FH,  # file handle of focal file
            	$chr,     # scaffold/chromosome ID
            	$counter, # counter for generating locus tags
            	$config,  # config object
            	$seq,     # string reference of raw sequence
            	$offset,  # start offset
            	$last_non_missing_index, # last true seq character
            );  
            $features->offset($offset);     
            write_features( $features, $tblFH );        
        }
        else {
            print $tblFH '>Features ', $chr, "\n";
        }
        
        # write the scaffold
        write_fasta( $chr, $seq, $scaffoldFH, $offset, $length, %info );        
    }   
}

1;

