#!/usr/bin/perl
use strict;
use warnings;
use URI::Escape;
use Getopt::Long;
use File::Path 'make_path';

my $Verbosity = 1;

# GFF3 column numbers
my $chr_idx    = 0;
my $source_idx = 1;
my $type_idx   = 2;
my $start_idx  = 3;
my $end_idx    = 4;
my $strand_idx = 6;
my $codon_idx  = 7;
my $meta_idx   = 8;

sub LOG ($$) {
    my ($msg,$method) = @_;
    my ( $package, $file1up, $line1up, $sub ) = caller( 2 );
    my ( $pack0up, $file, $line, $sub0up )    = caller( 1 );
    my $log = sprintf( "%s %s [%s %s] - %s\n", uc $method, $sub || '', $0, $line, $msg );
    print STDERR $log;
}
sub DEBUG ($) { LOG shift, 'DEBUG' if $Verbosity >= 3 }
sub INFO ($)  { LOG shift, 'INFO'  if $Verbosity >= 2 }
sub WARN ($)  { LOG shift, 'WARN'  if $Verbosity >= 1 }

sub check_args {
    my ( $fasta, $gff3, $info, $dir, $source, $prefix, $authority, $limit, $chunksize );
    GetOptions(
        'verbose+'    => \$Verbosity,
        'dir=s'       => \$dir,
        'source=s'    => \$source,
        'prefix=s'    => \$prefix,
        'fasta=s'     => \$fasta,
        'gff3=s'      => \$gff3,
        'info=s'      => \$info,
        'authority=s' => \$authority,
        'limit=i'     => \$limit,
        'chunksize=i' => \$chunksize,
    );
    die "Need -fasta argument" if not $fasta or not -e $fasta;
    die "Need -gff3 argument" if not $gff3 or not -e $gff3;
    return Config->new(
        'fasta'     => $fasta, 
        'gff3'      => $gff3, 
        'info'      => $info, 
        'dir'       => $dir, 
        'source'    => $source, 
        'prefix'    => $prefix, 
        'authority' => $authority,
        'limit'     => $limit,
        'chunksize' => $chunksize,
    );
}

sub read_info {
    my $info = shift;
    
    # an info file was provided
    if ( -e $info ) {
        INFO "going to read metadata from $info";
        my %result;
        
        # open a file handle
        open my $fh, '<', $info or die $!;
        while(<$fh>) {
            chomp;
            
            # assume simple INI-style key/value pairs
            my ( $key, $value ) = split /=/, $_;
            $result{$key} = $value;
        }
        return %result;
    }
    else {
        INFO "no additional metadata provided";
    }
}

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
    my ( $id, $seq, $fh, $offset, %info ) = @_;
    
    # print basic header
    print $fh '>', $id;
    
    # add extra info
    if ( %info ) {
        print $fh ' ', join(' ', map { "[$_=$info{$_}]" } keys %info), "\n";
    }
    
    # line break
    for my $line ( unpack "(a80)*", substr $$seq, $offset ) {
        print $fh $line, "\n";
    }
    INFO "wrote sequence $id";
}

sub read_features {
    my ( $fh, $chr, $counter, $config, $seq ) = @_;
    INFO "reading features for $chr";
    
    # instantiate new set
    my $set = FeatureSet->new( 'seqid' => $chr );
    
    # this is re-set for every gene
    my ( $gene, $cds, $mrna );
    
    LINE: while(<$fh>) {
        chomp;
        my @line = split /\t/, $_;
        my ( $start, $end ) = $line[$strand_idx] eq '+' ? @line[$start_idx,$end_idx] : @line[$end_idx,$start_idx];
        
        # initialize shared constructor args
        my %args = ( 'type'  => $line[$type_idx] );
        
        # found a new gene, start creating objects
        if ( $line[$type_idx] eq 'gene' ) {
            finalize_gene( $gene, $cds, $mrna, $seq ) if $gene;
            
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
            $gene = Feature->new(%args);
            $set->add($gene);
        }
        elsif ( $line[$type_idx] eq 'CDS' ) {
        
            # we only create one CDS and one mRNA object, which both span multiple ranges
            if ( not $cds and not $mrna ) {
                INFO "instantiating new CDS and mRNA";
            
                # create the CDS
                $args{'product'}       = $gene->product;
                $args{'protein_id'}    = $config->authority . $gene->locus_tag;
                $args{'transcript_id'} = $config->authority . $gene->locus_tag . '.mrna';
                $cds = Feature->new(
                    %args,
                    'note'     => $gene->note,
                    'db_xref' => $gene->db_xref,
                    'range'   => [ [ $start, $end ] ],
                );
                $set->add($cds);    
                
                # create the mRNA       
                $args{'type'} = 'mRNA';
                $mrna = Feature->new( %args,'range' => [ [ $start, $end ] ]);
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
    finalize_gene( $gene, $cds, $mrna, $seq );    
    return $set;
}

sub finalize_gene {
    my ( $gene, $cds, $mrna, $seq ) = @_;
    INFO "finalizing gene ".$gene->product; 
    
    # the 5' and 3' untranslated regions need to have been observed
    # in order to be sure about the gene (and, transitively, CDS and mRNA) coordinates
    my @five_prime_UTR  = $gene->five_prime_UTR;
    my @three_prime_UTR = $gene->three_prime_UTR;
    my ($generange) = $gene->range;
    my @mrna_ranges = $mrna->range;
    my @cds_ranges  = $cds->range;  
        
    # if no 5' UTR was seen, all we know is that the gene, CDS and mRNA may start
    # before the coordinate we have seen
    if ( not @five_prime_UTR ) {
        $generange->[0] = '<' . $generange->[0];
        $mrna_ranges[0]->[0] = '<' . $mrna_ranges[0]->[0];
        $cds_ranges[0]->[0]  = '<' .  $cds_ranges[0]->[0];
    }
    else {
    
        # extend the mRNA
        if ( $gene->strand eq '+' ) {
            $mrna_ranges[0]->[0] = $five_prime_UTR[0]->[0];     
        }
        else {
            $mrna_ranges[0]->[0] = $five_prime_UTR[0]->[1];
        }
    }
    
    # if no 3' UTR was seen we don't know where the gene and mRNA ended, but we
    # do know for the CDSs (because of the stop codons) - unless there isn't one!
    if ( not @three_prime_UTR ) {
        $generange->[1] = '>' . $generange->[1];
        $mrna_ranges[-1]->[1] = '>' . $mrna_ranges[-1]->[1];  
        
		# maybe here we need to check whether there is a stop codon?
		my $strand = $gene->strand;
		if ( $strand eq '-' ) {
		
			# get integer coordinates without previously introduced <> symbols
			my $start = $cds_ranges[-1]->[0];	
			my $end   = $cds_ranges[-1]->[1];
			$start =~ s/[<>]//;			
			$end =~ s/[<>]//;
			
			# get the subsequence
			( $start, $end ) = sort { $a <=> $b } ( $start, $end );
			my $length = $end - $start;
			my $subseq = substr $$seq, $start - 1, $length;
			my $stop_codon  = uc substr $subseq, 0, 3;
			
			# reverse complement
			$stop_codon = reverse $stop_codon;
			$stop_codon =~ tr/ACGT/TGCA/;
			
			# should be a stop codon
			if ( $stop_codon ne 'TAG' and $stop_codon ne 'TAA' and $stop_codon ne 'TGA' ) {
				WARN "no 3' UTR and no stop codon in ".$gene->product;
				$cds_ranges[-1]->[1] = '>' . $cds_ranges[-1]->[1];
			}
		}
    }
    else {
    
        # extend the mRNA
        if ( $gene->strand eq '+' ) {
            $mrna_ranges[-1]->[1] = $three_prime_UTR[-1]->[1];
        }
        else {
            $mrna_ranges[-1]->[1] = $three_prime_UTR[-1]->[0];
        }
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
        if ( $args->{'note'} =~ /: (.+?) \(/ ) {
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
            	WARN "[GENE] product contains 'gene': ".$args->{'product'};
            }
            
            # product may contain database identifier
            if ( $args->{'product'} =~ /\d{3}/ ) {
            	WARN "[DBID] product may contain database identifier: ".$args->{'product'};
            }
            
            # product contains organelle
            if ( $args->{'product'} =~ /mitochondrial/i ) {
            	WARN "[ORGANELLE] product contains organelle: ".$args->{'product'};
            }
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

sub main {
    my $config = check_args();
    
    # read the info file, if any
    my %info = read_info( $config->info );
    
    # open the fasta file handle
    open my $fastaFH, '<', $config->fasta or die $!;
    my $fastaPos = 0;   
    
    # counter to generate autoincrementing IDs for gene, mrna and protein
    my $counter = {};
    my $seq_counter;
    
    # we will re-use these handles
    my ( $tblFH, $scaffoldFH );
    
    # iterate over the genome
    while( not eof $fastaFH ) {
        last if $config->limit and ++$seq_counter == $config->limit;
    
        # advance to the next scaffold/chromosome 
        my ( $offset, $chr, $seq ) = ( 0 );
        ( $fastaPos, $chr, $seq ) = read_fasta( $fastaFH, $fastaPos );
        $chr =~ s/\|/-/;
        
        # compute offset, if any
        if ( $$seq =~ /^(N+)/ ) {
        	my $leading_gap = $1;
        	$offset = length($leading_gap);
        	WARN "leading ${offset}bp gap in $chr, will strip this and apply offset";
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
            my $features = read_features( $gff3FH, $chr, $counter, $config, $seq );       
            write_features( $features, $tblFH );        
        }
        else {
            print $tblFH '>Features ', $chr, "\n";
        }
        
        # write the scaffold
        write_fasta( $chr, $seq, $scaffoldFH, $offset, %info );        
    }   
}
main();

package Config;

sub new {
    my $class = shift;
    my %args = @_;
    my $self = \%args;
    return bless $self, $class;
}

sub prefix { shift->{'prefix'} || 'OPHA_' }

sub source { shift->{'source'} || 'maker' }

sub dir { shift->{'dir'} || 'submission' }

sub fasta { shift->{'fasta'} or die "Need FASTA file to operate on!" }

sub gff3 { shift->{'gff3'} or die "Need GFF3 file to operate on!" }

sub authority { shift->{'authority'} || 'gnl|NaturalisBC|' }

sub chunksize { shift->{'chunksize'} || 1 }

sub info { shift->{'info'} }

sub limit { shift->{'limit'} }

package Feature;

sub new {
    my $class = shift;
    my %args = @_;
    if ( $class eq __PACKAGE__ ) {
        my $subclass = ucfirst( lc $args{'type'} ) . $class;
        return $subclass->new( 'range' => [], %args);
    }   
    my $self = \%args;
    return bless $self, $class;
}

sub range { 
    my $self = shift;
    push @{ $self->{'range'} }, shift if @_;
    return @{ $self->{'range'} };
}

sub qualifiers { }

sub to_string {
    my $self = shift;
    my @range = $self->range;
    my $result = $range[0]->[0] . "\t" . $range[0]->[1] . "\t" . $self->{'type'} . "\n";
    for my $i ( 1 .. $#range ) {
        $result .= $range[$i]->[0] . "\t" . $range[$i]->[1] . "\n";
    }
    for my $q ( $self->qualifiers ) {
        my @values;
        if ( ref $self->{$q} and ref $self->{$q} eq 'ARRAY' ) {
            @values = @{ $self->{$q} };
        }
        else {
            @values = ( $self->{$q} );
        }
        for my $v ( @values ) {
            $result .= ( "\t" x 3 ) . $q . "\t" . $v . "\n" if $v;
        }
    }
    return $result;
}

package GeneFeature;
use base 'Feature';

sub new {
    my $class = shift;
    my %args = @_;
    die "need locus_tag" if not $args{'locus_tag'};
    $class->SUPER::new( 'five_prime_UTR' => [], 'three_prime_UTR' => [], %args );
}

sub note { shift->{'note'} }

sub db_xref { shift->{'db_xref'} }

sub strand { shift->{'strand'} }

sub gene { shift->{'gene'} }

sub five_prime_UTR {
    my $self = shift;
    push @{ $self->{'five_prime_UTR'} }, shift if @_;
    return @{ $self->{'five_prime_UTR'} };
}

sub three_prime_UTR {
    my $self = shift;
    push @{ $self->{'three_prime_UTR'} }, shift if @_;
    return @{ $self->{'three_prime_UTR'} };
}

sub product { shift->{'product'} }

sub locus_tag { shift->{'locus_tag'} }

sub qualifiers { qw(locus_tag gene) }

package MrnaFeature;
use base 'Feature';

sub new {
    my $class = shift;
    my %args = @_;
    if ( not $args{'product'} or not $args{'protein_id'} or not $args{'transcript_id'} ) {
        die "need product, protein_id and transcript_id, product_id";
    }
    $class->SUPER::new(%args);  
}

sub qualifiers { qw(product protein_id transcript_id) }

package CdsFeature;
use base 'MrnaFeature';

sub qualifiers { shift->SUPER::qualifiers, qw(note codon_start db_xref) }

package FeatureSet;

sub new {
    my $class = shift;
    my %args = @_;
    my $self = {
        '_seqid'    => $args{'seqid'},
        '_features' => $args{'features'} || [],
        '_offset'   => $args{'offset'} || 0,
    };
    return bless $self, $class;
}

sub offset {
	my $self = shift;
	$self->{'_offset'} = shift if @_;
	return $self->{'_offset'};
}

sub add {
    my ( $self, $feat ) = @_;
    push @{ $self->{'_features'} }, $feat;
}

sub focal_gene {
    my $self = shift;
    for ( my $i = $#{ $self->{'_features'} }; $i >= 0; $i-- ) {
        my $feat = $self->{'_features'}->[$i];
        return $feat if $feat->isa('GeneFeature');
    }
    return undef;
}

sub focal_mrnas {
    my $self = shift;
    my @mrna;
    my $i = $#{ $self->{'_features'} };
    while( not $self->{'_features'}->[$i]->isa('GeneFeature') ) {
        my $feat = $self->{'_features'}->[$i];
        push @mrna, $feat if $feat->isa('MrnaFeature');     
        $i--;
    }
    return @mrna;
}

sub focal_cdss {
    my $self = shift;
    my @cdss;
    my $i = $#{ $self->{'_features'} };
    while( not $self->{'_features'}->[$i]->isa('GeneFeature') ) {
        my $feat = $self->{'_features'}->[$i];
        push @cdss, $feat if $feat->isa('CdsFeature');      
        $i--;
    }
    return @cdss;
}

sub to_string {
    my $self = shift;
    my $result = '>Features ' . $self->{'_seqid'} . "\n";
    if ( my $offset = $self->offset ) {
    	$result .= "[offset=-${offset}]\n";
    }
    for my $feat ( @{ $self->{'_features'} } ) {
        $result .= $feat->to_string;
    }
    return $result;
}