package Bio::WGS2NCBI::FeatureSet;
use Bio::WGS2NCBI::Logger;
use strict;
use warnings;
use URI::Escape;

# GFF3 column numbers
my $chr_idx    = 0;
my $source_idx = 1;
my $type_idx   = 2;
my $start_idx  = 3;
my $end_idx    = 4;
my $strand_idx = 6;
my $codon_idx  = 7;
my $meta_idx   = 8;

=head1 NAME

Bio::WGS2NCBI::FeatureSet - container of sets of features

=head1 DESCRIPTION

Objects from this class maintain an internal set of all 
sequence features for a given sequence ID. Features can
are added (or removed) from the set while L<Bio::WGS2NCBI>
is reading a sequence with its annotations. In addition,
objects from this class maintain an (integer) offset by
which all features are shifted, e.g. in the case where
a scaffold starts with an assembly gap that needs to be
removed. Lastly, objects from this class implement the
logic for serializing all the annotations for a given
sequence in feature table format. 

=head1 METHODS

=over

=item new()

Returns a new Bio::WGS2NCBI::FeatureSet object. This
constructor has one required, named argument: 'seqid',
which is the identifier of the sequence for which this
set holds the features. In addition, optionally,
an array reference of 'features' and an 'offset' can be 
provided. For example:

 my $set = Bio::WGS2NCBI::FeatureSet->new(
 	seqid    => 'scaffold1.1',
 	features => [ $feat1, $feat2 ],
 	offset   => 33,
 );

=cut

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

=item offset()

Gets/sets the offset by which all feature coordinates
must be shifted.

=cut

sub offset {
    my $self = shift;
    $self->{'_offset'} = shift if @_;
    return $self->{'_offset'};
}

=item add()

Adds one or more features to the set.

=cut

sub add {
    my ( $self, @feat ) = @_;
    push @{ $self->{'_features'} }, @feat;
}

=item remove()

Removes one or more features from the set.

=cut

sub remove {
	my ( $self, @feat ) = @_;
	my %remove = map { $_->id => 1 } @feat;
	$self->{'_features'} = [ grep { ! $remove{ $_->id } } @{ $self->{'_features'} } ];

}

=item to_string()

Writes the set out in feature table format.

=cut

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

=item read_features

Reads and returns a FeatureSet. Arguments:

- the file handle to read from
- a counter hash to generate IDs
- the Config object
- the focal Seq object
- the offset at the start of the sequence (e.g. after leading NNNs)
- the index of the last non-N residue in the sequence

=cut

sub read_features {
    my ( $class, $fh, $counter, $config, $seq, $offset, $last_non_missing_index ) = @_;
    INFO "Reading features for ".$seq->id;
    
    # instantiate new set
    my $set = $class->new( 'seqid' => $seq->id );
    
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
            	$gene->annotate_short_introns($cds,$mrna,$set,$config);
            	$gene->annotate_partials($cds,$mrna,$seq);
            }
            
            # clear the caches
            $gene = $cds = $mrna = undef;
            
            # add additional constructor args
            $args{'strand'}    = $line[$strand_idx];
            $args{'locus_tag'} = $config->prefix . sprintf "%05s", ++$counter->{'gene'};
            $args{'range'}     = [ [ $start, $end ] ];
            
            # read the line
            _read_gene_line(\%args,@line);
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
		my $auth = $config->authority;
		if ( not $auth ) {
			WARN "No naming authority specified!";
		}
            
                # create the CDS                
                $args{'product'}       = $gene->product;
                $args{'protein_id'}    = $auth . $gene->locus_tag;
                $args{'transcript_id'} = $auth . $gene->locus_tag . '.mrna';
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
		$gene->annotate_short_introns($cds,$mrna,$set,$config);
		$gene->annotate_partials($cds,$mrna,$seq);
	}
            
    return $set;
}

# not a method, a helper function that is called from read_features
sub _read_gene_line {
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

=item write_features

Writes the FeatureSet to the provided file handle.

=cut

sub write_features {
    my ( $features, $fh ) = @_;
    print $fh $features->to_string;
}

=back

=cut

1;
