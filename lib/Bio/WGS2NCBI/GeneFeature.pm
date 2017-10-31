package Bio::WGS2NCBI::GeneFeature;
use Bio::WGS2NCBI::StrandedFeature;
use Bio::WGS2NCBI::Logger;
use Data::Dumper;
use base 'Bio::WGS2NCBI::StrandedFeature';

sub new {
    my $class = shift;
    my %args = @_;
    DEBUG "need locus_tag" if not $args{'locus_tag'};
    $class->SUPER::new( 'five_prime_UTR' => [], 'three_prime_UTR' => [], %args );
}

sub db_xref { 
	my $self = shift;
	$self->{'db_xref'} = shift if @_;
	return $self->{'db_xref'};
}

sub gene { 
	my $self = shift;
	$self->{'gene'} = shift if @_;
	return $self->{'gene'};
}

sub product {
	my $self = shift;	
	$self->{'product'} = shift if @_;
	return $self->{'product'};
}

sub locus_tag { 
	my $self = shift;
	$self->{'locus_tag'} = shift if @_; 
	return $self->{'locus_tag'};	
}

sub note {
	my $self = shift;
	$self->{'note'} = shift if @_;
	return $self->{'note'};
}

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

sub qualifiers { qw(locus_tag gene note) }

sub annotate_short_introns {
    my ( $gene, $cds, $mrna, $set, $config ) = @_;
    DEBUG "Checking for short introns in ".$gene->product;   
    if ( $cds ) {  
		my @cds_ranges = $cds->range;
		if ( @cds_ranges > 1 ) {
			for my $i ( 1 .. $#cds_ranges ) {
				my $intron_length = abs($cds_ranges[$i]->[0] - $cds_ranges[$i-1]->[1]);
				if ( $intron_length <= $config->minintron ) {
					$set->remove($cds,$mrna);
					$gene->note('nonfunctional due to frameshift');
					INFO "Found $intron_length nt intron in '".$gene->product."', marked pseudo-gene";
					return 0;
				}
			}
		}
		else {
			DEBUG "Gene had only one exon, no intron lengths to validate in ".$gene->product;
		}
		$gene->note('');
    }
    else {
    	WARN "No CDS for ".$gene->{'type'}." ".$gene->locus_tag;
    }
	return 1;
}

sub annotate_partials {
    my ( $gene, $cds, $mrna, $seq ) = @_;
    if ( $cds ) {
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
    else {
    	WARN "No CDS for ".$gene->{'type'}." ".$gene->locus_tag;
    }
}

1;