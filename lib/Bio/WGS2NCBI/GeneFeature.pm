package Bio::WGS2NCBI::GeneFeature;
use Bio::WGS2NCBI::StrandedFeature;
use base 'Bio::WGS2NCBI::StrandedFeature';

sub new {
    my $class = shift;
    my %args = @_;
    die "need locus_tag" if not $args{'locus_tag'};
    $class->SUPER::new( 'five_prime_UTR' => [], 'three_prime_UTR' => [], %args );
}

sub db_xref { shift->{'db_xref'} }

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

sub note {
	my $self = shift;
	if ( @_ ) {
		$self->{'note'} = shift;
	}
	return $self->{'note'};
}

1;