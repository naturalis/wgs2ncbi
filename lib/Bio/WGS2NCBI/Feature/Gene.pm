package Bio::WGS2NCBI::Fature::Gene;
use strict;
use base 'Bio::WGS2NCBI::Feature';

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

1;