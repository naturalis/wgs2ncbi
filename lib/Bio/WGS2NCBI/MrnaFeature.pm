package Bio::WGS2NCBI::MrnaFeature;
use Bio::WGS2NCBI::Logger;
use Bio::WGS2NCBI::StrandedFeature;
use base 'Bio::WGS2NCBI::StrandedFeature';

=head1 NAME

Bio::WGS2NCBI::MrnaFeature - mRNA feature

=head1 DESCRIPTION

Objects of this class represent an mRNA feature for a gene.

=head1 SEE ALSO

L<Bio::WGS2NCBI::StrandedFeature>

=head1 METHODS

=over

=item new()

Returns a new mRNA feature. Requires the arguments 'product', 'protein_id' and 
'transcript_id', for example:

 my $mrna = Bio::WGS2NCBI::MrnaFeature->new(
 	product       => $product,
 	protein_id    => $protein_id,
 	transcript_id => $transcript_id,
 );

=cut

sub new {
    my $class = shift;
    my %args = @_;
    if ( not $args{'product'} or not $args{'protein_id'} or not $args{'transcript_id'} ) {
        DEBUG "need product, protein_id and transcript_id, product_id";
    }
    $class->SUPER::new(%args);  
}

=item qualifiers()

Returns the feature qualifiers for mRNA features, i.e. 'product', 'protein_id' and 
'transcript_id'

=cut

sub qualifiers { qw(product protein_id transcript_id) }

sub product { 
	my $self = shift;
	$self->{'product'} = shift if @_;
	return $self->{'product'};
}

sub protein_id { 
	my $self = shift;
	$self->{'protein_id'} = shift if @_;
	return $self->{'protein_id'};
}

sub transcript_id { 
	my $self = shift;
	$self->{'transcript_id'} = shift if @_;
	return $self->{'transcript_id'};
}

=back

=cut

1;