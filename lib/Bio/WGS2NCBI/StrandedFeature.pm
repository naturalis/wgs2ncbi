package Bio::WGS2NCBI::StrandedFeature;
use Bio::WGS2NCBI::Feature;
use base 'Bio::WGS2NCBI::Feature';

=head1 NAME

Bio::WGS2NCBI::StrandedFeature - base class for features on an sequence strand

=head1 DESCRIPTION

Sequence features having to do with genes (i.e. CDS, mRNA, gene) are on one of the
two strands. Objects that inherit from this base class have a getter that identifies
this strand based on whatever is passed into their constructors, i.e. strand => '-' 
or strand => '+'.

=head1 SEE ALSO

L<Bio::WGS2NCBI::Feature>

=head1 METHODS

=over

=item strand()

Returns the strand, i.e. '+' or '-'.

=cut

sub strand { shift->{'strand'} }

=back

=cut

1;