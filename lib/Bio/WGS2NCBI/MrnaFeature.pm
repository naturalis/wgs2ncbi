package Bio::WGS2NCBI::MrnaFeature;
use Bio::WGS2NCBI::StrandedFeature;
use base 'Bio::WGS2NCBI::StrandedFeature';

sub new {
    my $class = shift;
    my %args = @_;
    if ( not $args{'product'} or not $args{'protein_id'} or not $args{'transcript_id'} ) {
        die "need product, protein_id and transcript_id, product_id";
    }
    $class->SUPER::new(%args);  
}

sub qualifiers { qw(product protein_id transcript_id) }

1;