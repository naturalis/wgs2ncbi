package Bio::WGS2NCBI::StrandedFeature;
use Bio::WGS2NCBI::Feature;
use base 'Bio::WGS2NCBI::Feature';

sub strand { shift->{'strand'} }

1;