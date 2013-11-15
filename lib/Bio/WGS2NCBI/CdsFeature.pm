package Bio::WGS2NCBI::CdsFeature;
use Bio::WGS2NCBI::MrnaFeature;
use base 'Bio::WGS2NCBI::MrnaFeature';

sub qualifiers { shift->SUPER::qualifiers, qw(note codon_start db_xref) }

sub phase { shift->{'phase'} }

1;