package Bio::WGS2NCBI::Feature::Cds;
use strict;
use base 'Bio::WGS2NCBI::Feature::Mrna';

sub qualifiers { shift->SUPER::qualifiers, qw(note codon_start db_xref) }

1;