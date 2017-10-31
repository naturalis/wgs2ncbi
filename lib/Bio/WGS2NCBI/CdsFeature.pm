package Bio::WGS2NCBI::CdsFeature;
use Bio::WGS2NCBI::MrnaFeature;
use base 'Bio::WGS2NCBI::MrnaFeature';

=head1 NAME

Bio::WGS2NCBI::CdsFeature - coding sequence feature

=head1 DESCRIPTION

An object of this class represents all CDS features
of a single gene.

=head1 SEE ALSO

L<Bio::WGS2NCBI::MrnaFeature>

=head1 METHODS

=over

=item qualifiers()

Returns the qualifiers to include in the feature table for a CDS. These
include all qualifiers used in a L<Bio::WGS2NCBI::MrnaFeature>, as well as
'note', 'codon_start', and 'db_xref'.

=cut

sub qualifiers { shift->SUPER::qualifiers, qw(note codon_start db_xref) }

sub note { 
	my $self = shift;
	$self->{'note'} = shift if @_;
	return $self->{'note'};
}

sub codon_start { 
	my $self = shift;
	$self->{'codon_start'} = shift if @_;
	return $self->{'codon_start'};
}

sub db_xref { 
	my $self = shift;
	$self->{'db_xref'} = shift if @_;
	return $self->{'db_xref'};
}

=item phase()

Returns the GFF3 'phase', i.e. the number of nucleotides after which the
next reading frame begins (0, 1 or 2).

=cut

sub phase { shift->{'phase'} }

=back

=cut

1;