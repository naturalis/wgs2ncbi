package Bio::WGS2NCBI::FeatureSet;

=head1 NAME

Bio::WGS2NCBI::FeatureSet - container of sets of features

=head1 DESCRIPTION

Objects from this class maintain an internal set of all 
sequence features for a given sequence ID. Features can
are added (or removed) from the set while L<Bio::WGS2NCBI>
is reading a sequence with its annotations. In addition,
objects from this class maintain an (integer) offset by
which all features are shifted, e.g. in the case where
a scaffold starts with an assembly gap that needs to be
removed. Lastly, objects from this class implement the
logic for serializing all the annotations for a given
sequence in feature table format. 

=head1 METHODS

=over

=item new()

Returns a new Bio::WGS2NCBI::FeatureSet object. This
constructor has one required, named argument: 'seqid',
which is the identifier of the sequence for which this
set holds the features. In addition, optionally,
an array reference of 'features' and an 'offset' can be 
provided. For example:

 my $set = Bio::WGS2NCBI::FeatureSet->new(
 	seqid    => 'scaffold1.1',
 	features => [ $feat1, $feat2 ],
 	offset   => 33,
 );

=cut

sub new {
    my $class = shift;
    my %args = @_;
    my $self = {
        '_seqid'    => $args{'seqid'},
        '_features' => $args{'features'} || [],
        '_offset'   => $args{'offset'} || 0,
    };
    return bless $self, $class;
}

=item offset()

Gets/sets the offset by which all feature coordinates
must be shifted.

=cut

sub offset {
    my $self = shift;
    $self->{'_offset'} = shift if @_;
    return $self->{'_offset'};
}

=item add()

Adds one or more features to the set.

=cut

sub add {
    my ( $self, @feat ) = @_;
    push @{ $self->{'_features'} }, @feat;
}

=item remove()

Removes one or more features from the set.

=cut

sub remove {
	my ( $self, @feat ) = @_;
	my %remove = map { $_->id => 1 } @feat;
	$self->{'_features'} = [ grep { ! $remove{ $_->id } } @{ $self->{'_features'} } ];

}

=item to_string()

Writes the set out in feature table format.

=cut

sub to_string {
    my $self = shift;
    my $result = '>Features ' . $self->{'_seqid'} . "\n";
    if ( my $offset = $self->offset ) {
        $result .= "[offset=-${offset}]\n";
    }
    for my $feat ( @{ $self->{'_features'} } ) {
        $result .= $feat->to_string;
    }
    return $result;
}

=back

=cut

1;