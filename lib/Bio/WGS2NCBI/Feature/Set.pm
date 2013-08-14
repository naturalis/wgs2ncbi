package Bio::WGS2NCBI::Feature::Set;

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

sub offset {
    my $self = shift;
    $self->{'_offset'} = shift if @_;
    return $self->{'_offset'};
}

sub add {
    my ( $self, $feat ) = @_;
    push @{ $self->{'_features'} }, $feat;
}

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

1;