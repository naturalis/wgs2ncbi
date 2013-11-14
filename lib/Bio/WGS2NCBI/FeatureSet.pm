package Bio::WGS2NCBI::FeatureSet;

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

sub focal_gene {
    my $self = shift;
    for ( my $i = $#{ $self->{'_features'} }; $i >= 0; $i-- ) {
        my $feat = $self->{'_features'}->[$i];
        return $feat if $feat->isa('GeneFeature');
    }
    return undef;
}

sub focal_mrnas {
    my $self = shift;
    my @mrna;
    my $i = $#{ $self->{'_features'} };
    while( not $self->{'_features'}->[$i]->isa('GeneFeature') ) {
        my $feat = $self->{'_features'}->[$i];
        push @mrna, $feat if $feat->isa('MrnaFeature');     
        $i--;
    }
    return @mrna;
}

sub focal_cdss {
    my $self = shift;
    my @cdss;
    my $i = $#{ $self->{'_features'} };
    while( not $self->{'_features'}->[$i]->isa('GeneFeature') ) {
        my $feat = $self->{'_features'}->[$i];
        push @cdss, $feat if $feat->isa('CdsFeature');      
        $i--;
    }
    return @cdss;
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