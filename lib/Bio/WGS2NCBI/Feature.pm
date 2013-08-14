package Bio::WGS2NCBI::Feature;

sub new {
    my $class = shift;
    my %args = @_;
    if ( $class eq __PACKAGE__ ) {
        my $subclass = $class . '::' . ucfirst lc $args{'type'};
        return $subclass->new( 'range' => [], %args);
    }   
    my $self = \%args;
    return bless $self, $class;
}

sub range { 
    my $self = shift;
    push @{ $self->{'range'} }, shift if @_;
    return @{ $self->{'range'} };
}

sub primary_tag { shift->{'type'} }

sub get_all_tags { }

sub get_tag_values {
	my ( $self, $tag ) = @_;
	return ref $self->{$tag} eq 'ARRAY' ? @{ $self->{$tag} } : $self->{$tag};
}

sub get_tagset_values {
	my $self = shift;
	return map { $self->get_tag_values } @_;
}

sub has_tag {
	my ( $self, $tag ) = @_;
	return grep { $_ eq $tag } $self->get_all_tags;
}

sub to_FTstring {
    my $self = shift;
    my @range = $self->range;
    my $result = $range[0]->[0] . "\t" . $range[0]->[1] . "\t" . $self->primary_tag . "\n";
    for my $i ( 1 .. $#range ) {
        $result .= $range[$i]->[0] . "\t" . $range[$i]->[1] . "\n";
    }
    for my $tag ( $self->get_all_tags ) {
        for my $v ( $self->get_tag_values($tag) ) {
            $result .= ( "\t" x 3 ) . $q . "\t" . $v . "\n" if $v;
        }
    }
    return $result;
}

1;