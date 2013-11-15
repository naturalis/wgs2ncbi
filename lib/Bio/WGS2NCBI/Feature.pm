package Bio::WGS2NCBI::Feature;

my $IDCOUNTER = 1;

sub new {
    my $class = shift;
    my %args = @_;
    if ( $class eq __PACKAGE__ ) {
    	my $type = ucfirst lc $args{'type'};
    	$class =~ s/::Feature$/::${type}Feature/;
    	eval "require $class";
        return $class->new( 'range' => [], %args );
    }
    my $self = \%args;
    $self->{'id'} = $IDCOUNTER++;
    return bless $self, $class;
}

sub range { 
    my $self = shift;
    push @{ $self->{'range'} }, shift if @_;
    return @{ $self->{'range'} };
}

sub id { shift->{'id'} }

sub qualifiers { }

sub to_string {
    my $self = shift;
    my @range = $self->range;
    my $result = $range[0]->[0] . "\t" . $range[0]->[1] . "\t" . $self->{'type'} . "\n";
    for my $i ( 1 .. $#range ) {
        $result .= $range[$i]->[0] . "\t" . $range[$i]->[1] . "\n";
    }
    for my $q ( $self->qualifiers ) {
        my @values;
        if ( ref $self->{$q} and ref $self->{$q} eq 'ARRAY' ) {
            @values = @{ $self->{$q} };
        }
        else {
            @values = ( $self->{$q} );
        }
        for my $v ( @values ) {
            $result .= ( "\t" x 3 ) . $q . "\t" . $v . "\n" if $v;
        }
    }
    return $result;
}

1;