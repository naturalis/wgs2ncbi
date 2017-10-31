package Bio::WGS2NCBI::Feature;

my $IDCOUNTER = 1;

=head1 NAME

Bio::WGS2NCBI::Feature - base class for sequence features

=head1 DESCRIPTION

This is the base class for all sequence feature classes for which 
L<Bio::WGS2NCBI> instantiates objects when reading a genome annotation.
In particular, this class manages the coordinates (1-based inclusive
start/stop) for the region(s) that the feature spans, the internal unique
identifiers for all features, and the serialization logic.

The API for features is deliberately sparse as it is only intended to be
used by L<Bio::WGS2NCBI> and may change in future versions (e.g. in order 
to become compatible with BioPerl features so that other annotation file 
formats can be used as well).

=head1 METHODS

=over

=item new()

Returns a new feature. This constructor is called with the 'type' 
argument to create any of the subclasses. For example:

 my $gene = Bio::WGS2NCBI::Feature->new( type => 'gene' );

In addition, the values for all methods except 'qualifiers', 'id'
and 'to_string' can be set here by using the method name as the
key name of a named argument. For example:

 my $gene = Bio::WGS2NCBI::Feature->new( type => 'gene', gene => 'CO1' );

=cut

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


=item seq

Getter/setter of sequence id

=cut

sub seq {
	my $self = shift;
	$self->{'seq'} = shift if @_;
	return $self->{'seq'};
}

=item range()

Returns a list of ranges, specified as start/stop coordinates, 
that the feature spans. Optional argument ([ $start, $stop ]) 
grows the set of ranges by one additional range.

=cut

sub range { 
    my $self = shift;
    push @{ $self->{'range'} }, shift if @_;
    return @{ $self->{'range'} };
}

=item lies_within

Tests whether the invocant lies within the provided coordinates

=cut

sub lies_within {
	my ( $self, $start, $stop ) = @_;
	$start =~ s/^[c<>]//;
	$stop  =~ s/^[c<>]//;
	my @r = $self->range;
	my $r1 = $r[0]->[0];
	my $r2 = $r[-1]->[-1];
	$r1 =~ s/^[c<>]//;
	$r2 =~ s/^[c<>]//;
	if ( $start <= $stop ) {
		if ( $r1 >= $start and $r2 <= $stop ) {
			return 1;
		}
		return 0;
	}
	
	# opposite strand
	if ( $start >= $stop ) {
		if ( $r1 <= $start and $r2 >= $stop ) {
			return 1;
		}
		return 0;
	}
}

=item id()

Returns a unique integer ID for the feature.

=cut

sub id { shift->{'id'} }

=item qualifiers()

Returns a list of qualifier names. This is an empty ("abstract") method
that the subclasses need to override.

=cut

sub qualifiers { }

=item to_string()

Returns a string representation of the feature in NCBI's tbl format.

=cut

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

=back

=cut

1;