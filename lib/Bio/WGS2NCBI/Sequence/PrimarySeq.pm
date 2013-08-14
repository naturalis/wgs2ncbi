package Bio::WGS2NCBI::Sequence::PrimarySeq;

sub new {
	my $class = shift;
	my %args  = @_;
	my $self  = {};
	for my $key ( keys %args ) {
		my $field = $key;
		$field =~ s/^-/_/;
		$self->{$field} = $args{$key};
	}
	return bless $self, $class;
}

sub seq { ${ shift->{'_seq'} } }

sub subseq {
	my ($self, $start, $end) = @_;
	return substr ${ shift->{'_seq'} }, $start - 1, ( $end - $start ) + 1;
}

sub display_id {
	my $self = shift;
	$self->{'_id'} = shift if @_;
	return $self->{'_id'};
}

sub trunc {
	my ($self, $start, $end) = @_;
	return __PACKAGE__->new( 
		'-id'  => $self->display_id, 
		'-seq' => \( my $sub = $self->subseq($start,$end) ),
	);
}


1;
