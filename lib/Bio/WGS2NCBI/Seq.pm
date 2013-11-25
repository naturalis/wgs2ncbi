package Bio::WGS2NCBI::Seq;
use strict;
use warnings;

# this sequence class implements the needed bits of the Bio::Seq interface while holding
# the raw sequence data as a scalar ref

sub new {
	my $class = shift;
	my %args = @_;
	my $self = bless {}, $class;
	for ( qw(id seq desc) ) {
		$self->$_( $args{"-$_"} ) if $args{"-$_"};
	}
	return $self;
}

sub seq {
	my $self = shift;
	if ( @_ ) {
		my $seq = shift;
		if ( ref $seq ) {
			$self->{'_seq'} = $seq;
		}
		else {
			$self->{'_seq'} = \$seq;
		}
	}
	if ( ref $self->{'_seq'} ) {
		return ${ $self->{'_seq'} };
	}
	else {
		return $self->{'_seq'};
	}
}

sub length {
	my $self = shift;
	return $self->{'_seq'} ? length ${ $self->{'_seq'} } : 0;
}

sub id {
	my $self = shift;
	if ( @_ ) {
		$self->{'_id'} = shift;
	}
	return $self->{'_id'};
}

sub desc {
	my $self = shift;
	if ( @_ ) {
		$self->{'_desc'} = shift;
	}
	return $self->{'_desc'};
}

sub subseq {
	my ( $self, $start, $stop ) = @_;
	my $subseq;
	if ( $self->{'_seq'} ) {
		$subseq = substr ${ $self->{'_seq'} }, $start - 1, $stop - $start + 1;
	}
	return $subseq;
}

sub trunc {
	my ( $self, $start, $stop ) = @_;
	my $substring;
	if ( $self->{'_seq'} ) {
		$substring = substr ${ $self->{'_seq'} }, $start - 1, $stop - $start + 1;
	}
	return ref($self)->new( 
		'-id'   => $self->id,
		'-desc' => $self->desc,
		'-seq'  => \$substring,
	);
}

sub revcom {
	my $self = shift;
	my $revcom;
	if ( $self->{'_seq'} ) {
		$revcom = reverse ${ $self->{'_seq'} };
		$revcom =~ tr/ACGTacgt/TGCAtgca/;
	}
	return ref($self)->new( 
		'-id'   => $self->id,
		'-desc' => $self->desc,
		'-seq'  => \$revcom,
	);	
}

1;