package Bio::WGS2NCBI::TableReader;
use strict;
use warnings;
use Bio::WGS2NCBI::Feature;
use Bio::WGS2NCBI::Logger;

sub new {
	my $class = shift;
	my %args = @_;
	my $self = { 
		'seq' => undef, 
		'fh'  => undef,
		'cb'  => sub { WARN "Starting @_ without a handler" }
	};
	
	# set callback
	if ( $args{'-cb'} ) {
		$self->{'cb'} = $args{'-cb'};
	}
	
	# open file
	if ( $args{'-file'} ) {
		INFO "Opening ".$args{'-file'};
		open $self->{'fh'}, '<', $args{'-file'} or die $!;
	}
	
	# store handle
	elsif ( $args{'-handle'} ) {
		INFO "Using provided handle ".$args{'-handle'};
		$self->{'fh'} = $args{'-handle'};
	}
	
	# need either -file or -handle
	else {
		ERROR "Need -file or -handle argument";
	}
	
	return bless $self, $class;
}

sub next_feature {
	my $self = shift;
	my $feat;
	FEAT: while( my $line = readline( $self->handle ) ) {
		chomp($line);

		# start of feature table (chunk)
		if ( $line =~ /^>Features (\S+)/ ) {
			my $seq = $1;
			
			# deal with >Features markers in two passes: first rewind, fall through
			# the loop, and return the (defined) $feat; second time store focal $seq
			if ( $feat ) {
				seek $self->handle, ( length($line) * -1 ) - length($/), 1;
				last FEAT;			
			}
			else {			
				$self->seq($seq);
				$self->callback->($seq);
				DEBUG "Starting features for sequence ".$self->seq;
			}
		}
			
		# should be five fields
		my @f = split /\t/, $line;
		DEBUG join "\t", @f;
		my ( $start, $stop, $type, $key, $value ) = @f;
		
		# start of new feature, don't have one yet
		if ( $type and defined $start and defined $stop and not $feat ) {
		
			# load appropriate class for feature type
			DEBUG "Instantiating new $type feature";
			$feat = Bio::WGS2NCBI::Feature->new( type => $type );
			
			DEBUG "Setting range $start => $stop";
			$feat->range([ $start => $stop ]);
		}
		
		# start of new feature, already have one
		elsif ( $type and defined $start and defined $stop and $feat ) {
			DEBUG "Finishing reading feature $feat";
			
			# seek back from current pos (3rd arg) to the start of the line,
			# taking into account the chomped line break
			seek $self->handle, ( length($line) * -1 ) - length($/), 1;
			last FEAT;
		}
		
		# line is a range extension
		elsif ( $start and $stop and not $type ) {
			DEBUG "Adding range $start => $stop";		
			$feat->range([ $start => $stop ]);
		}
	
		# line is an annotation
		elsif ( defined $key and defined $value ) {
			DEBUG "Setting '$key' to '$value'";
			$feat->$key($value);
		}		
		
		last FEAT if eof $self->handle;
	}
	return $feat;
}

sub callback {
	my ( $self, $cb ) = @_;
	if ( $cb ) {
		$self->{'cb'} = $cb;
	}
	return $self->{'cb'};
}

sub seq {
	my ( $self, $seq ) = @_;
	if ( $seq ) {
		$self->{'seq'} = $seq;
	}
	return $self->{'seq'};
}

sub handle {
	my ( $self, $handle ) = @_;
	if ( $handle ) {
		$self->{'fh'} = $handle;
	}
	return $self->{'fh'};
}


1;