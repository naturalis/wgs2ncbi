package Bio::WGS2NCBI::Seq;
use strict;
use warnings;
use Bio::WGS2NCBI::Logger;

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

sub get_non_missing_index {
	my ( $seq, $reverse ) = @_;
	
	# make a 1-based index array
	my @i = 1 .. $seq->length;
	
	# turn into decrementing if we read from the end
	@i = reverse @i if $reverse;
	
	# iterate over indices, return 0-based index of first non-missing residue
	for my $i ( @i ) {
		return $i - 1 if uc( $seq->subseq($i,$i) ) ne 'N';
	}
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

sub mask {
	my ( $seq, @masks ) = @_;
	my $raw = $seq->seq;
	my $id  = $seq->id;
	for my $m ( @masks ) {
		if ( $m =~ /(\d+)\.\.(\d+)/ ) {
			my ( $start, $stop ) = ( $1, $2 );
			my $index  = $start - 1;
			my $length = $stop - $index;
			substr $raw, $index, $length, ( 'N' x $length );
			WARN "Masked ${id}:${start}-${stop}";
		}		
	}
	$seq->seq( \$raw );
}

sub write_fasta {
    my ( $seq, $fh ) = @_;
    
    # print basic header
    print $fh '>', $seq->id;
    
    # add extra info
    if ( $seq->desc ) {
        print $fh ' ', $seq->desc, "\n";
    }
    
    # fold lines at 80 characters
    my $retval;
    for my $line ( unpack "(a80)*", $seq->seq ) {
        $retval = print $fh $line, "\n";
    }
    DEBUG "wrote sequence ".$seq->id;
    return $retval;
}

sub read_fasta {
    my ( $class, $fh, $pos ) = @_;
    $pos = 0 if not defined $pos;
    
    # advance to the current position
    seek $fh, $pos, 0;
    DEBUG "reading starts at $pos";
    
    # these will store what's on the def line and the subsequent seq data
    my ( $chr, $seq, $desc );
    LINE: while(<$fh>) {
        chomp;
        
        # we have a definition line, break if we've already processed a defline
        if ( /^>(.+)$/ ) {
            last LINE if $chr;
            my $defline = $1;
            ( $chr, $desc ) = split /\s+/, $defline, 2;
            DEBUG "going to concatenate sequence data for $chr";
        }
        else {
            $seq .= $_;
            
            # we store the current position here so that it is set to 
            # just before the defline when we encounter it
            $pos = tell $fh;
        }       
    }   
    DEBUG "reading ended at $pos";
    return $pos, $class->new( 
    	'-id'   => $chr, 
    	'-desc' => $desc,
    	'-seq'  => \$seq,     	
    );
}

1;