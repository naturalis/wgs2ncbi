package Bio::WGS2NCBI::Sequence::IO;
use strict;
use Bio::WGS2NCBI::Util::Logger;
use Bio::WGS2NCBI::Sequence::PrimarySeq;

sub new {
	my $class = shift;
	my %args = @_;
	my $self = {
		'fh'  => $args{'-fh'},
		'pos' => 0,
	};
	if ( $args{'-format'} and $args{'-format'} !~ /^fasta$/i ) {
		FATAL "Can only read fasta, not ".$args{'-format'};
	}
	if ( $args{'-file'} ) {
		open $self->{'fh'}, $args{'-file'} or FATAL "Can't open $args{-file}: $!";
	}
	return bless $self, $class;
}

sub next_seq {
	my $self = shift;
    
    # advance to the current position
    my $fh = $self->{'fh'};
    seek $fh, $self->{'pos'}, 0;
    
    # these will store what's on the def line and the subsequent seq data
    my ( $chr, $seq );
    LINE: while(<$fh>) {
        chomp;
        
        # we have a definition line, break if we've already processed a defline
        if ( />(.+)$/ ) {
            last LINE if $chr;
            $chr = $1;
        }
        else {
            $seq .= $_;
            
            # we store the current position here so that it is set to 
            # just before the defline when we encounter it
            $self->{'pos'} = tell $self->{fh};
        }       
    }   
    return Bio::WGS2NCBI::Sequence::PrimarySeq->new( '-id' => $chr, '-seq' => \$seq );
}

sub write_seq {
	my ( $self, $seq ) = @_;
	my $fh = $self->{'fh'};
    
    # print basic header
    print $fh '>', $seq->display_id, "\n";
    
    # fold lines at 80 characters
    for my $line ( unpack "(a80)*", $seq->seq ) {
        print $fh $line, "\n";
    }
	return 1;
}

1;
