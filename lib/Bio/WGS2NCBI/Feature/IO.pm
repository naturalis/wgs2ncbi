package Bio::WGS2NCBI::Feature::IO;
use strict;
use Bio::WGS2NCBI::Util::Logger;

sub new {
	my $class = shift;
	my %args = @_;
	my $self = { 'fh'  => $args{'-fh'} };
	if ( $args{'-format'} and $args{'-format'} !~ /^gff3$/i ) {
		FATAL "can only read gff3, not ".$args{'-format'};
	}
	if ( $args{'-file'} ) {
		open $self->{'fh'}, $args{'-file'} or FATAL "Can't open $args{-file}: $!";
	}
	return bless $self, $class;
}

sub next_feature {

}

sub write_feature {
	my ( $self, $feat ) = @_;
}

1;
