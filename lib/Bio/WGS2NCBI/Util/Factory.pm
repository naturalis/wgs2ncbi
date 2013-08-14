package Bio::WGS2NCBI::Util::Factory;
use strict;
use warnings;
use Bio::WGS2NCBI::Util::Config;

our $AUTOLOAD;
my $self;

sub new {
	my $class = shift;
	if ( not $self ) {
		$self = {};
		for ( qw(seqio_class featio_class) ) {
			my $class = $Config{'internal'}->{$_};
			s/_class//;
			$self->{$_} = $class;
		}
		bless $self, $class;
	}
	return $self;
}

sub AUTOLOAD {
	my $self = shift;
	my $method = $AUTOLOAD;
	$method =~ s/.+://;
	if ( $method =~ /create_(.+)/ ) {
		my $class = $1;
		return $self->{$class}->new(@_);
	}	
}

1;