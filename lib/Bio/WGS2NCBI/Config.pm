package Bio::WGS2NCBI::Config;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Path 'make_path';
use Bio::WGS2NCBI::Logger;

my $SINGLETON;

sub _file {
	my ( $key, $value, $self ) = @_;
	if ( $key ) {
		if ( not $value or not -e $value ) {
			die "argument -$key needs an existing file, not '$value'";
		}
		$self->{$key} = $value;
	}
	return 's';
}

sub _int {
	my ( $key, $value, $self ) = @_;
	if ( $key ) {
		if ( not defined $value or $value !~ /^\d+/ ) {
			die "argument -$key needs an integer, not '$value'";
		}
		$self->{$key} = $value;
	}
	return 'i';
}

sub _dir {
	my ( $key, $value, $self ) = @_;
	if ( $key ) {
		if ( not $value ) {
			die "argument -$key needs a name for a directory, not '$value'";
		}
		elsif ( not -d $value ) {
			INFO "will create directory $value";
			make_path( $value );
		}
		else {
			WARN "directory '$value' already exists, contents may be overwritten";
		}
		$self->{$key} = $value;
	}
	return 's';
}

sub _string {
	my ( $key, $value, $self ) = @_;
	if ( $key ) {
		if ( not $value ) {
			die "argument -$key needs a string, not '$value'";
		}
		$self->{$key} = $value;
	}
	return 's';
}

my %fields = (
	'prefix'    => \&_string,
	'template'  => \&_file,
	'info'      => \&_file,
	'products'  => \&_file,
	'masks'     => \&_file,
	'authority' => \&_string,
	'datadir'   => \&_dir,
	'datafile'  => \&_file,
	'gff3dir'   => \&_dir,
	'gff3file'  => \&_file,
	'source'    => \&_string,
	'chunksize' => \&_int,
	'verbosity' => \&_int,
	'minlength' => \&_int,
	'limit'     => \&_int,
	'minintron' => \&_int,
	'outdir'    => \&_dir,
	'discrep'   => \&_string,
);

sub _verbosity {
	my $v = shift;
	if ( $v and $v =~ /^(?:1|2|3)$/ ) {
		$Bio::WGS2NCBI::Logger::Verbosity = $v;
	}
}

{
	no strict 'refs';
	for my $key ( keys %fields ) {
		*$key = sub { shift->{$key} };
	}
}

sub new {
    my $class = shift;
	if ( $SINGLETON ) {
		return $SINGLETON;
	}
	else {
	
		# the location of a config ini file like wgs2ncbi.ini can be defined in an
		# environment variable called WGS2NCBI. This file will be read first.
		my %config  = $class->read_ini($ENV{'WGS2NCBI'}) if $ENV{'WGS2NCBI'};
		$SINGLETON  = \%config;
		_verbosity($config{'verbosity'});
		
		# the location of a config ini file can also be provided on the command line
		# using the -conf argument. Whether this overrides other command line arguments
		# or vice versa depends on the order in which the command line arguments are
		# given, as they are processed from left to right.
		my %options = (
			'conf=s' => sub {
				%config = $class->read_ini(pop);
				$SINGLETON = \%config;
			},
			'verbosity=i' => sub {
				my $v = pop;
				$SINGLETON->{'verbosity'} = $v;
				_verbosity($v);
			}
		);
		
		# make the other command line arguments in Getopt::Long style
		for my $key ( keys %fields ) {
			my $sub = $fields{$key};
			my $type = $sub->();
			if ( not exists $options{"${key}=${type}"} ) {
				$options{"${key}=${type}"} = sub { $sub->( @_, $SINGLETON ) };
			}
		}	
		
		# process command line arguments	
		GetOptions(%options);
    	return bless $SINGLETON, $class;
    }
}

sub read_ini {
	my ( $self, $file ) = @_;
	my %result;
	if ( $file and -e $file ) {
		open my $fh, '<', $file or die $!;
		while(<$fh>) {
			chomp;
			s/;.*$//; # strip comments
			if ( /^(.+?)=(.+)$/ ) {
				my ( $key, $value ) = ( $1, $2 );
				if ( $result{$key} and not ref $result{$key} ) {
					$result{$key} = [ $result{$key}, $value ];
				}
				elsif ( $result{$key} and ref $result{$key} eq 'ARRAY' ) {
					push @{ $result{$key} }, $value;
				}
				else {
					$result{$key} = $value;
				}
			}
			if ( /\[.*\]/ ) {
				INFO "ini-style headings are ignored: $_";
			}
		}
	}
	return %result;
}

1;