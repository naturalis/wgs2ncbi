package Bio::WGS2NCBI::Config;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::WGS2NCBI;
use Bio::WGS2NCBI::Logger;
use File::Path 'make_path';

my $SINGLETON;

sub _file {
	my ( $key, $value, $self ) = @_;
	if ( $key ) {
		if ( not $value or not -e $value ) {
			ERROR "Argument -$key needs an existing file, not '$value'";
			exit(1);
		}
		$self->$key( $value );
	}
	return 's';
}

sub _int {
	my ( $key, $value, $self ) = @_;
	if ( $key ) {
		if ( not defined $value or $value !~ /^\d+/ ) {
			ERROR "Argument -$key needs an integer, not '$value'";
			exit(1);
		}
		$self->$key( $value );
	}
	return 'i';
}

sub _dir {
	my ( $key, $value, $self ) = @_;
	if ( $key ) {
		if ( not $value ) {
			ERROR "Argument -$key needs a name for a directory, not '$value'";
			exit(1);
		}
		elsif ( not -d $value ) {
			INFO "Will create directory $value";
			make_path( $value );
		}
		else {
			INFO "Directory '$value' already exists, contents (if any) may be overwritten";
		}
		$self->$key( $value );
	}
	return 's';
}

sub _string {
	my ( $key, $value, $self ) = @_;
	if ( $key ) {
		if ( not $value ) {
			ERROR "Argument -$key needs a string, not '$value'";
			exit(1);
		}
		$self->$key( $value );
	}
	return 's';
}

sub _array {
	my ( $key, $value, $self ) = @_;
	if ( $key ) {
		if ( not $value ) {
			ERROR "Argument -$key needs a string, not '$value'";
			exit(1);
		}
		$self->$key( $value );
	}
	return 's@';
}

my %fields = (
	'prefix'     => \&_string,
	'template'   => \&_file,
	'info'       => \&_file,
	'products'   => \&_file,
	'masks'      => \&_file,
	'authority'  => \&_string,
	'datadir'    => \&_dir,
	'datafile'   => \&_file,
	'gff3dir'    => \&_dir,
	'gff3file'   => \&_file,
	'validation' => \&_file,
	'source'     => \&_string,
	'chunksize'  => \&_int,
	'verbosity'  => \&_int,
	'complexity' => \&_int,
	'minlength'  => \&_int,
	'limit'      => \&_int,
	'minintron'  => \&_int,
	'outdir'     => \&_dir,
	'discrep'    => \&_string,
	'feature'    => \&_array,
	'tbl2asn'    => \&_string,
	'archive'    => \&_string,
);

sub verbosity {
	my ( $self, $v ) = @_;
	if ( $v and $v =~ /^[0-3]$/ ) {
		$Bio::WGS2NCBI::Logger::Verbosity = $v;
	}
	return $Bio::WGS2NCBI::Logger::Verbosity;
}

sub complexity {
	my ( $self, $c ) = @_;
	if ( defined $c ) {
		$Bio::WGS2NCBI::Logger::Complexity = $c;
	}
	return $Bio::WGS2NCBI::Logger::Complexity;
}

{
	# populates the symbol table with a method to get (and 
	# optionally set) the value for a field in %fields
	no strict 'refs';
	for my $key ( keys %fields ) {
		next if $key eq 'verbosity' or $key eq 'complexity';
		my $type = $fields{$key}->();
		
		# method for configuration fields that are arrays
		if ( $type =~ /\@/ ) {
			*$key = sub {
				my $self = shift;
				$self->{$key} = [] if not $self->{$key};
				if ( @_ ) {										
					push @{ $self->{$key} }, map { ref $_ eq 'ARRAY' ? @$_ : $_ } @_;
					$self->{'_configured'}++;
				}
				return @{ $self->{$key} };
			};		
		}
		
		# method for configuration fields that are scalar
		else {
			*$key = sub {
				my $self = shift;
				if ( @_ ) {
					$self->{$key} = shift;
					$self->{'_configured'}++;
				}
				return $self->{$key};
			};
		}
	}
}

sub new {
    my $class = shift;
	if ( $SINGLETON ) {
		return $SINGLETON;
	}
	else {
		$SINGLETON = { '_configured' => 0 };
		bless $SINGLETON, $class;
	
		# the location of a config ini file like wgs2ncbi.ini can be defined in an
		# environment variable called WGS2NCBI. This file will be read first.
		if ( $ENV{'WGS2NCBI'} and -e $ENV{'WGS2NCBI'} ) {			
			$SINGLETON->_selfconfig($ENV{'WGS2NCBI'});
			INFO 'Read configuration INI file from $WGS2NCBI=' . $ENV{'WGS2NCBI'};
		}
		else {
			if ( not @ARGV ) {
				WARN 'No configuration from $WGS2NCBI and no command line arguments';
			}
		}
		
		# the location of a config ini file can also be provided on the command line
		# using the -conf argument. Whether this overrides other command line arguments
		# or vice versa depends on the order in which the command line arguments are
		# given, as they are processed from left to right.
		my ( $verbosity, $config_ini ) = $Bio::WGS2NCBI::Logger::Verbosity;
		my %options = (
			'conf=s'   => \$config_ini,
			'verbose+' => \$verbosity,
			'help|?'   => sub { Bio::WGS2NCBI->help }
		);
		
		# make the other command line arguments in Getopt::Long style
		for my $key ( keys %fields ) {
			my $sub = $fields{$key};
			my $type = $sub->();
			if ( not exists $options{"${key}=${type}"} ) {
				$options{"${key}=${type}"} = sub { 
					$sub->( @_, $SINGLETON );
					$SINGLETON->{'_configured'}++;
				};
			}
		}	
		
		# process command line arguments	
		GetOptions(%options);
		
		# process any config file provided on the command line
		if ( $config_ini and -e $config_ini ) {
			$SINGLETON->_selfconfig($config_ini);
			INFO 'Read configuration INI file from -conf=' . $config_ini;
		}
		else {
			if ( not $SINGLETON->{'_configured'} ) {
				WARN 'No configuration INI file provided as -conf=<config file>';
			}
		}
		
		# process any verbosity setting provided on the command line
		if ( $verbosity != $Bio::WGS2NCBI::Logger::Verbosity ) {
			$SINGLETON->verbosity($verbosity);
		}
		
		# check if we are configured
		if ( not $SINGLETON->{'_configured'} ) {
			ERROR 'No configuration info provided anywhere, quitting.';
			ERROR "Try 'wgs2ncbi help' for more info.";
			exit(1);
		}
		
    	return bless $SINGLETON, $class;
    }
}

sub _selfconfig {
	my ( $self, $file ) = @_;
	my %config = $self->read_ini($file);
	for my $key ( keys %config ) {
		my $sub = $fields{$key};
		$sub->( $key, $config{$key}, $self );
		$self->{'_configured'}++;
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
				DEBUG "ini-style headings are ignored: $_";
			}
		}
	}
	return %result;
}

1;
