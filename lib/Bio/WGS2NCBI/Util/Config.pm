package Bio::WGS2NCBI::Util::Config;
use strict;
use Getopt::Long;
use Bio::WGS2NCBI::Util::Logger;

our $AUTOLOAD;
our $ACTION;
my %Config;

sub import {
	my $ini;
	$ACTION = shift @ARGV;
	GetOptions(
		'config=s' => \$ini,
	);
	if ( -e $ini ) {
		read_ini($ini);
		validate();
	}
	else {
		FATAL "no config file!";
	}
}

sub read_ini {
	my $file = shift;
	open my $fh, '<', $file or FATAL "Can't open $file: $!";
	my $section = '_';
	$Config{$section} = {};
	bless $Config{$section}, __PACKAGE__;
	LINE: while(<$fh>) {
		chomp;
		next LINE if /^\s*$/; # skip blank lines
		next LINE if /^\s*;/; # skip comments
		if ( /^\s*\[(.+?)\]/ ) {
			$section = $1;
			$Config{$section} = {};
			bless $Config{$section}, __PACKAGE__;			
		}
		elsif ( /^([^=]+)=(.+)$/ ) {
			my ( $key, $value ) = ( $1, $2 );
			$Config{$section}->{$key} = $value;
		}
		else {
			ERROR "invalid INI line: $_";
		}
	}
}

sub open_handles {
	my ( $class, $section, $item, $mode ) = @_;
	$mode = '<' if not $mode;
	if ( my $hash = $Config{$section} ) {
		if ( my $files = $hash->{$item} ) {
			my @list = ref $files eq 'ARRAY' ? @{ $files } : $files;
			my @handles;
			for my $l ( @list ) {
				open my $fh, $mode, $l or FATAL "Can't open $l: $!";
				push @handles, $fh;				
			}
			return @handles;
		}
		else {
			FATAL "No item $item in section $section";
		}
	}
	else {
		FATAL "No such section: $section";
	}
}

sub validate {
	# check for verbosity setting and factories
	if ( my $internal = $Config{'internal'} ) {
		if ( exists $internal->{'verbosity'} ) {
			Bio::WGS2NCBI::Util::Logger->level($internal->{'verbosity'});
		}
		for my $key ( keys %{ $internal } ) {
			if ( $key =~ /_class$/ ) {
				my $class = $internal->{$key};
				eval "require $class";
				if ( $@ ) {
					FATAL "Can't load $key ($class): $@";
				}
			}
		}
	}
	else {
		FATAL "No internal settings supplied!";
	}

	# check for input files
	if ( my $project = $Config{'project'} ) {
		for my $key ( keys %{ $project } ) {
			if ( $key =~ /_file(s?)$/ ) {
				my $is_list = $1;
				my @values = split /,/, $project->{$key};
				my @clean;
				for my $file ( @values ) {
					if ( -e $file ) {
						INFO "$key: $file was found";
						push @clean, $file;
					}
					else {
						WARN "$key: $file was not found";
					}
				}
				$project->{$key} = $is_list ? \@clean : shift @clean;
			}
		}
	}
	else {
		FATAL "No project metadata supplied!";
	}
	
	# check for output directories
	if ( my $output = $Config{'output'} ) {
		for my $key ( keys %{ $output } ) {
			if ( $key =~ /_(?:dir|file)$/ ) {
				my $entry = $output->{$key};
				if ( -e $entry ) {
					WARN "$key ($entry) already exists!";
				}
			}
		}
	}
	else {
		FATAL "No output settings supplied!";
	}
}

sub keys {
	my ( $self, $key ) = @_;
	return ref $self ? keys %{ $self->{$key} } : keys %Config;
}

sub AUTOLOAD {
	my $self = shift;
	my $method = $AUTOLOAD;
	$method =~ s/.+://;
	
	# it's a section
	if ( ref $self ) {
		return $self->{$method};
	}
	elsif ( exists $Config{$method} ) {
		return $Config{$method};
	}
	else {
		DEBUG "ignoring $method";
	}
}

1;