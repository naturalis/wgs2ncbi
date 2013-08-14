package Bio::WGS2NCBI::Util::Logger;
use strict;
use base 'Exporter';

our @EXPORT  = qw(DEBUG INFO WARN ERROR FATAL);
our $VERBOSE = 3;
my %level = (
	'debug' => 5,
	'info'  => 4,
	'warn'  => 3,
	'error' => 2,
	'fatal' => 1,
);

sub LOG ($$) {
    my ($msg,$method) = @_;
    my ( $package, $file1up, $line1up, $sub ) = caller( 2 );
    my ( $pack0up, $file, $line, $sub0up )    = caller( 1 );
    my $log = sprintf( "%s %s [%s %s] - %s\n", uc $method, $sub || '', $0, $line, $msg );
    print STDERR $log;
}
sub DEBUG ($) { LOG shift, 'DEBUG' if $VERBOSE >= 5 }
sub INFO ($)  { LOG shift, 'INFO'  if $VERBOSE >= 4 }
sub WARN ($)  { LOG shift, 'WARN'  if $VERBOSE >= 3 }
sub ERROR ($) { LOG shift, 'ERROR' if $VERBOSE >= 2 }
sub FATAL ($) { LOG shift, 'FATAL' if $VERBOSE >= 1; exit(1) }

sub level {
	my ( $class, $newlevel ) = @_;
	if ( $newlevel ) {
		if ( exists $level{lc $newlevel} ) {
			$VERBOSE = $level{lc $newlevel};
		}
		else {
			ERROR "Can't set verbosity to $newlevel";
		}
	}
	return $VERBOSE;
}

1;