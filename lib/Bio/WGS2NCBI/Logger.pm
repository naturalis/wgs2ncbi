package Bio::WGS2NCBI::Logger;
use Term::ANSIColor;
require Exporter;
use base 'Exporter';

our $Complexity = 0;
our $Verbosity  = 0;
our @EXPORT      = qw(DEBUG INFO WARN ERROR);
our %COLORS      = (
	'DEBUG' => 'blue', 
	'INFO'  => 'green',
	'WARN'  => 'yellow',
	'ERROR' => 'red',
);

sub LOG ($$) {
    my ($msg,$method) = @_;

	# compose the message. Complex messages do some amount of introspection
	# to figure out where the message originated.
    my $log;
    if ( $Complexity ) {
		my ( $package, $file1up, $line1up, $sub ) = caller( 2 );
		my ( $pack0up, $file, $line, $sub0up )    = caller( 1 );
		$log = sprintf( "%s %s [%s %s] - %s\n", uc $method, $sub || '', $file, $line, $msg );
    }
    else {
    	$log = $method . ' ' . $msg . "\n";    	
    }
    
    # check if we are actually writing to the terminal rather than a file
    if ( -t STDERR ) {
    	print STDERR colored( $log, $COLORS{$method} );
    }
    else {
    	print STDERR $log;
    }
    return $Verbosity;
}
sub DEBUG ($) { LOG shift, 'DEBUG' if $Verbosity >= 3 }
sub INFO ($)  { LOG shift, 'INFO'  if $Verbosity >= 2 }
sub WARN ($)  { LOG shift, 'WARN'  if $Verbosity >= 1 }
sub ERROR ($) { LOG shift, 'ERROR' if $Verbosity >= 0 }

1;
