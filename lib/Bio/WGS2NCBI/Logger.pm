package Bio::WGS2NCBI::Logger;
require Exporter;
use base 'Exporter';

our $Verbosity = 1;
our @EXPORT = qw(DEBUG INFO WARN);

sub LOG ($$) {
    my ($msg,$method) = @_;
    my ( $package, $file1up, $line1up, $sub ) = caller( 2 );
    my ( $pack0up, $file, $line, $sub0up )    = caller( 1 );
    my $log = sprintf( "%s %s [%s %s] - %s\n", uc $method, $sub || '', $0, $line, $msg );
    print STDERR $log;
    return $Verbosity;
}
sub DEBUG ($) { LOG shift, 'DEBUG' if $Verbosity >= 3 }
sub INFO ($)  { LOG shift, 'INFO'  if $Verbosity >= 2 }
sub WARN ($)  { LOG shift, 'WARN'  if $Verbosity >= 1 }

1;