use Config;
my $perl   = $Config{'perlpath'};
my $script = shift @ARGV;
open my $fh, '>', $script;
print $fh "#!$perl\n";
print $fh <<'LITERAL';
use Bio::WGS2NCBI;
my $action = shift;
if ( $action and Bio::WGS2NCBI->can($action) ) {
    Bio::WGS2NCBI->$action;
    exit 0;
}
Bio::WGS2NCBI->help;
LITERAL
close $fh;
chmod 0755, $script;
