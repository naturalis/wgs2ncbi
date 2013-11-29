use strict;
use warnings;
use Test::More 'no_plan';

BEGIN { 
	use Bio::WGS2NCBI::Config;
	push @ARGV, '-verbosity' => 2;
	my $config = Bio::WGS2NCBI::Config->new;
	ok( $config->verbosity == 2, 'incremented verbosity' );
}

BEGIN {
	use_ok('Bio::WGS2NCBI::Logger');
}

ok( INFO "*** info msg here is normal ***" == 2 );