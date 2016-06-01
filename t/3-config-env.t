use strict;
use warnings;
use Test::More 'no_plan';

BEGIN { 
	use_ok('Bio::WGS2NCBI::Config'),
	$ENV{'WGS2NCBI'} = 'share/wgs2ncbi.ini';
}

my $config = Bio::WGS2NCBI::Config->new;

ok( $config->prefix eq 'L345_', 'prefix' );
ok( $config->template eq 'share/template.sbt', 'template' );
ok( $config->info eq 'share/info.ini', 'info' );
ok( $config->products eq 'share/products.ini', 'products' );
ok( $config->masks eq 'share/adaptors.ini', 'masks' );
ok( $config->authority eq 'gnl|NaturalisBC|', 'authority' );
ok( $config->datadir eq 'share/tblfasta', 'datadir' );
ok( $config->datafile eq 'share/king_cobra_scaffolds_spring_2011-cleaned.fasta', 'datafile' );
ok( $config->gff3dir eq 'share/gff3', 'gff3dir' );
ok( $config->gff3file eq 'share/cobra.functional-cleaned.gff', 'gff3file' );
ok( $config->source eq 'maker', 'source' );
ok( $config->chunksize == 5000, 'chunksize' );
ok( $config->minlength == 200, 'minlength' );
ok( $config->limit == 0, 'limit' );
ok( $config->minintron == 10, 'minintron' );
ok( $config->outdir eq 'share/asn1val', 'outdir' );
ok( $config->discrep eq 'share/discrep.txt', 'discrep' );
ok( $config->verbosity == 1, 'verbosity' );

my @features = $config->feature;
ok( @features == 4, 'feature' );
