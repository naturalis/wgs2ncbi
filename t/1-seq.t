#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use Bio::WGS2NCBI;
use Bio::WGS2NCBI::Seq;
use Bio::WGS2NCBI::Logger;

# this is just so that intentional warnings 
# are hidden from users running the unit tests
$Bio::WGS2NCBI::Logger::Verbosity = 0;

my $raw = 'AGCTCGATCG';
my $seq = Bio::WGS2NCBI::Seq->new(
	'-id'   => 'foo',
	'-seq'  => $raw,
	'-desc' => 'bar',
);

ok( $seq->seq eq $raw );
ok( $seq->length == 10 );
ok( $seq->id eq 'foo' );
ok( $seq->desc eq 'bar' );

my $trunc = $seq->trunc(1,3);
ok( $trunc->seq eq 'AGC' );
ok( $trunc->revcom->seq eq 'GCT' );

# test masking
my $orig = Bio::WGS2NCBI::Seq->new(
	'-id'   => 'foo',
	'-seq'  => 'AGCTCGATCG',
	'-desc' => 'bar',
);

$orig->mask( '2..3' );
ok( $orig->seq eq 'ANNTCGATCG' );

# test writing, when executed outside 'make test' this will show in the console
ok( $orig->write_fasta(\*STDOUT) );

# test reading
my $fasta = <<HERE;
>foo
AGCATGACATAGCGA
HERE
open my $fh, '<', \$fasta;
my $dataseq;
my $pos;
ok( ( $pos, $dataseq ) = Bio::WGS2NCBI::Seq->read_fasta( $fh ) );
ok( $dataseq->seq eq 'AGCATGACATAGCGA' );
ok( $pos == length($fasta) );
ok( $dataseq->get_non_missing_index == 0 );
ok( $dataseq->get_non_missing_index(1) == ( length($dataseq->seq) - 1 ) );

