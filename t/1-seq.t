#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use Bio::WGS2NCBI;
use Bio::WGS2NCBI::Seq;
use Bio::WGS2NCBI::Logger;

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

WARN "*** warning is expected here ***";
Bio::WGS2NCBI::mask_seq( $orig, '2..3' );
ok( $orig->seq eq 'ANNTCGATCG' );