#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use Bio::WGS2NCBI::Seq;

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

my $trunc = $seq->trunc(1,4);
ok( $trunc->seq eq 'AGC' );
ok( $trunc->revcom->seq eq 'GCT' );