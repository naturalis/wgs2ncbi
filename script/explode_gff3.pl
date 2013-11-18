#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my ( $source, $gff3, @features, $dir );
GetOptions(
	'source=s'  => \$source,
	'gff3=s'    => \$gff3,
	'feature=s' => \@features,
	'dir=s'     => \$dir,
);

warn "source: $source\n";
warn "gff3: $gff3\n";
warn "feat: @features\n";
warn "dir: $dir\n";

my %types = map { $_ => 1 } @features;
my $chr_idx  = 0;
my $src_idx  = 1;
my $type_idx = 2;
my @queue;
my %handles;

open my $fh, '<', $gff3 or die;
LINE: while(<$fh>) {
	next LINE if /^#/;
	my @line = split /\t/, $_;
	next LINE if not $line[$src_idx] or $line[$src_idx] ne $source;
	next LINE if not $line[$type_idx] or not $types{$line[$type_idx]};
	my $chr = $line[$chr_idx];
	
	# we keep up to 100 file handles (which is max on POSIX). once we
	# reach that number we close existing ones on a first in, first out basis
	if ( not $handles{$chr} ) {
		if ( scalar @queue == 100 ) {
			my $fifo = shift @queue;
			close $handles{$fifo};
			delete $handles{$fifo};
		}
		open $handles{$chr}, '>', "${dir}/${chr}.gff3" or die $!;
		push @queue, $chr;
	}
	my $fh = $handles{$chr};
	print $fh $_;
}