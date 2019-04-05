#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my $Verbosity = 1;

sub LOG ($$) {
    my ($msg,$method) = @_;
    my ( $package, $file1up, $line1up, $sub ) = caller( 2 );
    my ( $pack0up, $file, $line, $sub0up )    = caller( 1 );
    my $log = sprintf( "%s %s [%s %s] - %s\n", uc $method, $sub || '', $0, $line, $msg );
    print STDERR $log;
}
sub DEBUG ($) { LOG shift, 'DEBUG' if $Verbosity >= 3 }
sub INFO ($)  { LOG shift, 'INFO'  if $Verbosity >= 2 }
sub WARN ($)  { LOG shift, 'WARN'  if $Verbosity >= 1 }

sub check_args {
	my $infile;
	my $evidence = 'paired-ends';
	my $gap_type = 'scaffold';
	my $linkage  = 'yes';
	my $organism = 'Ophiophagus hanna';
	my $taxid    = 8665;
	my $name     = 'PRJNA201683';
	my $center   = 'NaturalisBC';
	GetOptions(
		'verbose+'   => \$Verbosity,
		'infile=s'   => \$infile,
		'evidence=s' => \$evidence,
		'gap_type=s' => \$gap_type,
		'organism=s' => \$organism,
		'taxid=i'    => \$taxid,
		'name=s'     => \$name,
		'center=s'   => \$center,
	);
	return 
		'infile'   => $infile, 
		'evidence' => $evidence, 
		'gap_type' => $gap_type,
		'linkage'  => $linkage,
		'organism' => $organism,
		'taxid'    => $taxid,
		'name'     => $name,
		'center'   => $center;
}

sub read_fasta {
    my ( $fh, $pos ) = @_;
    
    # advance to the current position
    seek $fh, $pos, 0;
    DEBUG "reading starts at $pos";
    
    # these will store what's on the def line and the subsequent seq data
    my ( $chr, $seq );
    LINE: while(<$fh>) {
        chomp;
        
        # we have a definition line, break if we've already processed a defline
        if ( />(.+)$/ ) {
            last LINE if $chr;
            $chr = $1;
            DEBUG "going to concatenate sequence data for $chr";
        }
        else {
            $seq .= $_;
            
            # we store the current position here so that it is set to 
            # just before the defline when we encounter it
            $pos = tell $fh;
        }       
    }   
    DEBUG "reading ended at $pos";
    return $pos, $chr, \$seq;
}

sub write_agp_header {
	my ( $organism, $taxid, $name, $center ) = @_;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my @month = qw(January February March April May June July August September October November December);
	my $date = sprintf("%02d",$mday) . '-' . $month[$mon] . '-' . ( 1900 + $year );
	print <<"HEADER";
##agp-version	2.0
# ORGANISM: $organism
# TAX_ID: $taxid
# ASSEMBLY NAME: $name
# ASSEMBLY DATE: $date
# GENOME CENTER: $center
HEADER
}

sub write_gaps {
	my ( $chr, $seq, %args ) = @_;
	my $pos = 0;
	my $max = length($$seq) - 1;
	my $counter;
	my ( $type, $linkage, $evidence ) = @args{qw(gap_type linkage evidence)};
	my $template = "$chr\t%d\t%d\t%d\tN\t%d\t$type\t$linkage\t$evidence\n";
	while( ( my $i = index( $$seq, 'N', $pos ) ) != -1 ) {
		my $j = $i;
		STRETCH: while ( $j <= $max ) {
			last STRETCH if 'N' ne substr $$seq, $j, 1;
			$j++;
		}
		printf $template, $i + 1, $j, ++$counter, $j - $i;
		$pos = $j;
	}
}

sub splice_abutting_gaps {
	my $seqref = shift;
	my $i = 0;
	LEADING: while( $i < length($$seqref) ) {
		last LEADING if substr( $$seqref, $i, 1 ) ne 'N';
		$i++;
	}
	my $j = length($$seqref) - 1;
	TRAILING: while( $j >= 0 ) {
		last TRAILING if substr( $$seqref, $j, 1 ) ne 'N';
		$j--;
	}
	return \( substr( $$seqref, $i, ( $j - $i ) + 1 ) );
}

sub main {
	my %args = check_args();
    
    # open the fasta file handle
    open my $fastaFH, '<', $args{'infile'} or die $!;
    my $fastaPos = 0;
    
    write_agp_header(@args{qw(organism taxid name center)});
    
    # iterate over the genome
    while( not eof $fastaFH ) {
    
        # advance to the next scaffold/chromosome 
        my ( $chr, $seq );
        ( $fastaPos, $chr, $seq ) = read_fasta( $fastaFH, $fastaPos );
        $chr =~ s/\|/-/;

		# write the AGP rows
        write_gaps( $chr, splice_abutting_gaps($seq), %args );
    }   
}
main();
