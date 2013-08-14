#!/usr/bin/perl
use strict;
use warnings;
use Bio::WGS2NCBI::Util::Logger;
use Bio::WGS2NCBI::Util::Config;
use Bio::WGS2NCBI::Util::Factory;

my $conf = 'Bio::WGS2NCBI::Util::Config';
my $fac  = Bio::WGS2NCBI::Util::Factory->new;

sub get_non_missing_index {
	my ( $seq, $reverse, $missing ) = @_;
	$missing = 'N' if not defined $missing;
	my $seq = $self->{'seq'};
	if ( $reverse ) {
		for ( my $i = $seq->length; $i >= 1; $i-- ) {
			return $i - 1 if $seq->subseq($i,$i) ne $missing;
		}
	}
	else {
		for my $i ( 1 .. $seq->length ) {
			return $i - 1 if $seq->subseq($i,$i) ne $missing;
		}	
	}
}

sub make_defline {
	my $seq = shift;
	my $line = $seq->display_id;
	for my $key ( $conf->defline->keys ) {
		$line .= ' [' . $key . '=' . $conf->defline->$key . ']';
	}
	return $line;
}

sub annotations {

}

sub validate {

}

sub data {
    
    # open the genome file handles
	HANDLE: for my $seq_handle ( $conf->open_handles( 'project' => 'genome_files' ) ) {
		my $seq_reader = $fac->create_seqio( 
			'-format' => $conf->project->genome_format, 
			'-fh'     => $seq_handle,
		);
	
		# counts processed sequences
		my $i;
	
		# we will re-use these handles
		my ( $tblFH, $seq_writer );
	
		# iterate over the genome
		SEQ: while( my $seq = $seq_reader->next_seq ) {
			$i++;
			last HANDLE if $conf->internal->limit and $i == $conf->internal->limit;
	
			# get length and id
			my $length = $seq->length;
			my $id     = $seq->display_id;
			$id =~ s/\|/-/; # XXX deleteme
					
			# compute starting offset, if any
			my $offset;
			if ( $offset = get_non_missing_index($seq) ) {
				INFO "leading ${offset}bp gap in $id, will strip this and apply offset";
			}
		
			# check what we have left
			my $last_non_missing_index = get_non_missing_index($seq,'reverse');
			$length -= ( $length - 1 - $last_non_missing_index ) + $offset;
			if ( $length < $conf->filter->minlength ) {
				WARN "remaining seq $id is too short ($length bp), skipping";
				$i--;
				next SEQ;
			}
		
			# open output handles
			if ( ( $i % $conf->internal->chunksize ) == 1 ) {
				INFO "starting next chunk of sequences and feature tables";
								
				# generate a new name indicating the range
				my $upper = $i + $conf->internal->chunksize - 1;
				my $stem = $conf->output->tblfasta_dir . "/combined_${i}-${upper}";
				
				# open a new feature table handle
				open $tblFH, '>', "${stem}.tbl" or FATAL "Can't write ${stem}.tbl: $!"; 
				
				# open a new FASTA writer  
				open my $fh, '>', "${stem}.fsa" or FATAL "Can't write ${stem}.fsa: $!";
				$seq_writer = $fac->create_seqio( 
					'-format' => 'fasta', 
					'-fh'     => $fh,
				);
			}
		
			# get the features for that scaffold/chromosome, if we have them
			my $gff3 = $conf->output->gff3_dir . "/${id}.gff3";
			if ( -e $gff3 ) {
				open my $gff3FH, '<', $gff3 or FATAL $!;  
				my $features = read_features( 
					$gff3FH,  # file handle of focal file
					$id,      # scaffold/chromosome ID
					$counter, # counter for generating locus tags
					$config,  # config object
					$seq,     # string reference of raw sequence
					$offset,  # start offset
					$last_non_missing_index, # last true seq character
				);  
				$features->offset($offset);     
				write_features( $features, $tblFH );        
			}
			else {
				print $tblFH '>Features ', $id, "\n";
			}
		
			# write the truncated sequence with expanded definition line
			my $truncated = $seq->trunc($offset,$last_non_missing_index);
			$truncated->display_id(make_defline($truncated));
			$seq_writer->write_seq($truncated);      
		} 
    }  
}
