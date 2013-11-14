package Bio::WGS2NCBI::Config;

sub new {
    my $class = shift;
    my %args = @_;
    my $self = \%args;
    return bless $self, $class;
}

sub prefix { shift->{'prefix'} || 'OPHA_' }

sub source { shift->{'source'} || 'maker' }

sub dir { shift->{'dir'} || 'submission' }

sub fasta { shift->{'fasta'} or die "Need FASTA file to operate on!" }

sub gff3 { shift->{'gff3'} or die "Need GFF3 file to operate on!" }

sub authority { shift->{'authority'} || 'gnl|NaturalisBC|' }

sub chunksize { shift->{'chunksize'} || 1 }

sub info { shift->{'info'} }

sub limit { shift->{'limit'} }

sub minlength { shift->{'minlength'} || 200 }

1;