use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'          => 'Bio::WGS2NCBI',
    'VERSION_FROM'  => 'lib/Bio/WGS2NCBI.pm',
    'ABSTRACT_FROM' => 'lib/Bio/WGS2NCBI.pm',
    'AUTHOR'        => 'Rutger Vos <Rutger.Vos@naturalis.nl>',
    'LICENSE'       => 'perl',
    'EXE_FILES'     => [ 'script/wgs2ncbi' ],
);
