use strict;
use warnings;

use vars qw(%options);

push @{$options{crosstab}}, (
    {           
        long        => 'cell-ranger-dir',
        short       => 'R',
        type        => 'str',            
        required    => 1,
        default     => undef,
        message     => "directory with prior Cell Ranger DNA output",
    },                 
    {           
        long        => 'bin-size',
        short       => 'z',
        type        => 'int',            
        required    => 0,
        default     => 1e6,
        message     => "the size of the genome bins",
    },
    {           
        long        => 'manifest-file',
        short       => 'F',
        type        => 'str',            
        required    => 1,
        default     => undef,
        message     => "Michigan Advanced Genomics Core manifest file path",
    },
);

1;

#export GENOME_PREFIX="/home/wilsonte_lab/clubhouse/genomes/$GENOME/$GENOME"
#export BAD_REGIONS_DIR="/treehouse/wilsonte_lab/ssd/genomes/Blacklist/lists"
#export GAP_FILE="$GENOME_PREFIX.gap.bed"
#export MAPPABILITY_FILE="$GENOME_PREFIX.kmer_50.bin_1000.bed.gz"
#export GC_FILE="$GENOME_PREFIX.gc5Base.bin_1000.bed.gz"
#export BAD_REGIONS_FILE="$BAD_REGIONS_DIR/$GENOME-blacklist.v2.bed"