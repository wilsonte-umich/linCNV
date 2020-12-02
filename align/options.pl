use strict;
use warnings;

use vars qw(%options);

push @{$options{align}}, (
    {           
        long        => 'input-dir',
        short       => 'i',
        type        => 'str',            
        required    => 1,
        default     => undef,
        message     => "every sub-directory of input-dir will be searched for cell FASTQ files",
    },
    {           
        long        => 'genome',
        short       => 'g',
        type        => 'str',            
        required    => 1,
        default     => undef,
        message     => "the reference genome to which reads were aligned",
    },
    {           
        long        => 'genome-dir',
        short       => 'G',
        type        => 'str',            
        required    => 1,
        default     => undef,
        message     => "the directory in which reference genome files can be found",
    },    
    {           
        long        => 'min-mapq',
        short       => 'Q',
        type        => 'int',            
        required    => 0,
        default     => 5,
        message     => "reject read alignments with lower than this MAPQ",
    },
);

