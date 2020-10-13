use strict;
use warnings;

use vars qw(%options);


push @{$options{genome}}, ( 
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
        long        => 'gap-file',
        short       => 'X',
        type        => 'str',            
        required    => 0,
        default     => undef,
        message     => "BED file with genomic gap regions (i.e. N base spans)",
    },
    {           
        long        => 'bad-regions-file',
        short       => 'B',
        type        => 'str',            
        required    => 1,
        default     => undef,
        message     => "BED file with genomic regions to exclude from analysis",
    },
    {           
        long        => 'gc-file',
        short       => 'm',
        type        => 'str',            
        required    => 1,
        default     => undef,
        message     => "BED file with fraction GC bases per fixed bin (e.g. 1kb)",
    },    
    {           
        long        => 'mappability-file',
        short       => 'M',
        type        => 'str',            
        required    => 1,
        default     => undef,
        message     => "BED file with fraction mappability score per fixed bin (e.g. 1kb)",
    },
);

push @{$options{shared}}, ( 
    {           
        long        => 'min-mapq',
        short       => 'Q',
        type        => 'int',            
        required    => 0,
        default     => 5,
        message     => "reject read alignments with lower than this MAPQ",
    },
    {           
        long        => 'ploidy',
        short       => 'P',
        type        => 'int',            
        required    => 0,
        default     => '2',
        message     => "expected autosome copy number per cell (e.g. 2=diploid)",
    },
);

1;

