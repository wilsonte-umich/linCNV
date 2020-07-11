use strict;
use warnings;

use vars qw(%options);

push @{$options{bin}}, (
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
        long        => 'blacklist-file',
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
    {           
        long        => 'cell-ranger-dir',
        short       => 'R',
        type        => 'str',            
        required    => 1,
        default     => undef,
        message     => "directory with prior Cell Ranger DNA output",
    },  
    {           
        long        => 'min-mapq',
        short       => 'Q',
        type        => 'int',            
        required    => 0,
        default     => 5,
        message     => "reject read alignments with lower than this MAPQ",
    },
    {           
        long        => 'include-Y',
        short       => 'Y',
        type        => '',            
        required    => 0,
        default     => undef,
        message     => "include the Y chromosome in the output [chrY omitted]",
    },    
    {           
        long        => 'weight-per-cell',
        short       => 'w',
        type        => 'int',            
        required    => 0,
        default     => 10,
        message     => "average number of DNA fragments per bin per cell",
    },
);

1;

