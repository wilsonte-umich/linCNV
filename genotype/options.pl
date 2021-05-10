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
);

push @{$options{genotype}}, ( 
    {           
        long        => 'will-impute',
        short       => 'w',
        type        => 'str',            
        required    => 0,
        default     => "TRUE",
        message     => "set to FALSE if genotypes will _not_ be manually subjected to imputation phasing",
    },
);

push @{$options{shared}}, (
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
        default     => 10,
        message     => "reject read alignments with lower than this MAPQ",
    }
);

1;

