use strict;
use warnings;

use vars qw(%options);

push @{$options{shared}}, ( 
    {           
        long        => 'genome',
        short       => 'g',
        type        => 'str',            
        required    => 1,
        default     => undef,
        message     => "the reference genome to which reads were aligned",
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

