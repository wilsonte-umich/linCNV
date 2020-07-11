use strict;
use warnings;

use vars qw(%options);

push @{$options{normalize}}, ( 
    {           
        long        => 'min-modal-cn',
        short       => 'c',
        type        => 'dbl',            
        required    => 0,
        default     => 0.25,
        message     => "only use bins with at least this fractional modal copy number",
    },
    {           
        long        => 'min-mappability',
        short       => 'b',
        type        => 'dbl',            
        required    => 0,
        default     => 0.25,
        message     => "only use bins with at least this fractional mappability",
    },
    {           
        long        => 'max-excluded-bases',
        short       => 'x',
        type        => 'int',            
        required    => 0,
        default     => 1000,
        message     => "only use bins containing fewer than this many excluded bases",
    },
    {           
        long        => 'min-allele-depth',
        short       => 'a',
        type        => 'dbl',            
        required    => 0,
        default     => 2,
        message     => "only use cells with at least this many average reads per allele",
    },  
);

1;

