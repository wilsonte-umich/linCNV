use strict;
use warnings;

use vars qw(%options);

push @{$options{segment}}, ( 
    {           
        long        => 'transition-prob',
        short       => 's',
        type        => 'dbl',            
        required    => 0,
        default     => 1e-6,
        message     => "HMM transition probability",
    },
);

