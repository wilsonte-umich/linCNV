use strict;
use warnings;

use vars qw(%options);

push @{$options{scan}}, ( 
    {           
        long        => 'n-scan-bins',
        short       => 'S',
        type        => 'int',            
        required    => 0,
        default     => 100,
        message     => "sc0an the genome in windows comprising this many bins",
    },
);

1;

