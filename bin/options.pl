use strict;
use warnings;

use vars qw(%options);

push @{$options{bin}}, (
    {           
        long        => 'cell-ranger-dir',
        short       => 'R',
        type        => 'str',            
        required    => 1,
        default     => undef,
        message     => "directory with prior Cell Ranger DNA output",
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

