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
);

