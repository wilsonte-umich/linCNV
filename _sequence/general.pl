use strict;
use warnings;

#----------------------------------------------------------
# sequence manipulations
#----------------------------------------------------------
# reverse complement reads in place (i.e. by reference)
sub rc { 
    my ($seq) = @_;
    $$seq = reverse $$seq;
    $$seq =~ tr/ATGC/TACG/;
}

# get rightmost mapped read position in reference genome from POS and CIGAR
sub getEnd { 
    my ($start, $cigar) = @_;
    $cigar =~ s/\d+S//g;
    my $end = $start - 1;
    while ($cigar =~ (m/(\d+)(\w)/g)) {
        $2 eq "I" or $end += $1;
    }
    return $end;
}

1;

