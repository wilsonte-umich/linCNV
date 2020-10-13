use strict;
use warnings;

# prepare alignments in bam file for crosstab by parsing to BED

# constants
use constant {
    QNAME => 0, # sam fields
    FLAG => 1,
    RNAME => 2,
    POS => 3,
    MAPQ => 4,
    CIGAR => 5,
    RNEXT => 6,
    PNEXT => 7,
    TLEN => 8,
    SEQ => 9,
    QUAL => 10,
	TAGS => 11
};

# bam2bed
while(my $line = <STDIN>){
	my @f = split("\t", $line, 12);
	$f[RNAME] =~ m/_/ and next; # only use canonical chromosomes
	$f[TAGS] =~ m/CB:Z:(\S+)/;
    
    # need to re-run linCNV align after bug fix, should never get negative TLEN
    $f[TLEN] <= 0 and next;
    
    print join("\t", $f[RNAME], $f[POS]-1, $f[POS] + $f[TLEN] - 1, $1, 0, '.'), "\n";
}

