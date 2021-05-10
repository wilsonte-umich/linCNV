use strict;
use warnings;

#----------------------------------------------------------------
# only pass reads from accepted cells (good or wavy)
# based on cell list established by filter_reads.pl
# caller or receiver must enforce any needed chromosome restrictions
#----------------------------------------------------------------

# load dependencies
use File::Basename;
my $scriptDir = dirname(__FILE__);
require "$scriptDir/../common/utilities.pl";

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
    #-------------
    TAGS => 11, # non-standard fields, also held in @aln
};

# load resource data
use vars qw(%cells @isAccepted);         
loadCells10X();
loadAcceptedCells();

while(my $line = <STDIN>){
    
    # pass header lines as is
    unless($line =~ m/^\@/){

        # only allow good and wavy cells
        my @aln = split("\t", $line, TAGS + 1);        
        $aln[TAGS] =~ m/CB:Z:(\S+)/ or next; # the 10X error-corrected, validated cell barcode
        ($cells{$1} and $isAccepted[$cells{$1}]) or next;
    }
    
    print $line;     
}

1;

