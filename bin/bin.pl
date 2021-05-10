use strict;
use warnings;

#----------------------------------------------------------------
# declare variable-length bins by considering reads in all cells
# each bin has a fixed read count when summed over all cells
#----------------------------------------------------------------

# initialize reporting
our $script = $ENV{FORCE_SCRIPT} || 'bin';
our $error  = "$script error";

# load dependencies
use File::Basename;
my $scriptDir = dirname(__FILE__);
require "$scriptDir/../_workflow/workflow.pl";
require "$scriptDir/../_numeric/utilities.pl";
require "$scriptDir/../common/utilities.pl";
resetCountFile();

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
#    #-------------
    TAGS => 11, # non-standard fields, also held in @aln
    CELL_ID => 0, # 0..nCells
    END_POS => 1,
};

# set flags to recover forward read of a quality proper pair,
# not secondary, supplemental or duplicate
my $minMapQ = $ENV{MIN_MAPQ} || 5;
my $f = 0x001 + # the read is paired in sequencing
        0x002;  # the read is mapped in a proper pair
my $F = 0x004 + # the query sequence itself is unmapped
        0x008 + # the mate is unmapped
        0x010 + # strand of the query (1 for reverse)
        0x100 + # the alignment is not primary
        0x200 + # the read fails platform/vendor quality checks
        0x400 + # the read is either a PCR or an optical duplicate
        0x800;  # supplementary alignment

# determine if input bam reheader action is needed to add "chr"
# script assumes that resource, exclusion and other files are "chr1" format
my $testChrom = qx/samtools view $ENV{BAM_FILE} | head -n 1 | cut -f 3/;
my $reheader = '';
if($testChrom !~ m/^chr/){
    my $bamHeader = qx/samtools view -H $ENV{BAM_FILE}/;
    $bamHeader =~ s/(\@SQ\tSN:)/$1chr/g; # OK, since will only use canonical chroms
    my $reheaderFile = "$ENV{DATA_PREFIX}.bin.reheader.sam";
    open my $headerH, ">", $reheaderFile or die "$error: could not open $reheaderFile: $!\n";
    print $headerH $bamHeader;
    close $headerH;
    $reheader = " | samtools reheader $reheaderFile - "; # will swap chrom names on the fly
}

# load resource data
use vars qw(%cells @allCells $nCells @acceptedCellIds @isAccepted $nAcceptedCells
            @chroms %chroms %chromNs $nChroms);
loadChromList();            
my $binWeightAllCells = loadCells10X();
$ENV{ACCEPTED_CELLS} and $binWeightAllCells = loadAcceptedCells();
my $ALL = @allCells; # index for holding the count of all cells (one higher than highest cell_id)

# process data by chromosome over multiple parallel threads
runChildThreads_array(\&bin_chromosome, @chroms);

# print summary information
printCount($nChroms, 'nChroms', 'chromosomes in output');
printCount($nCells,  'nCells',  'total cells as determined by CellRanger');
if($nAcceptedCells){
    printCount($nAcceptedCells, 'nAcceptedCells', 'cells marked as accepted and used during binning');
}

# child process to parse Cell Ranger BAM file on one chromosome
sub bin_chromosome {
    my ($childI) = @_;
    my $chrom = $chroms[$childI];

    # working variables
    my (@alns, @nextBinCount, @binSizes) = ();
    my ($wrkCount, $maxEnd,
        $prevBinEnd, $brkStart, $brkEnd,
        $nCountedPairs, $nBins) = (0) x 20;
    
    # output file handles
    my $paddedChildI = sprintf("%02d", $childI);
    my $binsFile = "$ENV{BIN_PREFIX}.counts.$paddedChildI.gz";
    open my $binFileH, "|-", "gzip -c > $binsFile" or die "$error: could not open $binsFile: $!\n";
    print $binFileH "#".join("\t", qw(chrom start end), @allCells), "\n"; # cell ids range from 0..nCells-1
    
    # open input bam stream
    my $viewChrom = $chrom;
    if($reheader){ $viewChrom =~ s/^chr// }
    my $bamStream = "samtools view -b -q $ENV{MIN_MAPQ} -f $f -F $F $ENV{BAM_FILE} $viewChrom ".
                    " $reheader ".
                    " | bedtools intersect -v -a stdin -b $ENV{EXCLUSIONS_FILE} ".
                    " | samtools view - ".
                    ""; #" | head -n 1000000";   
    open my $bamH, "-|", $bamStream or die "$error: could not open bam stream: $bamStream";
    while(my $line = <$bamH>){
    
        # finish filtering to countable reads
        my @aln = split("\t", $line, TAGS + 1);
        $aln[TLEN] > 0 or next; # why wouldn't it be at this point?
        
        # only allow true cells to continue
        $aln[TAGS] =~ m/CB:Z:(\S+)/ or next; # the 10X error-corrected, validated cell barcode
        defined $cells{$1} or next;  # true if Cell Ranger said it was a cell      
        $nAcceptedCells and ($isAccepted[$cells{$1}] or next); # ignore reads from cells not accepted by user
        $aln[CELL_ID] = $cells{$1};  # convert sequence barcode to CellRanger 0-based numeric cell_id
        $nCountedPairs++;
    
        # parse additional fragment information
        $aln[END_POS] = $aln[POS] + $aln[TLEN] - 1; # the end coordinate of the fragment
        $maxEnd >= $aln[END_POS] or $maxEnd = $aln[END_POS];
        $wrkCount++;              

        # as needed, handle variable-width/equal-weight bin processing
        if($wrkCount >= $binWeightAllCells){ # this fragment is guaranteed to create a break condition

            # this is the _first_ fragment that will create a break condition (index fragment)
            if(!$brkEnd){
                $brkStart = $aln[POS]; # break cannot occur before this fragment starts
                $brkEnd = $maxEnd;     # but could break as far as the end of any fragment until now        

            # this is the first fragment that CANNOT contribute to the previous bin
            } elsif($aln[POS] > $brkEnd){        

                # initialize counts
                my @binCount = @nextBinCount; # the count carried over from the last bin's break
                @nextBinCount = map { 0 } @allCells, $ALL; # the count to carry over to the next bin
                my ($brkCount, $binEnd, %brkSpans) = (0, 0);

                # examine every fragment that might contribute to the bin
                foreach my $x(@alns){
        
                    # these fragments are entirely within this bin, before the break
                    if($$x[END_POS] < $brkStart){
                        $binCount[$ALL]++; # just count them
                        $binCount[$$x[CELL_ID]]++;

                    # these fragments are likely to be split by the break
                    } else {
                        my $inc = 1 / $$x[TLEN];
                        my $CI = $$x[CELL_ID];
                        my ($minPos, $maxPos, $brkWeight);

                        # some fragments must have been encountered BEFORE the index fragment
                        if($$x[POS] < $brkStart){
                            my $bin = ($brkStart - $$x[POS]) * $inc;
                            $binCount[$ALL] += $bin;
                            $binCount[$CI]  += $bin; # these bases always in this bin
                            $brkCount += 1 - $bin;
                            ($minPos, $maxPos, $brkWeight) = ($brkStart, $$x[END_POS], 1 - $bin); # these bases in break region
  
                        # some fragments must have been encountered AFTER the index fragment
                        } elsif($$x[END_POS] > $brkEnd){
                            my $next = ($$x[END_POS] - $brkEnd) * $inc;
                            $nextBinCount[$ALL] += $next; # these bases always in next bin
                            $nextBinCount[$CI]  += $next;
                            $brkCount += 1 - $next;
                            ($minPos, $maxPos, $brkWeight) = ($$x[POS], $brkEnd, 1 - $next); # these bases in break region

                        # some fragments may reside entirely in the break region
                        } else {
                            $brkCount++;
                            ($minPos, $maxPos, $brkWeight) = ($$x[POS], $$x[END_POS], 1);
                        }
                        push @{$brkSpans{$CI}}, [$minPos, $maxPos, $brkWeight];
                    }
                }
        
                # interpolate the break position based on actual and needed break region counts
                # this provides an estimate of the true break position, but it's good enough
                # the error results in a slight variance in total bin counts
                # but each sample count in the final bin is accurate
                $binEnd = roundCount($brkStart + ($brkEnd - $brkStart) * ($binWeightAllCells - $binCount[$ALL]) / $brkCount, 1);
                $binEnd - $prevBinEnd < 0 and die "$error: parsing resulted in negative bin size\n";

                # add the portions of fragments within the break region to this bin or next bin
                foreach my $CI(keys %brkSpans){
                    foreach my $brkSpan(@{$brkSpans{$CI}}){
                        my ($minPos, $maxPos, $brkWeight) = @$brkSpan;
                        if ($maxPos <= $binEnd){
                            $binCount[$ALL] += $brkWeight;
                            $binCount[$CI]  += $brkWeight;
                        } elsif($minPos > $binEnd){
                            $nextBinCount[$ALL] += $brkWeight;
                            $nextBinCount[$CI]  += $brkWeight;
                        } else {
                            my $bin = $brkWeight * ($binEnd - $minPos) / ($maxPos - $minPos);
                            $binCount[$ALL] += $bin;
                            $binCount[$CI]  += $bin;
                            $nextBinCount[$ALL] += 1 - $bin;
                            $nextBinCount[$CI]  += 1 - $bin;
                        }
                    }
                }
        
                # commit the newly finished bin
                print $binFileH join("\t", $chrom, $prevBinEnd, $binEnd);
                foreach my $CI(@allCells){
                    my $count = roundCount2($binCount[$CI] || 0);
                    $nAcceptedCells and !$isAccepted[$CI] and $count = "NA"; # override unaccepted cells to NA, but do print them
                    print $binFileH "\t".$count;
                }
                print $binFileH "\n";
                push @binSizes, $binEnd - $prevBinEnd;
                $prevBinEnd = $binEnd;
        
                # handle rare regions with extremely high read density
                # "bad" genome regions that have not been excluded by user
                # not common if using the ENCODE genome bad region files
                while($nextBinCount[$ALL] >= $binWeightAllCells){
                    my $badStart = min($binEnd + 1, $brkEnd);
                    my $badFrac = $binWeightAllCells / $nextBinCount[$ALL];
                    $binEnd = roundCount($badStart + ($brkEnd - $badStart) * $badFrac, 1);
                    $binEnd - $prevBinEnd < 0 and die "$error: parsing resulted in negative bin size while fixing a bad region\n";
                    print $binFileH join("\t", $chrom, $prevBinEnd, $binEnd);
                    foreach my $CI($ALL, @allCells){
                        my $count = ($nextBinCount[$CI] || 0) * $badFrac;
                        $nextBinCount[$CI] -= $count;                        
                        unless($CI == $ALL){
                            $count = roundCount2($count);
                            $nAcceptedCells and !$isAccepted[$CI] and $count = "NA"; 
                            print $binFileH "\t".$count;
                        }
                    }
                    print $binFileH "\n";
                    $prevBinEnd = $binEnd;
                }
        
                # reset for the next bin
                ($wrkCount, $maxEnd, $brkEnd) = ($nextBinCount[$ALL], $aln[END_POS], 0);
                @alns = ();
                $nBins++;
            }
        }
        
        # collect the set of working fragments
        push @alns, \@aln;        
    }
    
    # clean up
    close $bamH;        
    close $binFileH;
    
    # report chromosome-specific counts
    my $medianBinSize = roundCount(median(@binSizes),  1);
    printCount($nCountedPairs, 'nCountedPairs', "$chrom: proper read pairs accepted for counting");
    printCount($nBins,         'nBins',         "$chrom: output bins");
    printCount($medianBinSize, 'medianBinSize', "$chrom: median bin size");
}

1;

# both reads of a pair are marked as duplicates
#A00228:250:H3TVLDSXX:1:1110:25138:24925 1024,1024
#A00228:250:H3TVLDSXX:1:1139:28248:11710 1024,1024
#A00228:250:H3TVLDSXX:1:1146:8477:8469   1024,1024
#A00228:250:H3TVLDSXX:1:1151:9643:32659  1024,1024

# the two reads for one of those pairs
#A00228:250:H3TVLDSXX:1:1139:28248:11710 1123    22      16057320        27      84M     =       16057594        374    GTAAAATGTAGATTTTGGCACACGGGGACCTAAGACAATCCTGCCTGAAACAGGGCCTGGCATCCCATCAGTGCTCAATAAACA     FFFFF:F:FFFFFFFFF:FFFF,,FF:FFFFFFFFFF:FF:FFF,FFFFF:FFFFF,FFF,FFF,FF,,FF,FFFFFFFFFFFF    NM:i:0  MD:Z:84 AS:i:84 XS:i:84 XA:Z:14,-19785603,84M,0;       CR:Z:AAGCCGCGTTGACGGA    CY:Z:FFFFFFFFFFFFFFFF   CB:Z:AAGCCGCGTTGACGGA-1 BC:Z:ACCCTCCT   QT:Z:,F::FFF:   GP:i:698302392 MP:i:698302766   RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1
#A00228:250:H3TVLDSXX:1:1139:28248:11710 1171    22      16057594        27      100M    =       16057320        -374   TGAGAAGCTGTCCTGGTCTAAGTATCTTTCTCCCATTTTACATAAAGGAATACACAGTGTCAGAAGGAGGACCTGTGTCCAGCCCCTGTGTTCCCCACTT     FFFFFFFF,FF:FFFF:FFFFFF:FFFFF:FFFFFF:,,,FF,::FF:FF,FFFFFF:FF:,FF:FFFFFFF:FFF,:FFFFF:FF,FFF:FFFFFFF:F    NM:i:1  MD:Z:99C0       AS:i:99XS:i:94  XA:Z:14,+19785313,100M,2;       CR:Z:AAGCCGCGTTGACGGA   CY:Z:FFFFFFFFFFFFFFFF   CB:Z:AAGCCGCGTTGACGGA-1 BC:Z:ACCCTCCT   QT:Z:,F::FFF:   GP:i:698302766  MP:i:698302392  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1

# three different reads with the same barcode and start position
# 1123 = 99 | 1024
# thus ONE of the duplicate is NOT marked as duplicate (as should be, but documentation was badly worded)
#A00228:250:H3TVLDSXX:1:1139:27434:9925  99      22      16057320        27      84M     =       16057594        374    GTAAAATGTAGATTTTGGCACACGGGGACCTAAGACAATCCTGCCTGAAACAGGGCCTGGCATCCCATCAGTGCTCAATAAACA     FFF:FFFFF:FFF:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF,FFFF:F,,F,:FF,
#NM:i:0  MD:Z:84 AS:i:84 XS:i:84 XA:Z:14,-19785603,84M,0;       CR:Z:AAGCCGCGTTGACGGA    CY:Z:FFFFFFFFFFFFFFFF
#CB:Z:AAGCCGCGTTGACGGA-1 BC:Z:ACCCTCCT   QT:Z:FFFF,F,,   GP:i:698302392 MP:i:698302766   RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1

#A00228:250:H3TVLDSXX:1:1139:28248:11710 1123    22      16057320        27      84M     =       16057594        374    GTAAAATGTAGATTTTGGCACACGGGGACCTAAGACAATCCTGCCTGAAACAGGGCCTGGCATCCCATCAGTGCTCAATAAACA     FFFFF:F:FFFFFFFFF:FFFF,,FF:FFFFFFFFFF:FF:FFF,FFFFF:FFFFF,FFF,FFF,FF,,FF,FFFFFFFFFFFF
#NM:i:0  MD:Z:84 AS:i:84 XS:i:84 XA:Z:14,-19785603,84M,0;       CR:Z:AAGCCGCGTTGACGGA    CY:Z:FFFFFFFFFFFFFFFF
#CB:Z:AAGCCGCGTTGACGGA-1 BC:Z:ACCCTCCT   QT:Z:,F::FFF:   GP:i:698302392 MP:i:698302766   RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1

#A00228:250:H3TVLDSXX:4:1315:24505:4507  1123    22      16057320        27      84M     =       16057594        374    GTAAAATGTAGATTTTGGCACACGGGGACCTAAGACAATCCTGCCTGAAACAGGGCCTGGCATCCCATCAGTGCTCAATAAACA     FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
#NM:i:0  MD:Z:84 AS:i:84 XS:i:84 XA:Z:14,-19785603,84M,0;       CR:Z:AAGCCGCGTTGACGGA    CY:Z:FFFFFFFFFFFFFFFF
#CB:Z:AAGCCGCGTTGACGGA-1 BC:Z:GTTGCAGC   QT:Z:FFFFFFFF   GP:i:698302392 MP:i:698302766   RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:4


# a different set: here, a position duplicate from another CELL is NOT marked duplicate
#A00228:250:H3TVLDSXX:1:1146:8477:8469   1171    22      16053829        27      100M    =       16053921        -8     TACCAGAGGCATGGGGTGAGGCATGTTGCAGGTCAAGGACCAGGGCCATCTCACTGCCTGAGCCCATGGACTGGCTCAGGGGTCTGTCAGATGATTCTAG     FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF    NM:i:1  MD:Z:34G65      AS:i:95XS:i:90  XA:Z:14,+19789081,100M,2;       CR:Z:ACAGCTAAGGCATTGG   CY:Z:FFFFFFFFFFFFFFFF
#CB:Z:ACAGCTAAGGCATTGG-1 BC:Z:ACCCTCCT   QT:Z:FFFFFFFF   GP:i:698299001  MP:i:698298993  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1
#
#A00228:250:H3TVLDSXX:1:1146:8160:8484   147     22      16053829        27      100M    =       16053921        -8     TACCAGAGGCATGGGGTGAGGCATGTTGCAGGTCAAGGACCAGGGCCATCTCACTGCCTGAGCCCATGGACTGGCTCAGGGGTCTGTCAGATGATTCTAG     FFFFFFF,FFFFFFFF:FFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF    NM:i:1  MD:Z:34G65      AS:i:95XS:i:90  XA:Z:14,+19789081,100M,2;       CR:Z:ACAGCTAAGGCATTGG   CY:Z:FFFFFFFFFFFFFFFF
#CB:Z:ACAGCTAAGGCATTGG-1 BC:Z:ACCCTCCT   QT:Z:FFFFFFFF   GP:i:698299001  MP:i:698298993  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1
#
#A00228:250:H3TVLDSXX:1:2628:1479:16297  1171    22      16053829        27      100M    =       16053921        -8     TACCAGAGGCATGGGGTGAGGCATGTTGCAGGTCAAGGACCATGGCCATCTCACTGCCTGAGCCCATGGACTGGCTCAGGGGTCTGTCAGATGATTCTAG     FF,F:FFFF:FFFF:FF,,FF,FFFFFFFF,FFFFFFFFFFF,:::FFFFFFFF,,FF::FFF:FFFFFFFFFFFFFFFFFFFFFFFFF,,FFFF:FFFF    NM:i:2  MD:Z:34G7G57    AS:i:90XS:i:85  XA:Z:14,+19789081,100M,3;       CR:Z:ACAGCTAAGGCATTGG   CY:Z:FFFF,FFFFF,FF:FF
#CB:Z:ACAGCTAAGGCATTGG-1 BC:Z:ACCCTCCT   QT:Z:FFFF,F,F   GP:i:698299001  MP:i:698298993  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1
#
#A00228:250:H3TVLDSXX:2:1373:28854:34428 147     22      16053829        50      100M    =       16053669        -260   TACCAGAGGCATGGGGTGAGGCATGTTGCAGGTCGAGGACCAGGGCCATCTCACTGCCTGAGCCCATGGACTGGCTCAGGGGTCTGTCAGATGATTCTAG     :FFF:FFFFFF::FFF:FFF,FFF:FFFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF    NM:i:0  MD:Z:100        AS:i:10XS:i:95  XA:Z:14,+19789081,100M,1;       CR:Z:CGGAGCTAGATCCGAG   CY:Z:FFF:FFFFFFFFFFFF
#CB:Z:CGGAGCTAGATCCGAG-1 BC:Z:CAATGGAG   QT:Z::::FF,,,   GP:i:698299001  MP:i:698298741  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:2

