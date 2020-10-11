use strict;
use warnings;

#----------------------------------------------------------------
# align a series of individual cells to genome and write output to a single stream
#----------------------------------------------------------------

# initialize reporting
our $script = 'align';
our $error  = "$script error";
my ($nCells) = (0) x 10;

# load dependencies
use File::Basename;
my $scriptDir = dirname(__FILE__);
require "$scriptDir/../_workflow/workflow.pl";
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
    #-------------
    _RNAME => 0, # grouping stream columns
    _GROUP_POS1 => 1,
    _GROUP_POS2 => 2,
    _POS => 3,
    _MAPQ => 4,
    _TLEN => 5,
    _DUP_RATE => 6,
    #-------------
    _PAIRED => 1,
    _PROPER_PAIR => 2,
    _PROPER_PAIR_BITS => 1 + 2,
};

# initialize outputs
my $bamListFile = "$ENV{LOGS_PREFIX}.$script.bam.list";
open my $bamListH, ">", $bamListFile or die "could not open $bamListFile for writing\n";
open my $cellH, ">", $ENV{CELL_LIST_FILE} or die "could not open $ENV{CELL_LIST_FILE} for writing\n";
print $cellH join(",", qw(barcode cell_id)), "\n";
my $logFile = "$ENV{LOGS_PREFIX}.$script.bwa.log";
unlink $logFile;

# get inputs
my @cellDirs = split(/\s+/, $ENV{CELL_DIRS});

# collect a common bam header for all cell output files
my $cellName = $cellDirs[0];
my $fastQs = join(" ", glob("$cellName/*.fastq.gz"));
my $bwa = "bwa mem -Y -t $ENV{N_CPU} $ENV{GENOME_DIR}/$ENV{GENOME}.fa $fastQs 2>/dev/null";
my $samtools = "samtools view -H";
my $bamHeader = qx/$bwa | $samtools/;

# align one cell at a time with parallelization
foreach my $cellName(@cellDirs){
    
    # add cell to cell list in abbreviated CellRanger format    
    $nCells++;
    my $cellN = $nCells - 1;
    print STDERR "  $cellName (cell #$cellN)\n";
    print $cellH join(",", $cellName, $cellN), "\n";
    
    # set files
    my @fastQs = glob("$cellName/*.fastq.gz");
    my $tmpBam = "$ENV{TMP_DIR_WRK}/cell_$cellN.sorted.bam";
    print $bamListH "$tmpBam\n";
    -e $tmpBam and next; # here to help with job recovery without remapping prior successes
    
    # align reads and remove duplicates
    # output does not retain alignments, just the needed information about unique source molecules
    my $logPrefix = "$ENV{LOGS_PREFIX}.$cellName.fastp";
    my $fastp = "fastp --thread $ENV{N_CPU} --in1 $fastQs[0] --in2 $fastQs[1] --stdout ".
                "--merge --include_unmerged --disable_quality_filtering ".
                "--html $logPrefix.html --json $logPrefix.json --report_title \"$cellName\" 2>/dev/null";
    #my $head = "head -n 100000";
    my $bwa = "bwa mem -p -Y -t $ENV{N_CPU} $ENV{GENOME_DIR}/$ENV{GENOME}.fa - 2>>$logFile";    
    my $parseBWA = "perl $ENV{PIPELINE_DIR}/align/parse_bwa.pl";
    my $sort = "sort --parallel=$ENV{N_CPU} -T $ENV{TMP_DIR_WRK} -S $ENV{MAX_SORT_RAM} --compress-program=pigz -k1,1 -k2,2n -k3,3n";
    my $group = "groupBy -g 1,2,3 -c 4,5,6,6 -o first,max,first,count";
    my $inStream = "$fastp | $bwa | $parseBWA | $sort | $group";
    open my $inH, "-|", $inStream or die "could not open bwa stream: $!\n";

    # process to bam for subsequent merging
    my $toBam = "samtools view --threads $ENV{N_CPU} -b -";
    my $bamSort = "samtools sort --threads $ENV{N_CPU} -m $ENV{RAM_PER_CPU} ".
                  "-T $ENV{TMP_DIR_WRK}/align_sort - 2>/dev/null"; # yes, we DO need to resort by POS (not groupPos)         
    my $slurp = "slurp -s 100M -o $tmpBam";
    my $outStream = "$toBam | $bamSort | $slurp";    
    open my $outH, "|-", $outStream or die "could not open samtools stream: $!\n";

    # process to sam with 10x-compatible CB:Z: tag
    my $molId = 1;
    print $outH $bamHeader;
    while (my $line = <$inH>) {
        chomp $line;
        my @f = split("\t", $line);
        print $outH join("\t",
            $molId, 
            _PROPER_PAIR_BITS,
            @f[_RNAME, _POS, _MAPQ],
            '100M', # samtools demands a CIGAR string, doesn't matter what it is, we will never use it
            '=',
            1,
            $f[_TLEN],
            '*',
            '*',
            "XC:i:$f[_DUP_RATE]",
            "CB:Z:$cellName"
        ), "\n";        
        $molId++;
    }
    
    # clean up
    close $inH;
    close $outH;
}
    
# report counts
printCount($nCells, 'nCells', 'number of cells aligned to genome');

# clean up
close $cellH;
close $bamListH;

