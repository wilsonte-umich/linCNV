use strict;
use warnings;

#----------------------------------------------------------------
# align a series of individual cells to genome and write output to a single stream
#----------------------------------------------------------------

# initialize reporting
our $script = 'align';
our $error  = "$script error";
my ($nCells, $nAlignments) = (0) x 10;

# load dependencies
use File::Basename;
my $scriptDir = dirname(__FILE__);
require "$scriptDir/../_workflow/workflow.pl";
require "$scriptDir/../common/utilities.pl";
resetCountFile();

# initialize outputs
open my $outH,  ">", $ENV{NAMED_PIPE}     or die "could not open $ENV{NAMED_PIPE} for writing\n";
open my $cellH, ">", $ENV{CELL_LIST_FILE} or die "could not open $ENV{CELL_LIST_FILE} for writing\n";
print $cellH join(",", qw(barcode cell_id)), "\n";

# get inputs
my @cellDirs = split(/\s+/, $ENV{CELL_DIRS});

# align one cell at a time with parallelization
foreach my $cellName(@cellDirs){
    
    # add cell to cell list in abbreviated CellRanger format    
    $nCells++;
    my $cellN = $nCells - 1;
    print STDERR "  $cellName (cell #$cellN)\n";
    print $cellH join(",", $cellName, $cellN), "\n";
    
    # align using BWA    
    my $fastQs = join(" ", glob("$cellName/*.fastq.gz"));
    my $bwa = "bwa mem -t $ENV{N_CPU} $ENV{GENOME_DIR}/$ENV{GENOME}.fa $fastQs 2>$ENV{LOGS_PREFIX}.align.bwa.log";
    open my $inH, "-|", $bwa or die "could not open bwa stream\n";
    
    # add the cell name as a pretend barcode in 10x scCNV format
    my $printHeader = ($nCells == 1);    
    while (my $line = <$inH>) {
        $nAlignments++;
        if ($line =~ m/^\@/) {
            $printHeader and print $outH $line; # only print the header once, at the first sample
        } else {
            chomp $line;
            print $outH "$line\tCB:Z:$cellName\n";
        }  
    }
    close $inH;
}

# report counts
printCount($nCells,      'nCells',      'number of cells aligned to genome');
printCount($nAlignments, 'nAlignments', 'number of aligned segments over all cells');

# clean up
close $cellH;
close $outH;

