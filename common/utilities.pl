use strict;
use warnings;

# utility functions common to linCNV scripts

# working variables
use vars qw($error $nCells);

#----------------------------------------------------------
# load the list of chromosomes in use
#----------------------------------------------------------
our (@chroms, %chroms, %chromNs, $nChroms);
sub loadChromList {
    my $faiFile = "$ENV{GENOME_DIR}/$ENV{GENOME}.fa.fai";
    open my $inH, "<", $faiFile or die "$error: could not open $faiFile for reading: $!\n";
    while(my $line = <$inH>){    
        chomp $line;
        my @f = split("\t", $line);
        $f[0] =~ m/_/ and next; # just use canonical chromosomes
        $f[0] =~ m/M$/ and next;
        $f[0] =~ m/EBV$/ and next;
        $ENV{INCLUDE_Y} or ($f[0] =~ m/Y$/ and next); # and possibly exclude Y
        $chroms{$f[0]} = \@f;
    }
    close $inH;
    my $chromN = 1;
    foreach my $chr(1..100, qw(X Y)){ # ensure that the chromosome list is properly sorted
        my $chrom = "chr$chr";
        $chroms{$chrom} or next;
        push @chroms, $chrom;
        $chromNs{$chrom} = $chromN;
        $chromN++;
    }
    $nChroms = @chroms;
}

#----------------------------------------------------------
# load sample-specific data resources
#----------------------------------------------------------
# cells previously designated as "true" by CellRanger
#   barcode,cell_id,total_num_reads,num_unmapped_reads,num_lowmapq_reads,num_duplicate_reads,num_mapped_dedup_reads,frac_mapped_duplicates,effective_depth_of_coverage,effective_reads_per_1Mbp,raw_mapd,normalized_mapd,raw_dimapd,normalized_dimapd,mean_ploidy,ploidy_confidence,is_high_dimapd,is_noisy,est_cnv_resolution_mb
#   AAACCTGCACCACACG-1,0,1401732,6661,222057,163321,1009693,0.11651371303501667,0.0700481200042838,370,0.11524713366463572,0.11524713366463572,1.0562023782600596,1.0562023782600596,1.949662924883774,-2,0,0,0.9989746093750012
# NB: cell ranger cell_id ranges from 0 to nCells-1, while our R-based numbers index cells from 1 to nCells
# the cell order is maintained, so adding or subtracting 1 to convert works
our (%cells, @allCells, $nCells, @acceptedCellIds, @isAccepted, $nAcceptedCells);
sub loadCells10X {
    my $cellsFile = "$ENV{CELL_RANGER_DIR}/per_cell_summary_metrics.csv";
    open my $inH, "<", $cellsFile or die "$error: could not open $cellsFile for reading: $!\n";
    my $header = <$inH>;
    while(my $line = <$inH>){
        my ($CB, $cellId) = split(",", $line);
        $cells{$CB} = $cellId;
        push @allCells, $cellId;
    }
    close $inH;
    $nCells = @allCells;      
    my $binWeightAllCells = $nCells * $ENV{WEIGHT_PER_CELL};
    return $binWeightAllCells;
}
sub loadAcceptedCells {
    @acceptedCellIds = split(' ', $ENV{ACCEPTED_CELLS});
    $nAcceptedCells = @acceptedCellIds;
    @isAccepted = (0) x $nCells; # a boolean lookup
    map { $isAccepted[$_] = 1 } @acceptedCellIds;
    my $binWeightAllCells = $nAcceptedCells * $ENV{WEIGHT_PER_CELL}; 
    return $binWeightAllCells;  
}

1;

