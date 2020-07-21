use strict;
use warnings;

# this script is not used by the pipeline itself
# it is a helper script to assist with construction of the required 1 kbp bin GC bed file

# example
# zcat hg38.gc5Base.wigVarStep.gz | perl /treehouse/wilsonte_lab/ssd/pipelines/wilson/action_scripts/executables/linCNV/bin/gc5base2Bed.pl | gzip -c > hg38.gc5Base.bin_1000.bed.gz

my ($span, $chrom, $start, $end);
my ($binSum, $binN) = (0, 0);
my $outSpan = 1000;

while(my $line = <STDIN>){
    chomp $line;
    if($line =~ m/variableStep chrom=(\w+) span=(\d+)/){
        printBin(undef);
        $chrom = $1;
        print STDERR "$chrom\n";
        $span = $2;
        $span != 5 and die "not all spans are 5!\n";
    } else {
        my ($pos, $gc) = split("\t", $line);
        if(defined $start){
            $pos > $end and printBin($pos);
        } else {
            setBin($pos);
        }
        $binSum += $gc / 100;
        $binN++;
    }   
}

sub printBin {
    my ($pos) = @_;
    defined $start and print join("\t", $chrom, $start, $end, $binN, int($binSum / $binN * 1000 + 0.5) / 1000, "."), "\n";
    setBin($pos);
}
sub setBin {
    my ($pos) = @_;
    ($binSum, $binN) = (0, 0);
    if($pos){
        $start = $pos - 1;
        $end = $start + $outSpan;          
    } else {
        $start = undef;
    }
}

#variableStep chrom=chr1 span=5
#10001   40
#10006   40
#10011   40
#10016   60
#10021   60
#10026   60
#10031   40
#10036   40
#10041   40


