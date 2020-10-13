use strict;
use warnings;

# constants
use constant {
	CHROM => 0, # BED fields
	START => 1,
	END_ => 2,
	NAME => 3,
	SCORE => 4,
	STRAND => 5,
	IS_BAD_REGION => 6
};

# working variables
my (%counts, %cells);
my $binSize = $ENV{BIN_SIZE};

# create a crosstab of genome bin x cell, values are read pair counts
while(my $line = <STDIN>){
	chomp $line;
	my @f = split("\t", $line);
	
	# reject read pairs in bad genome regions
	$f[IS_BAD_REGION] and next; 
	
	my $bin = int($f[START] / $binSize) * $binSize;	
	$counts{$f[CHROM]}{$bin}{$f[NAME]}++;
	$counts{$f[CHROM]}{$bin}{ALL}++;
	$cells{$f[NAME]}++;
}

# initialize the ouput
my @cells = sort {$a cmp $b} keys %cells;
print join("\t", qw(chrom start end binN all_cells strand), @cells), "\n";

# commit the crosstab
my $binN = 1;
foreach my $chrom(sort {$a cmp $b} keys %counts){
	foreach my $bin(sort {$a <=> $b} keys %{$counts{$chrom}}){
		print join("\t", $chrom, $bin, $bin + 1e6, $binN, $counts{$chrom}{$bin}{ALL}, '.');
		foreach my $cell(@cells){
			my $count = $counts{$chrom}{$bin}{$cell} || 0;
			print "\t$count";
		}
		print "\n";
		$binN++;	
	}	
}

