
# thread through sorted, high-quality read pairs and parse into bins
# output is one file per chromosome
perl $PIPELINE_DIR/bin/bin.pl
checkPipe

# merge chromosome bin files
cat \
<(zcat $BIN_PREFIX.counts.*.gz | head -n1) \
<(zcat $BIN_PREFIX.counts.*.gz | grep -v start) |
pigz -p $N_CPU -c |
$SLURP -s 100M -o $BIN_PREFIX.counts.gz
checkPipe

# determine the number of excluded bases per bin
cat \
<(echo -e "#chrom\tstart\tend\texcluded") \
<(zcat $BIN_PREFIX.counts.gz |
  cut -f 1-3 |
  grep -v start |
  bedtools intersect -wao -a stdin -b $EXCLUSIONS_FILE |
  bedtools groupby -g 1,2,3 -c 8 -o sum) | 
pigz -p $N_CPU -c |
$SLURP -s 100M -o $BIN_PREFIX.exclusions.gz
checkPipe

# establish mappability values for each variable bin
# to improve mappability reference, could see https://www.biostars.org/p/181014/
cat \
<(echo -e "#chrom\tstart\tend\texcluded\tmappability") \
<(bedtools intersect -v -a $MAPPABILITY_FILE -b $EXCLUSIONS_FILE |
  bedtools intersect -wao -a $BIN_PREFIX.exclusions.gz -b stdin |
  bedtools groupby -g 1,2,3,4 -c 9 -o mean) | 
pigz -p $N_CPU -c |
$SLURP -s 100M -o $BIN_PREFIX.mappability.gz # has both excluded base count and ~mappability fraction of non-excluded bases
checkPipe

# establish GC fractions for each variable bin
cat \
<(echo -e "#chrom\tstart\tend\texcluded\tmappability\tgc_fraction") \
<(bedtools intersect -v -a $GC_FILE -b $EXCLUSIONS_FILE |
  bedtools intersect -wao -a $BIN_PREFIX.mappability.gz -b stdin |
  bedtools groupby -g 1,2,3,4,5 -c 10 -o mean) | 
pigz -p $N_CPU -c |
$SLURP -s 100M -o $BIN_PREFIX.gc_fraction.gz # has excluded, ~mappability and GC fraction
checkPipe

# clean up
rm $BIN_PREFIX.exclusions.gz
rm $BIN_PREFIX.mappability.gz

# process the bin x cell matrix in R, adjusting for ploidy and CNVs
Rscript $PIPELINE_DIR/bin/adjust.R
checkPipe

