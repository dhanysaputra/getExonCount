samtools view -F 4 $1 | awk -v mapq=$2 '{ if ($5>=mapq && $7=="=") print $0 } ' | sort -k2,2 -t '.' > deleteme.sam
for i in {1..$3}; do grep -P "\tchr$i\t" deleteme.sam > $4.chr$i
grep -P '\tchrX\t' deleteme.sam > $4.chrX
grep -P '\tchrY\t' deleteme.sam > $4.chrY
rm -f deleteme.sam