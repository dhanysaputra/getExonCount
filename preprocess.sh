samtools view -F 4 <yourfile.bam> | awk '{ if ($5>=30 && $7=="=") print $0 } ' | sort -k2,2 -t '.' > deleteme.sam
grep -P '\tchr1\t' deleteme.sam > human1.chr1
grep -P '\tchr2\t' deleteme.sam > human1.chr2
grep -P '\tchr3\t' deleteme.sam > human1.chr3
grep -P '\tchr4\t' deleteme.sam > human1.chr4
grep -P '\tchr5\t' deleteme.sam > human1.chr5
grep -P '\tchr6\t' deleteme.sam > human1.chr6
grep -P '\tchr7\t' deleteme.sam > human1.chr7
grep -P '\tchr8\t' deleteme.sam > human1.chr8
grep -P '\tchr9\t' deleteme.sam > human1.chr9
grep -P '\tchr10\t' deleteme.sam > human1.chr10
grep -P '\tchr11\t' deleteme.sam > human1.chr11
grep -P '\tchr12\t' deleteme.sam > human1.chr12
grep -P '\tchr13\t' deleteme.sam > human1.chr13
grep -P '\tchr14\t' deleteme.sam > human1.chr14
grep -P '\tchr15\t' deleteme.sam > human1.chr15
grep -P '\tchr16\t' deleteme.sam > human1.chr16
grep -P '\tchr17\t' deleteme.sam > human1.chr17
grep -P '\tchr18\t' deleteme.sam > human1.chr18
grep -P '\tchr19\t' deleteme.sam > human1.chr19
grep -P '\tchr20\t' deleteme.sam > human1.chr20
grep -P '\tchr21\t' deleteme.sam > human1.chr21
grep -P '\tchr22\t' deleteme.sam > human1.chr22
grep -P '\tchrX\t' deleteme.sam > human1.chrX
grep -P '\tchrY\t' deleteme.sam > human1.chrY
rm -f deleteme.sam