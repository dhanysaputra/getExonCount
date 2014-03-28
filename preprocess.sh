module load bioinfo-tools
module load samtools/0.1.18

# Simple command: bash preprocess.sh file=yourfile.bam
# example: bash preprocess.sh srr28-29.bam 30 22 tumor

chr=""
file=""
parallel="-1"
for ARG in $*
do
        if [[ $ARG == file=* ]];
        then
                myparam=$( cut -d '=' -f 2- <<< "$ARG" )
                file=$myparam
        elif [[ $ARG == hseparator=* ]];
        then
                myparam=$( cut -d '=' -f 2- <<< "$ARG" )
                hseparator=$myparam
        elif [[ $ARG == hdirection=* ]];
        then
                myparam=$( cut -d '=' -f 2- <<< "$ARG" )
                hdirection=$myparam
        elif [[ $ARG == parallel=* ]];
        then
                myparam=$( cut -d '=' -f 2- <<< "$ARG" )
                parallel=$myparam
        elif [[ $ARG == batchsize=* ]];
        then
                myparam=$( cut -d '=' -f 2- <<< "$ARG" )
                batchsize=$myparam
        elif [[ $ARG == chr=*  ]];
        then
                myparam=$( cut -d '=' -f 2- <<< "$ARG" )
                if [ $myparam != "ALL" -a $myparam != "" ];
                then
                        chr=$myparam
                fi
        fi
done

if [[ $file == "" ]];
then
        echo 'Please provide the bamfile to process. file=yourfile.bam'
else:
        if [[ $parallel == "-1" ]];
        then
                echo 'parallel was not specified, 8 threads will be used.'
                parallel='8'
        fi
        date
        echo '1. SAM -> BAM'
        samtools view -F 4 $file | awk ' { if ( $7=="=" && $3 !~ /M/ ) print $0 } ' > $file.filtered
        date
        echo '2. Sort'
        /proj/b2012036/INBOX/Dhany/tryout/coreutils-8.13/src/sort --parallel=$parallel $file.filtered > $file.sorted;
        date
        echo '3. SortMA'
        python /proj/b2012036/INBOX/Dhany/tryout/sortMA.py $file.sorted > $file.sortMA;
        date
fi

for (( i=1; i<=22; i++ )); do awk -v j=$i '{ if($3==j) print $0 }' $file.sortMA > $file.chr$i & done
awk '{ if($3==X) print $0 }' accepted_hits.bam.sortMA > accepted_hits.chrX &
awk '{ if($3==Y) print $0 }' accepted_hits.bam.sortMA > accepted_hits.chrY &
wait
for (( i=1; i<=22; i++ )); do python getPairCount.py $file output.$i hg19.sqlite chr$i & done
python getPairCount.py $file output.X hg19.sqlite chrX &
python getPairCount.py $file output.Y hg19.sqlite chrY &
wait
cat output.* > finaloutput.txt
echo 'Your count can be read at finaloutput.txt'
