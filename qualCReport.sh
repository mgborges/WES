#run from raw
files=`ls ../qc/*/*_data* | awk -F " " '{print $1}'`
outpath='../qc/'

i=0

for file in $files;
do
if [ $i = 0 ]; then
   printf "Filename\tTotal Sequences\tPer base sequence quality\tPer sequence quality scores\tPer base sequence content\tGC Percent\tPer base GC content\tPer sequence GC content\tPer base N content\tSequence Length Distribution\tSequence Duplication Levels\tTotal Duplicate Percentage\tOverrepresented sequences\tKmer Content\tFile\n" >> $outpath/qualityCReport.txt
fi

grep "Filename" $file | awk -F "\t" '{print $2}' - >> tmp.out;

grep "Total Sequences" $file | awk -F "\t" '{print $2}' - >> tmp.out;

grep "Per base sequence quality" $file | awk -F "\t" '{print $2}' - >> tmp.out; 

grep "Per sequence quality scores" $file | awk -F "\t" '{print $2}' - >> tmp.out;

grep "Per base sequence content" $file | awk -F "\t" '{print $2}' - >> tmp.out;

grep "^%GC" $file | awk -F "\t" '{print $2}' - >> tmp.out;

grep "Per base GC content" $file | awk -F "\t" '{print $2}' - >> tmp.out;

grep "Per sequence GC content" $file | awk -F "\t" '{print $2}' - >> tmp.out;

grep "Per base N content" $file | awk -F "\t" '{print $2}' - >> tmp.out;

grep "Sequence Length Distribution" $file | awk -F "\t" '{print $2}' - >> tmp.out;

grep "Sequence Duplication Levels" $file | awk -F "\t" '{print $2}' - >> tmp.out;

grep "Total Duplicate Percentage" $file | awk -F "\t" '{print $2}' - >> tmp.out;

grep "Overrepresented sequences" $file | awk -F "\t" '{print $2}' - >> tmp.out;

grep "Kmer Content" $file | awk -F "\t" '{print $2}' - >> tmp.out;

grep "Filename" $file | awk -F "\t|.fastq.gz" '{OFS = ""} {print "=HIPERLINK(\"",$2,"_fastqc","\\","fastqc_report.html\"",";\"Link\")"}' >> tmp.out;


awk '
{ 
    for (i=1; i<=NF; i++)  {
a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
str=a[1,j]
for(i=2; i<=NR; i++){
    str=str"\t"a[i,j];
}
print str
    }
}' tmp.out >> $outpath/qualityCReport.txt;

rm tmp.out
i=1
done
