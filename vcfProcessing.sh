function vcfProcessing {
  
grep -v ^# $1 > tmp.fir.$1
grep -v ^# $2 > tmp.sec.$2

cut -f 8 tmp.sec.$2 | sed -e 's/EFF=/\t/g' -e 's/|/\t./g' -e 's/(/\t/g' - | awk '{print $2 "\t" $3 "\t" $4 "\t" $8}' | sed 's/\t\./\t/g' - > tmp.$VCFFile

paste tmp.fir.$1 tmp.$VCFFile | sed -e 's/\//|/g' > tmp5.$VCFFile

grep -n ^ tmp5.$VCFFile | awk '{print $1 "\t" $2 "\t" $4 "\t" $5}' > tmp5.$VCFFile.lines.sep

sed -e 's/:/\t/g' tmp5.$VCFFile.lines.sep > tmp5.$VCFFile.lines

uniq -f 1 tmp5.$VCFFile.lines > tmp5.$VCFFile.lines.uniq

cut -f 1 tmp5.$VCFFile.lines.uniq > tmp5.$VCFFile.lines.uniq.lines

while read line
do

line=`echo -e "$line"`
sed ''"$line"'q;d' tmp5.$VCFFile

done < tmp5.$VCFFile.lines.uniq.lines > $1.EFF.vcf

rm tmp*$VCFFile*
}
