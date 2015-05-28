function mergegz
{
folders=`ls`
for f in $folders
do
{
cd $f
filename=`ls *gz | awk -F "_R" '{print $1}' | uniq`
for file in $filename
do
zcat $file*R1* | gzip > ../$file.R1.fastq.gz
zcat $file*R2* | gzip > ../$file.R2.fastq.gz
done
cd ..
} &
done
}



# Global variables
FRAME_LIST=(R1 R2)
REF="/home/bioinf/ref/GRCh38.p1/wholeGenome/genome.fa"
DBSNP="/home/bioinf/ref/var/All.vcf"
GATK="java -jar -Xmx20G /home/bioinf/exoma/ell/bin/GenomeAnalysisTK.jar"
PICARD="java -jar -Xmx3G /home/bioinf/exoma/ell/bin/picard-tools-1.124/picard.jar"
THREADS="12"

## samples to be analysed
SAMPLE_LIST=`ls *R1* | awk -F "." '{print $1}'`

## Creation of folders
log="../aln/log/"
mkdir ../qc ../aln ../results $log

function aln {
    sample=$1

    ## input fastq files
    fasta1="${sample}.${FRAME_LIST[0]}.fastq.gz"
    fasta2="${sample}.${FRAME_LIST[1]}.fastq.gz"
    bam="../aln/${sample}.bam";
    
    ## id from fasta files
    id=$(zcat ${fasta1} | head -n 1 | grep "@" | cut -f3-4 -d: | uniq)
    SM=`echo $sample | awk -F _ '{print $1}'`
    FLOWCELL=`zcat ${fasta1} | head -n 1 | grep "@" | cut -f 3 -d:`
    BARCODE=`zcat ${fasta1} | head -n 1 | grep "@" | cut -f 10 -d:`
    LANE=`zcat ${fasta1} | head -n 1 | grep "@" | cut -f 4 -d:`
    ## alignment
    bwa mem -t $THREADS $REF $fasta1 $fasta2 | $PICARD SortSam I=/dev/stdin O=$bam SO=coordinate &&
    ## markdup
    $PICARD MarkDuplicates INPUT=$bam OUTPUT=/dev/stdout METRICS_FILE=$bam.metrics | $PICARD AddOrReplaceReadGroups I=/dev/stdin O=$bam.tmp RGID=$id RGSM=$SM RGCN=LaCTAD RGPL=ILLUMINA RGPU=$FLOWCELL-$BARCODE.$LANE RGLB=Nextera &&
    rm $bam; mv $bam.tmp $bam &&
    ## build index and validate to GATK
    $PICARD BuildBamIndex INPUT=$bam OUTPUT=$bam.bai &&
    $PICARD ValidateSamFile INPUT=$bam OUTPUT=$bam.vald VALIDATE_INDEX=true MODE=SUMMARY &&

    $GATK -T RealignerTargetCreator -nt $THREADS -R $REF -I $bam -o $bam".intervals" -known $DBSNP -log $bam".intervals.log" &&
    $GATK -T IndelRealigner -R $REF -I $bam -targetIntervals $bam".intervals" -o $bam".realn.bam" -known $DBSNP -LOD 0.4 -model USE_READS -compress 0 --disable_bam_indexing -log $bam".realn.bam.log" &&
    rm $bam &&
    $PICARD BuildBamIndex INPUT=$bam".realn.bam" OUTPUT=$bam".realn.bam.bai" &&
    $GATK -T BaseRecalibrator -nct $THREADS -R $REF -I $bam".realn.bam" -o $bam".recal.csv" -knownSites $DBSNP -l INFO -log $bam".recal.csv.log" &&
    $GATK -T PrintReads -nct $THREADS -R $REF -I $bam".realn.bam" -BQSR $bam".recal.csv" -o $bam".recal.bam" -log $bam".recal.bam.log" &&
    
    rm $bam".realn.bam" &&
    rm $bam".realn.bam.bai" &&
    rm $bam.bai &&
    mv $bam".recal.bam" $bam &&
    mv $bam".recal.bai" "../aln/${sample}.bai"
    ## post alignment qc
    $PICARD CollectMultipleMetrics INPUT=$bam REFERENCE_SEQUENCE=$REF OUTPUT="../qc/"$bam PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics  PROGRAM=QualityScoreDistribution  PROGRAM=MeanQualityByCycle  PROGRAM=CollectBaseDistributionByCycle 
}

function mergeBAM {
cd ../aln &&
samples=`ls *bam | awk -F "_" '{print $1}' | uniq`
for sample in $samples
do
files=`ls $sample* | awk -F "." '{print "INPUT=" $1 ".bam"}' | uniq`
bam=`ls $sample* | awk -F "_" '{print $1 ".bam"}' | uniq`

cd ../aln &&

$PICARD MergeSamFiles $files OUTPUT=$bam &&
#    $PICARD ValidateSamFile INPUT=$bam OUTPUT=$bam.vald VALIDATE_INDEX=true MODE=SUMMARY &&
$PICARD MarkDuplicates INPUT=$bam OUTPUT=$bam.tmp METRICS_FILE=$bam.metrics &&
rm $bam && mv $bam.tmp $bam &&
$PICARD BuildBamIndex INPUT=$bam OUTPUT=$bam.bai &&
$GATK -T RealignerTargetCreator -nt $THREADS -R $REF -I $bam -o $bam".intervals" -known $DBSNP -log $bam".intervals.log" &&
$GATK -T IndelRealigner -R $REF -I $bam -targetIntervals $bam".intervals" -o $bam".realn.bam" -known $DBSNP -LOD 0.4 -model USE_READS -compress 0 --disable_bam_indexing -log $bam".realn.bam.log" --maxReadsForRealignment 40000 &&
rm $bam $bam.bai &&
$PICARD BuildBamIndex INPUT=$bam".realn.bam" OUTPUT=$bam".realn.bam.bai" &&
mv $bam".realn.bam" $bam &&
mv $bam".realn.bam.bai" "../aln/${sample}.bai" &&
$PICARD CollectMultipleMetrics INPUT=$bam REFERENCE_SEQUENCE=$REF OUTPUT="../qc/"$bam PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics  PROGRAM=QualityScoreDistribution  PROGRAM=MeanQualityByCycle  PROGRAM=CollectBaseDistributionByCycle 
done
}

## QC pre aln
fastqc -o ../qc/ -t $THREADS *.fastq.gz &&

for sample in $SAMPLE_LIST
do 
    aln $sample
done

cd ../aln

bams=`ls *bam`

for bam in $bams
do 
$PICARD CalculateHsMetrics I=$bam O=$bam.HsMetrics  R=/home/bioinf/ref/GRCh38.p1/wholeGenome/genome.fa TARGET_INTERVALS=/home/bioinf/ref/GRCh38.p1/nexteraRapidCaputureRegions.bed  BAIT_INTERVALS=/home/bioinf/ref/GRCh38.p1/nexteraRapidCaputureRegions.bed   
done

mv * log/
mv log/*bai .
mv log/*bam .

mergeBAM

mv * log/
mv log/*bai .
mv log/*bam .
#aumentar memória se necessário
#$GATK -T HaplotypeCaller -R $REF -I monica1.bam -I monica2.bam -I monica3.bam -I monica4.bam -o ../results/monica.raw.vcf -stand_call_conf 40 -stand_emit_conf 10 -minPruning 3 -log ../results/monica.raw.log -D $DBSNP -maxAltAlleles 18

#perl /home/bioinf/exoma/eltm_38/bin/ensembl-tools-release-77/scripts/variant_effect_predictor/variant_effect_predictor.pl -i ../results/monica.raw.vcf -o ../results/monica.VEP.vcf --cache --format vcf --fasta /home/bioinf/ref/GRCh38.p1/wholeGenome/genome.fa --sift b --polyphen b --regulatory --symbol --gmaf --pubmed --vcf
