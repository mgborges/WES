#######################################
#### Global variables (may be adjusted)

FRAME_LIST=(R1 R2)
REF="/home/bioinf/ref/GRCh38.p1/wholeGenome/genome.fa"
DBSNP="/home/bioinf/ref/var/All.vcf"
GATK="java -jar -Xmx20G /home/bioinf/exoma/ell/bin/GenomeAnalysisTK.jar"
PICARD="java -jar -Xmx3G /home/bioinf/exoma/ell/bin/picard-tools-1.124/picard.jar"
THREADS="10"

#######################################
#### Functions

function mergegz {
## Function to merge separated .fastq.gz files from HiSeq 2500

f=$1
cd $f
filename=`ls *gz | awk -F "_R" '{print $1}' | uniq`
for file in $filename
do
zcat $file*R1* | gzip > ../$file.R1.fastq.gz
zcat $file*R2* | gzip > ../$file.R2.fastq.gz
done
cd ..
}

function aln {
# Function implements alignment against a genome with Burrows-Wheeler algorithm
sample=$1
fasta1="${sample}.${FRAME_LIST[0]}.fastq.gz"
fasta2="${sample}.${FRAME_LIST[1]}.fastq.gz"
bam="../aln/${sample}.bam"

bwa mem -t $THREADS $REF $fasta1 $fasta2 > $bam.mem 2> $bam.mem.log
}


function post_aln {
# Implements intemediate .bam processing steps with GATK
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
    
    $PICARD SortSam I=$bam.mem O=$bam SO=coordinate &&
    rm $bam.mem &&
    $PICARD MarkDuplicates INPUT=$bam OUTPUT=/dev/stdout METRICS_FILE=$bam.metrics | $PICARD AddOrReplaceReadGroups I=/dev/stdin O=$bam.tmp RGID=$id RGSM=$SM RGCN=LaCTAD RGPL=ILLUMINA RGPU=$FLOWCELL-$BARCODE.$LANE RGLB=Nextera 2> $bam.AddOrReplaceReadGroups.log &&
    rm $bam; mv $bam.tmp $bam &&
    ## build index and validate to GATK
    $PICARD BuildBamIndex INPUT=$bam OUTPUT=$bam.bai 2> $bam.bai.log &&
    $PICARD ValidateSamFile INPUT=$bam OUTPUT=$bam.vald VALIDATE_INDEX=true MODE=SUMMARY 2> $bam.vald.log &&
    ## Indel realignment
    $GATK -T RealignerTargetCreator -R $REF -I $bam -o $bam".intervals" -known $DBSNP -log $bam".intervals.log" &&
    $GATK -T IndelRealigner -R $REF -I $bam -targetIntervals $bam".intervals" -o $bam".realn.bam" -known $DBSNP -LOD 0.4 -model USE_READS -compress 0 --disable_bam_indexing -log $bam".realn.bam.log" &&
    rm $bam &&
    $PICARD BuildBamIndex INPUT=$bam".realn.bam" OUTPUT=$bam".realn.bam.bai" 2> $bam".realn.bam.log"&&
    ## Quality recalibration
    $GATK -T BaseRecalibrator -R $REF -I $bam".realn.bam" -o $bam".recal.csv" -knownSites $DBSNP -l INFO -log $bam".recal.csv.log" &&
    $GATK -T PrintReads -R $REF -I $bam".realn.bam" -BQSR $bam".recal.csv" -o $bam".recal.bam" -log $bam".recal.bam.log" &&
    
    rm $bam".realn.bam" &&
    rm $bam".realn.bam.bai" &&
    rm $bam.bai &&
    mv $bam".recal.bam" $bam &&
    mv $bam".recal.bai" "../aln/${sample}.bai" &&
    ## post alignment qc
    $PICARD CollectMultipleMetrics INPUT=$bam REFERENCE_SEQUENCE=$REF OUTPUT="../qc/"$bam PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics  PROGRAM=QualityScoreDistribution  PROGRAM=MeanQualityByCycle  PROGRAM=CollectBaseDistributionByCycle 2> "../qc/"$bam.CollectMultipleMetrics.log
    $PICARD CalculateHsMetrics I=$bam O=$bam.HsMetrics  R=$REF TARGET_INTERVALS=/home/bioinf/ref/GRCh38.p1/nexteraRapidCaputureRegions.bed  BAIT_INTERVALS=/home/bioinf/ref/GRCh38.p1/nexteraRapidCaputureRegions.bed
}

function mergeBAM {
## Used in case of multiplexing (e.g. a same sample sequenced in multiple lanes or flowcells)
sample=$1
cd ../raw
files=`ls $sample*.gz | awk -F "." '{print "INPUT=" $1 ".bam"}' | uniq`
bam=`ls $sample* | awk -F "_" '{print $1 ".bam"}' | uniq`

cd ../aln

$PICARD MergeSamFiles $files OUTPUT=$bam 2> $bam.log &&

filesToRemove=`echo $files | sed 's/INPUT=//g'` &&
rm $filesToRemove &&
rm *bai &&

$PICARD MarkDuplicates INPUT=$bam OUTPUT=$bam.tmp METRICS_FILE=$bam.metrics 2> $bam.duplic.log &&
mv $bam.tmp $bam &&

$PICARD BuildBamIndex INPUT=$bam OUTPUT=$bam.bai 2> $bam.bai.log
$GATK -T RealignerTargetCreator -R $REF -I $bam -o $bam".intervals" -known $DBSNP -log $bam".intervals.log"
$GATK -T IndelRealigner -R $REF -I $bam -targetIntervals $bam".intervals" -o $bam".realn.bam" -known $DBSNP -LOD 0.4 -model USE_READS -compress 0 --disable_bam_indexing -log $bam".realn.bam.log" --maxReadsForRealignment 40000 &&
rm $bam $bam.bai
$PICARD BuildBamIndex INPUT=$bam".realn.bam" OUTPUT=$bam".realn.bam.bai" 2> $bam".realn.bam.log"
mv $bam".realn.bam" $bam
mv $bam".realn.bam.bai" "../aln/${sample}.bai"
$PICARD CollectMultipleMetrics INPUT=$bam REFERENCE_SEQUENCE=$REF OUTPUT="../qc/"$bam PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics  PROGRAM=QualityScoreDistribution  PROGRAM=MeanQualityByCycle  PROGRAM=CollectBaseDistributionByCycle 2> "../qc/"$bam.CollectMultipleMetrics.log
}

##MERGE GZ###################

SAMPLE_FOLDER=`ls`
for i in $SAMPLE_FOLDER
do
mergegz $s
done

### ALN #####################

cd ../raw
SAMPLE_LIST=`ls -Sr *R1* | awk -F "." '{print $1}'`
for i in $SAMPLE_LIST
do
aln $i
done

### POST ALN ################

cd ../raw
SAMPLE_LIST=`ls -Sr *R1* | awk -F "." '{print $1}'`
for i in $SAMPLE_LIST
post_aln $s;
done

### MULTIPLEXING ############

cd ../raw
samples=`ls -rS *R1*gz | awk -F "_" '{print $1}' | uniq`
for i in $samples
mergeBAM $s;
done
### VARITANT CALLING ########

cd ../aln
bams=`ls *.bam`
for bam in $bams
do
$GATK -T HaplotypeCaller -R $REF -I $bam -o ../results/$bam.g.vcf -log ../results/$bam.g.vcf.log -D $DBSNP -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000
done

### MERGE CALLINGS ###########
### Uncomment and adjust the names within "<>"

# $GATK -T GenotypeGVCFs -R $REF -o ../results/<STUDY NAME>.g.raw.vcf -V <FILE1> -V <FILE2> -D $DBSNP 
# cd ../results
# /home/bioinf/exomapqm/pqmAnaPaula/bin/ensembl-tools-release-79/scripts/variant_effect_predictor/variant_effect_predictor.pl -i <STUDY NAME>.g.raw.vcf -o <STUDY NAME>.g.VEP.vcf --format vcf --fasta $REF --sift b --polyphen b --regulatory --symbol --gmaf --vcf --cache 

