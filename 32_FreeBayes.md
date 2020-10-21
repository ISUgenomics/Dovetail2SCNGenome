# Need to call variants from the different populations to see how some may influence protein structure
```
Attempting to run a freebayes pipeline for the 15 populations  from the draft genome paper, the TN10 illumina reads provided by Matt Hudsons lab, and whatever other sequencing I can find.   
```

### File setup
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/43_FreeBayes


# All of the reads, renamed to something interpretable
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/G3/FCC7GPJANXX-CHKPEI15080015_L1_1.fq.gz G3_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/G3/FCC7GPJANXX-CHKPEI15080015_L1_2.fq.gz G3_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/LY1/LY1_S6_L006_R1_001.fastq.gz LY1_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/LY1/LY1_S6_L006_R2_001.fastq.gz LY1_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/OP20/OP20_S5_L006_R1_001.fastq.gz OP20_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/OP20/OP20_S5_L006_R2_001.fastq.gz OP20_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/OP25/FCC7GPJANXX-CHKPEI15080011_L1_1.fq.gz OP25_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/OP25/FCC7GPJANXX-CHKPEI15080011_L1_2.fq.gz OP25_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/OP50/FCC7GPJANXX-CHKPEI15080012_L1_1.fq.gz OP50_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/OP50/FCC7GPJANXX-CHKPEI15080012_L1_2.fq.gz OP50_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/PA3/FCC7GPJANXX-CHKPEI15080014_L1_1.fq.gz PA3_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/PA3/FCC7GPJANXX-CHKPEI15080014_L1_2.fq.gz PA3_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN13/FCC7GPJANXX-wHAIPI022664-37_L1_1.fq.gz TN13_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN13/FCC7GPJANXX-wHAIPI022664-37_L1_2.fq.gz TN13_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN15/TN15_S2_L006_R1_001.fastq.gz TN15_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN15/TN15_S2_L006_R2_001.fastq.gz TN15_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN16/FCC7GPJANXX-wHAIPI022665-46_L1_1.fq.gz TN16_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN16/FCC7GPJANXX-wHAIPI022665-46_L1_2.fq.gz TN16_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN19/FCC7GPJANXX-CHKPEI15080013_L1_1.fq.gz TN19_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN19/FCC7GPJANXX-CHKPEI15080013_L1_2.fq.gz TN19_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN1/FCC7GPJANXX-wHAIPI022662-14_L1_1.fq.gz TN1_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN1/FCC7GPJANXX-wHAIPI022662-14_L1_2.fq.gz TN1_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN21/TN21_S3_L006_R1_001.fastq.gz TN21_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN21/TN21_S3_L006_R2_001.fastq.gz TN21_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN22/TN22_S4_L006_R1_001.fastq.gz TN22_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN22/TN22_S4_L006_R2_001.fastq.gz TN22_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN7/FCC7GPJANXX-wHAIPI022663-35_L1_1.fq.gz TN7_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN7/FCC7GPJANXX-wHAIPI022663-35_L1_2.fq.gz TN7_R2_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN8/TN8_S1_L006_R1_001.fastq.gz TN8_R1_fastq.gz
ln -s ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN8/TN8_S1_L006_R2_001.fastq.gz TN8_R2_fastq.gz
ln -s ../07_TN10IlluminaReads/27_27799_CAGATC_L002_R2_001.fastq.gz TN10_R2_fastq.gz
ln -s ../07_TN10IlluminaReads/27_27799_CAGATC_L002_R1_001.fastq.gz TN10_R1_fastq.gz

#The genome with the names for publication
ln -s ../49_RenameChromosomes/01_Transfer2Box/SCNgenome.fasta

```


### Prepare Gsnap alignment, as hisat2 did not do such a great time last time.
```
#build the genome database
ml gmap-gsnap/2018-07-04-gtu46xu
gmap_build -d SCNgenome -D SCNgenome SCNgenome.fasta


#run_gsnapPE.sh
##############################################################################################################################################
#!/bin/bash

module load gmap-gsnap/2018-07-04-gtu46xu
dbname=$1
dbloc=$2
dbfasta=$3
read1=$4
read2=$5
#gmap_build -d $dbname  -D $dbloc $dbfasta
gsnap -D $dbloc -d $dbname  --orientation=RF --pairmax-dna=10000 --pairdev=500 --gunzip -B 5 -t 16 -N 1 --input-buffer-size=1000000 --output-buffer-size=1000000 -A sam  $read1 $read2 >${dbname%.*}.${read1%.*}.sam
#############################################################################################################################################


# Build the script
paste <(ls *R1_fastq.gz ) <(ls *R2_fastq.gz) |while read line; do echo "sh run_gsnapPE.sh SCNgenome /work/GIF/remkv6/Baum/04_Dovetail2Restart/43_FreeBayes/SCNgenome/ SCNgenome.fasta " $line; done >align.sh

# I actually cut the tn10 read alignment short, as it was running significantly longer with 2-3 times the largest file size in the group
#just clipping the last line in the sam file.
mv SCNgenome.TN10_R1_fastq.sam UnfinishedSCNgenome.TN10_R1_fastq.sam
head -n -1 UnfinishedSCNgenome.TN10_R1_fastq.sam >SCNgenome.TN10_R1_fastq.sam &


# Sort, convert to bam
for f in *sam; do echo "ml samtools; samtools view --threads 16 -b "$f"|samtools sort -o "${f%.*}"_sorted.bam";done >samtools.sh

# get alignment stats
for f in *sam; do echo "ml samtools; samtools flagstat "$f" >"${f%.*}"stats"; done >stats.sh

for f in *stats; do awk '{if(NR==1){print FILENAME,$1}else if(NR==5) {print FILENAME,$1,$5}}' $f;done |sed 's/(//g' |tr "\n" " " |sed 's/%/%\n/g' |awk '{print $1,$4,$2,$5}' |tr " " "\t" |less
#######################################################################################################
SCNgenome.G3_R1_fastqstats      548128146       549172676       99.81%
SCNgenome.LY1_R1_fastqstats     782108051       791927572       98.76%
SCNgenome.OP20_R1_fastqstats    457081693       473106024       96.61%
SCNgenome.OP25_R1_fastqstats    184538022       186233955       99.09%
SCNgenome.OP50_R1_fastqstats    34072149        34857042        97.75%
SCNgenome.PA3_R1_fastqstats     742460927       743852256       99.81%
SCNgenome.TN10_R1_fastqstats    11289582        11417284        98.88%
SCNgenome.TN13_R1_fastqstats    335808379       336772437       99.71%
SCNgenome.TN15_R1_fastqstats    1258134071      1262720552      99.64%
SCNgenome.TN16_R1_fastqstats    5266615 5356162 98.33%
SCNgenome.TN19_R1_fastqstats    274092409       275638595       99.44%
SCNgenome.TN1_R1_fastqstats     451279049       456705114       98.81%
SCNgenome.TN21_R1_fastqstats    895577919       903894102       99.08%
SCNgenome.TN22_R1_fastqstats    1768922367      1774394578      99.69%
SCNgenome.TN7_R1_fastqstats     244247401       248253921       98.39%
SCNgenome.TN8_R1_fastqstats     1980377810      1987499308      99.64%
#######################################################################################################


```

### Assign read groups
```

for f in *_sorted.bam; do echo "sh run_PicardReadGroups.sh "$f;done >addReadGroups.sh


#check to make sure picard script will create proper names

#run_PicardReadGroups.sh
##################################################################################
#!/bin/bash
ulimit -c unlimited

BAM="$1"
ml picard


RGID=$(basename $BAM |sed 's/_R1/\t/1' |cut -f 1)
RGSM=$(basename $BAM| sed 's/_R1/\t/1' |cut -f 1 |awk '{print $0"SM"}'  )
RGLB="${RGSM}-L001"
RGPU=001
RGPL=ILLUMINA
echo -e "$RGID\t$RGSM\t$RGLB\t$RGPU"
java -Djava.io.tmpdir=$TMPDIR -Xmx180G -jar $PICARD/picard.jar AddOrReplaceReadGroups \
      I=${BAM} \
      O=${BAM%.*}_new.bam \
      RGID=$RGSM \
      RGLB=$RGLB \
      RGPL=$RGPL \
      RGPU=$RGPU \
      RGSM=$RGSM
module load samtools
samtools index ${bam%.*}_new.bam
##################################################################################

```

### Remove duplicate reads
```

for f in *bam; do echo "ml picard; ml jdk; java -jar $PICARD MarkDuplicates I="$f" O="${f%.*}"Dedup.bam M="${f%.*}"Dedup.metrics.txt REMOVE_SEQUENCING_DUPLICATES=true";done >dedup.sh
```
