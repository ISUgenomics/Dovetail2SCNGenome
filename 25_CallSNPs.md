### Would like to find differences in the genes within populations

### Gather data
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/43_FreeBayes

for f in ../../../../archive1/Baum/101415_SNPSCN/RawReads/clean_reads/*/*gz; do ls $f;done |awk '{print $0"\t"$0}'|sed 's|/|\t|10' |sed 's|/|\t|9' |cut -f 2,4,6 |awk '{print "ln -s "$2,$1" _fastq.gz"}' |awk '{print NR%2,$0}'|awk '{if($1==1){print $2,$3,$4,$5"_R1"$6} else {print $2,$3,$4,$5"_R2"$6}}'  >SOFTLINKER.SH

# matt's reads, 300bp paired end.  (not 150bp)
ln -s ../07_TN10IlluminaReads/27_27799_CAGATC_L002_R2_001.fastq.gz TN10_R2_fastq.gz
ln -s ../07_TN10IlluminaReads/27_27799_CAGATC_L002_R1_001.fastq.gz TN10_R1_fastq.gz

ln -s ../10_RepeatModeler/SCNgenome.fasta

hisat2-build SCNgenome SCNgenome.fasta


paste <(ls -1 *R1_fastq.gz ) <(ls -1 *R2_fastq.gz) |awk '{print "sh runFeatureCountsNBigwig.sh "$1,$2,"/work/GIF/remkv6/Baum/04_Dovetail2Restart/43_FreeBayes SCNgenome.fasta" }' >align.sh


################################################################################################################
#runFeatureCounts.sh
#note: I premade the database so it would not be made again and again
###############################################################################
#!/bin/bash

#Working model.  Things to consider: the script needs different settings for 2 analyses. dna-seq (hisat2), as this is set up for RNA.  free node settings vs regular node settings (procand memory).
#You must provide the following. Note variable DBDIR does not need a "/" at the end.
# sh runFeatureCounts.sh 28 sequence_1.fastq sequence_2.fastq /work/GIF/remkv6/files genome.fa genome.gff3


PROC=16
R1_FQ="$1"
R2_FQ="$2"
DBDIR="$3"
GENOME="$4"
GFF="$5"

module load trimgalore
module load py-setuptools/35.0.2-py2-hqrsjff
trim_galore --paired ${R1_FQ} ${R2_FQ}

module load hisat2
#hisat2-build ${GENOME} ${GENOME%.*}
hisat2 -p ${PROC} --no-spliced-alignment -x ${GENOME%.*} -1 ${R1_FQ%%.*}_val_1.fq.gz -2 ${R2_FQ%%.*}_val_2.fq.gz  -S ${R1_FQ%.*}.sam

module load samtools
samtools view --threads ${PROC} -b -o ${R1_FQ%.*}.bam ${R1_FQ%.*}.sam
samtools sort  -o ${R1_FQ%.*}_sorted.bam -T ${R1_FQ%.*}_temp --threads ${PROC} ${R1_FQ%.*}.bam

module load GIF2/java/1.8.0_25-b17
module load picard/2.17.0-ft5qztz
java -jar /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/picard-2.17.0-ft5qztzntoymuxiqt3b6yi6uqcmgzmds/bin/picard.jar CollectAlignmentSummaryMetrics  REFERENCE_SEQUENCE=${DBDIR}/${GENOME} INPUT=${R1_FQ%.*}_sorted.bam OUTPUT=${R1_FQ%.*}.bam_alignment.stats

########################################################################################################
```
