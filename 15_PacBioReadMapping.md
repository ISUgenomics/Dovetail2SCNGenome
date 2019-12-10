# Assess genome assembly quality by read mapping statistics


### map subreads
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/30_ReadMapping/02_Subread
ln -s ../../10_RepeatModeler/SCNgenome.fasta
cat /work/GIF/archive1/Purcell/SeriolaRivoliana/pacbio/*/*/Analysis_Results/*fastq >AllSubreads.fastq



sh runPilon.sh AllSubreads.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/30_ReadMapping/02_Subread SCNgenome.fasta


#runPilon.sh
###############################################################################
#!/bin/bash

#You must provide the following. Note variable DBDIR does not need a "/" at the end.
# sh runPilon.sh LongReads.fastq /work/GIF/remkv6/files genome.fa ShortReadsR1.fq ShortReadsR2.fq


PBReadsFq="$1"
DIR="$2"
GENOME="$3"
R1_FQ="$4"
R2_FQ="$5"

#module load hisat2
#hisat2-build ${GENOME} ${GENOME%.*}
#hisat2 -p 16 -x ${GENOME%.*} -1 $R1_FQ -2 $R2_FQ -S ${R1_FQ%.*}.sam
#module load samtools
#samtools view --threads 16 -b -o ${GENOME%.*}.${R1_FQ%.*}.bam ${GENOME%.*}.${R1_FQ%.*}.sam
#samtools sort -m 7G -o ${GENOME%.*}.${R1_FQ%.*}_sorted.bam -T Round3PilonPB_temp --threads 16 ${GENOME%.*}.${R1_FQ%.*}.bam
#samtools index ${GENOME%.*}.${R1_FQ%.*}_sorted.bam


module load minimap2
minimap2 -L -ax asm10 ${GENOME} ${PBReadsFq}  >${GENOME%.*}.${PBReadsFq%.*}.sam

module load samtools
samtools view --threads 16 -b -o ${GENOME%.*}.${PBReadsFq%.*}.bam ${GENOME%.*}.${PBReadsFq%.*}.sam
samtools sort -m 7G -o ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam -T Round3PilonPB_temp --threads 16 ${GENOME%.*}.${PBReadsFq%.*}.bam
samtools index ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam

module load pilon
java -Xmx120g -Djava.io.tmpdir=/scratch/remkv6 -jar /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/pilon-1.22-s7zrot6o5yqjh6oxpdxsxcdiswpjioyy/bin/pilon-1.22.jar --genome ${GENOME} --unpaired ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam --output ${GENOME%.*}.Pilon --outdir ${DIR} --changes --fix local,gaps --threads 16 --mingap 0 --chunksize 300000
################################################################################################################################################



```

### map ccsreads

```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/30_ReadMapping/01_CCS
ln -s ../../10_RepeatModeler/SCNgenome.fasta
ln -s /work/GIF/archive1/Baum/091615_SCN/PacBio_Jobs/016/016450/data/reads_of_insert.fastq

 echo sh  runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/30_ReadMapping/01_CCS SCNgenome.fasta >align.sh

 #runPilon.sh is the same as above
```
