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
# sh runPilon.sh LongReads.fastq /work/GIF/remkv6/files genome.fa


PBReadsFq="$1"
DIR="$2"
GENOME="$3"

module load minimap2
minimap2 -L -ax asm10 ${GENOME} ${PBReadsFq}  >${GENOME%.*}.${PBReadsFq%.*}.sam

module load samtools
samtools view --threads 16 -b -o ${GENOME%.*}.${PBReadsFq%.*}.bam ${GENOME%.*}.${PBReadsFq%.*}.sam
samtools sort -m 7G -o ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam -T Round3PilonPB_temp --threads 16 ${GENOME%.*}.${PBReadsFq%.*}.bam
samtools index ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam
samtools flagstat ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam >flagstat.out
##############################################################################

8305392 + 0 in total (QC-passed reads + QC-failed reads)
4331864 + 0 secondary
1590664 + 0 supplementary
0 + 0 duplicates
8083125 + 0 mapped (97.32% : N/A)



```

### map ccsreads

```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/30_ReadMapping/01_CCS
ln -s ../../10_RepeatModeler/SCNgenome.fasta
ln -s /work/GIF/archive1/Baum/091615_SCN/PacBio_Jobs/016/016450/data/reads_of_insert.fastq

 echo sh  runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/30_ReadMapping/01_CCS SCNgenome.fasta >align.sh

#modified from above, not the same.
 #runPilon.sh
 ###############################################################################
 #!/bin/bash

 #You must provide the following. Note variable DBDIR does not need a "/" at the end.
 # sh runPilon.sh LongReads.fastq /work/GIF/remkv6/files genome.fa


 PBReadsFq="$1"
 DIR="$2"
 GENOME="$3"

 module load minimap2
 minimap2 -L -ax map-pb ${GENOME} ${PBReadsFq}  >${GENOME%.*}.${PBReadsFq%.*}.sam

 module load samtools
 samtools view --threads 16 -b -o ${GENOME%.*}.${PBReadsFq%.*}.bam ${GENOME%.*}.${PBReadsFq%.*}.sam
 samtools sort -m 7G -o ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam -T Round3PilonPB_temp --threads 16 ${GENOME%.*}.${PBReadsFq%.*}.bam
 samtools index ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam
 samtools flagstat ${GENOME%.*}.${PBReadsFq%.*}_sorted.bam >flagstat.out
 ##############################################################################

 290491 + 0 in total (QC-passed reads + QC-failed reads)
 124146 + 0 secondary
 76924 + 0 supplementary
 0 + 0 duplicates
 282382 + 0 mapped (97.21% : N/A)

```
### 260bp illumina reads
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/14_FinalizePilonShortReads

samtools flagstat MergedPE_sorted.bam >flagstat.out

371609941 + 0 in total (QC-passed reads + QC-failed reads)
88336631 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
273149405 + 0 mapped (73.50% : N/A)
283273310 + 0 paired in sequencing
141636655 + 0 read1
141636655 + 0 read2
145554472 + 0 properly paired (51.38% : N/A)
157444528 + 0 with itself and mate mapped
27368246 + 0 singletons (9.66% : N/A)
5516540 + 0 with mate mapped to a different chr
3523393 + 0 with mate mapped to a different chr (mapQ>=5)

## Remapping to polished genome.

for f in /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/14_FinalizePilonShortReads/Split*fq; do ln -s  $f;done
cp ~/common_scripts/runFeatureCountsNBigwig.sh .

paste <(ls -1 *R1*fq) <(ls -1 *R2*fq) |awk '{print "sh runFeatureCountsNBigwig.sh "$0" /work/GIF/remkv6/Baum/04_Dovetail2Restart/30_ReadMapping/03_MattsIllumina SCNgenome.fasta"}' |tr "\t" " " >align.sh

../SCNgenome.fasta
ml hisat2
hisat2-build SCNgenome.fasta SCNgenome

################################################################################
#runFeatureCounts.sh
#note: I premade the database so it would not be made again and again
###############################################################################
#!/bin/bash

#Working model.  Things to consider: the script needs different settings for 2 analyses. dna-seq (hisat2), as this is set up for RNA.  free node settings vs regular node settings (proc and memory).
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
hisat2 -p ${PROC} -x ${GENOME%.*} -1 ${R1_FQ%%.*}_val_1.fq -2 ${R2_FQ%%.*}_val_2.fq -S ${R1_FQ%.*}.sam

module load samtools
samtools view --threads ${PROC} -b -o ${R1_FQ%.*}.bam ${R1_FQ%.*}.sam
samtools sort  -o ${R1_FQ%.*}_sorted.bam -T ${R1_FQ%.*}_temp --threads ${PROC} ${R1_FQ%.*}.bam
samtools index  ${R1_FQ%.*}_sorted.bam
################################################################################

#number of reads mapped
less flagstat.out |grep "mapped (" |awk 'substr($1,1,1)!=0 {print $1}' |summary.sh
Total:  255,147,730
Count:  25
Mean:   10,205,909
Median: 10,249,030
Min:    9,703,494
Max:    10,474,570

#number of total reads passing QC
less flagstat.out |grep "total" |awk 'substr($1,1,1)!=0 {print $1}' |summary.sh
Total:  345,079,489
Count:  25
Mean:   13,803,179
Median: 13,837,573
Min:    13,437,925
Max:    13,932,134


#math
255,147,730/345,079,489=.7394 == 73.94% of reads mapped.  

#why such a low mapping rate?  run kraken against possible sources of contamination
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/30_ReadMapping/03_MattsIllumina/01_Kraken
ln -s ../../../../../Purcell/01_WhiteAbalone/05_Kraken/BactFungVirHum/
ln -s ../../../07_TN10IlluminaReads/01_TrimSequences/27_27799_CAGATC_L002_R1_001_val_1.fq
ln -s ../../../07_TN10IlluminaReads/01_TrimSequences/27_27799_CAGATC_L002_R2_001_val_2.fq

module load GIF/kraken2;kraken2 -db BactFungVirHum --threads 16 --report 27_27799_CAGATC_L002_R1_001_val_1.report --unclassified-out 27_27799_CAGATC_L002_R1_001_val_1#unclassified.fq --classified-out 27_27799_CAGATC_L002_R1_001_val_1#classified.fq --paired 27_27799_CAGATC_L002_R1_001_val_1.fq 27_27799_CAGATC_L002_R2_001_val_2.fq  > 27_27799_CAGATC_L002_R1_001_val_1.Kraken.out

19,647,012//345,079,489=.05693==5.69%
Hmm. so still low at ~80%

#try kraken to nematodes
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/30_ReadMapping/03_MattsIllumina/01_Kraken/01_nematodesDB
ln -s ../27_27799_CAGATC_L002_R1_001_val_1.fq
ln -s ../27_27799_CAGATC_L002_R2_001_val_2.fq
ln -s ../kraken_0.sub
ln -s ../../../../../../USDA/21_kraken/NematodeViral/

module load GIF/kraken2;kraken2 -db NematodeViral --threads 16 --report 27_27799_CAGATC_L002_R1_001_val_1.report --unclassified-out 27_27799_CAGATC_L002_R1_001_val_1#unclassified.fq --classified-out 27_27799_CAGATC_L002_R1_001_val_1#classified.fq --paired 27_27799_CAGATC_L002_R1_001_val_1.fq 27_27799_CAGATC_L002_R2_001_val_2.fq  > 27_27799_CAGATC_L002_R1_001_val_1.Kraken.out

```

## How much of the transcript libraries were aligned?

### EST's
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA
samtools flagstat SCNgenome.H.glycinesEST.bam

45915 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
44337 + 0 mapped (96.56% : N/A)

```
### Isoseq reads
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA

samtools flagstat SCNgenome.consensus_isoforms.bam

28490 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
26215 + 0 mapped (92.01% : N/A)
```

### RNAseq reads
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA

samtools flagstat -@ 16  AllRNASEQ_sorted.bam

1938458257 + 0 in total (QC-passed reads + QC-failed reads)
988104297 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
1859206493 + 0 mapped (95.91% : N/A)
950353960 + 0 paired in sequencing
475176980 + 0 read1
475176980 + 0 read2
830953842 + 0 properly paired (87.44% : N/A)
842316868 + 0 with itself and mate mapped
28785328 + 0 singletons (3.03% : N/A)
2419198 + 0 with mate mapped to a different chr
1555323 + 0 with mate mapped to a different chr (mapQ>=5)


```
