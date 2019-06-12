# Close gaps in the Juicebox assembly


#First round # technically the second.
```
module load GIF/pilon
java -Xmx120g -Djava.io.tmpdir=/scratch/remkv6 -jar /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/pilon-1.22-s7zrot6o5yqjh6oxpdxsxcdiswpjioyy/bin/pilon-1.22.jar --genome Genome.fasta --unpaired AllCCS.bam --output PilonGenome.fa --outdir /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/12_Pilon2ndRun/01_SecondRound --changes --vcf --tracks --fix all --threads 16 --mingap 0 --chunksize 100000

```

# Second Round
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/11_AlignMattsTN10Reads

#split the fastq reads for alignment
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/11_AlignMattsTN10Reads/01_hisat2
split -l 21000000 27_27799_CAGATC_L002_R1_001_val_1.fq SplitR1 &
split -l 21000000 27_27799_CAGATC_L002_R2_001_val_2.fq SplitR2 &

while read a b; do echo "sh runFeatureCountsNBigwig.sh "$a" "$b" /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/11_AlignMattsTN10Reads/02_Round2Align Genome.fasta" > align.sh

################################################################################
PROC=16
R1_FQ="$1"
R2_FQ="$2"
DBDIR="$3"
GENOME="$4"
GFF="$5"

#module load trimgalore
#module load py-setuptools/35.0.2-py2-hqrsjff
#trim_galore --paired ${R1_FQ} ${R2_FQ}

module load hisat2
#hisat2-build ${GENOME} ${GENOME%.*}
hisat2 -p ${PROC} -x ${GENOME%.*} -1 $R1_FQ -2 $R2_FQ -S ${R1_FQ%.*}.sam

module load samtools
samtools view --threads ${PROC} -b -o ${R1_FQ%.*}.bam ${R1_FQ%.*}.sam
samtools sort -m 7G -o ${R1_FQ%.*}_sorted.bam -T ${R1_FQ%.*}_temp --threads ${PROC} ${R1_FQ%.*}.bam

################################################################################

ls -1 *sorted.bam |tr "\n" " " |awk '{print "samtools merge -@ 16 MergedPE_sorted.bam "$0}' >samtoolsMerge.sh

################################################################################
samtools merge -@ 16 MergedPE_sorted.bam SplitR1aa_sorted.bam SplitR1ab_sorted.bam SplitR1ac_sorted.bam SplitR1ad_sorted.bam SplitR1ae_sorted.bam SplitR1af_sorted.bam SplitR1ag_sorted.bam SplitR1ah_sorted.bam SplitR1ai_sorted.bam SplitR1aj_sorted.bam SplitR1ak_sorted.bam SplitR1al_sorted.bam SplitR1am_sorted.bam SplitR1an_sorted.bam SplitR1ao_sorted.bam SplitR1ap_sorted.bam SplitR1aq_sorted.bam SplitR1ar_sorted.bam SplitR1as_sorted.bam SplitR1at_sorted.bam SplitR1au_sorted.bam SplitR1av_sorted.bam SplitR1aw_sorted.bam SplitR1ax_sorted.bam SplitR1ay_sorted.bam SplitR1az_sorted.bam SplitR1ba_sorted.bam
################################################################################


java -Xmx120g -Djava.io.tmpdir=/scratch/remkv6 -jar /work/GIF/software/programs/pilon/1.23/pilon-1.23.jar --genome Round3Pilon.fasta --frags MergedPE_sortedSorted.bam --unpaired Round3PilonPB_sorted.bam  --output PilonGenome.fa --outdir /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/12_Pilon2ndRun/03_ThirdRoundPilon --changes --vcf --tracks --fix all --threads 16 --mingap 0 --chunksize 100000
```

# Rounds 5-38 of pilon
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/05_Round4Pilon

sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/05_Round4Pilon Pilon4Round4.fa
################################################################################
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome11.fasta ;mv Genome11.Pilon.fasta Genome12.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome12.fasta ;mv Genome12.Pilon.fasta Genome13.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome13.fasta ;mv Genome13.Pilon.fasta Genome14.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome14.fasta ;mv Genome14.Pilon.fasta Genome15.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome15.fasta ;mv Genome15.Pilon.fasta Genome16.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome16.fasta ;mv Genome16.Pilon.fasta Genome17.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome17.fasta ;mv Genome17.Pilon.fasta Genome18.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome18.fasta ;mv Genome18.Pilon.fasta Genome19.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome19.fasta ;mv Genome19.Pilon.fasta Genome20.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome20.fasta ;mv Genome20.Pilon.fasta Genome21.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome21.fasta ;mv Genome21.Pilon.fasta Genome22.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome22.fasta ;mv Genome22.Pilon.fasta Genome23.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome23.fasta ;mv Genome23.Pilon.fasta Genome24.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome24.fasta ;mv Genome24.Pilon.fasta Genome25.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome25.fasta ;mv Genome25.Pilon.fasta Genome26.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome26.fasta ;mv Genome26.Pilon.fasta Genome27.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome27.fasta ;mv Genome27.Pilon.fasta Genome28.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome28.fasta ;mv Genome28.Pilon.fasta Genome29.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome29.fasta ;mv Genome29.Pilon.fasta Genome30.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome30.fasta ;mv Genome30.Pilon.fasta Genome31.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome31.fasta ;mv Genome31.Pilon.fasta Genome32.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome32.fasta ;mv Genome32.Pilon.fasta Genome33.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome33.fasta ;mv Genome33.Pilon.fasta Genome34.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome34.fasta ;mv Genome34.Pilon.fasta Genome35.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome35.fasta ;mv Genome35.Pilon.fasta Genome36.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome36.fasta ;mv Genome36.Pilon.fasta Genome37.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome37.fasta ;mv Genome37.Pilon.fasta Genome38.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome38.fasta ;mv Genome38.Pilon.fasta Genome39.fasta
sh runPilon.sh reads_of_insert.fastq /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/13_AlignCCSReads/13_Round11Pilon Genome39.fasta ;mv Genome39.Pilon.fasta Genome40.fasta

################################################################################

################################################################################
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
##############################################################################
```

#Round 39
```

# align split fastq reads to genome with hisat2, change to bam and sorted bam with samtools.

samtools merge -@ 16 MergedPE_sorted.bam SplitR1aa_sorted.bam SplitR1ab_sorted.bam SplitR1ac_sorted.bam SplitR1ad_sorted.bam SplitR1ae_sorted.bam SplitR1af_sorted.bam SplitR1ag_sorted.bam SplitR1ah_sorted.bam SplitR1ai_sorted.bam SplitR1aj_sorted.bam SplitR1ak_sorted.bam SplitR1al_sorted.bam SplitR1am_sorted.bam SplitR1an_sorted.bam SplitR1ao_sorted.bam SplitR1ap_sorted.bam SplitR1aq_sorted.bam SplitR1ar_sorted.bam SplitR1as_sorted.bam SplitR1at_sorted.bam SplitR1au_sorted.bam SplitR1av_sorted.bam SplitR1aw_sorted.bam SplitR1ax_sorted.bam SplitR1ay_sorted.bam SplitR1az_sorted.bam SplitR1ba_sorted.bam

module load GIF/pilon
java -Xmx110g -Djava.io.tmpdir=/scratch/remkv6 -jar /work/GIF/software/programs/pilon/1.23/pilon-1.23.jar --genome Genome40.fasta --frags MergedPE_sorted.bam --output Genome40.Pilon --outdir /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/14_FinalizePilonShortReads --changes --fix all --threads 16 --mingap 0 --chunksize 20000

```


#Results

```

---------------- Information for assembly 'Genome40.Pilon.fasta' ----------------


                                         Number of scaffolds          9
                                     Total size of scaffolds  157982452
                                            Longest scaffold   23986085
                                           Shortest scaffold   10884008
                                 Number of scaffolds > 1K nt          9 100.0%
                                Number of scaffolds > 10K nt          9 100.0%
                               Number of scaffolds > 100K nt          9 100.0%
                                 Number of scaffolds > 1M nt          9 100.0%
                                Number of scaffolds > 10M nt          9 100.0%
                                          Mean scaffold size   17553606
                                        Median scaffold size   16794675
                                         N50 scaffold length   17908190
                                          L50 scaffold count          4
                                                 scaffold %A      31.34
                                                 scaffold %C      18.12
                                                 scaffold %G      18.15
                                                 scaffold %T      31.33
                                                 scaffold %N       1.06
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs     100.0%
              Percentage of assembly in unscaffolded contigs       0.0%
                      Average number of contigs per scaffold      234.3
Average length of break (>25 Ns) between contigs in scaffold        797

                                           Number of contigs       2109
                              Number of contigs in scaffolds       2109
                          Number of contigs not in scaffolds          0
                                       Total size of contigs  156301043
                                              Longest contig     834407
                                             Shortest contig        221
                                   Number of contigs > 1K nt       2094  99.3%
                                  Number of contigs > 10K nt       1916  90.8%
                                 Number of contigs > 100K nt        483  22.9%
                                   Number of contigs > 1M nt          0   0.0%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size      74111
                                          Median contig size      41011
                                           N50 contig length     140610
                                            L50 contig count        320
                                                   contig %A      31.68
                                                   contig %C      18.31
                                                   contig %G      18.34
                                                   contig %T      31.66
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0

# The inital assembly prior to pilon

---------------- Information for assembly 'MisAssFixed.Pilon.FINAL.fasta' ----------------


                                         Number of scaffolds          9
                                     Total size of scaffolds  158541936
                                            Longest scaffold   24162073
                                           Shortest scaffold   10911090
                                 Number of scaffolds > 1K nt          9 100.0%
                                Number of scaffolds > 10K nt          9 100.0%
                               Number of scaffolds > 100K nt          9 100.0%
                                 Number of scaffolds > 1M nt          9 100.0%
                                Number of scaffolds > 10M nt          9 100.0%
                                          Mean scaffold size   17615771
                                        Median scaffold size   16842261
                                         N50 scaffold length   17946823
                                          L50 scaffold count          4
                                                 scaffold %A      31.32
                                                 scaffold %C      18.10
                                                 scaffold %G      18.13
                                                 scaffold %T      31.31
                                                 scaffold %N       1.13
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs     100.0%
              Percentage of assembly in unscaffolded contigs       0.0%
                      Average number of contigs per scaffold      264.0
Average length of break (>25 Ns) between contigs in scaffold        756

                                           Number of contigs       2376
                              Number of contigs in scaffolds       2376
                          Number of contigs not in scaffolds          0
                                       Total size of contigs  156744176
                                              Longest contig     811002
                                             Shortest contig        221
                                   Number of contigs > 1K nt       2344  98.7%
                                  Number of contigs > 10K nt       2140  90.1%
                                 Number of contigs > 100K nt        462  19.4%
                                   Number of contigs > 1M nt          0   0.0%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size      65970
                                          Median contig size      38009
                                           N50 contig length     116317
                                            L50 contig count        379
                                                   contig %A      31.68
                                                   contig %C      18.31
                                                   contig %G      18.34
                                                   contig %T      31.67
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0

```
