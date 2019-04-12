#  Received second SCN genome from dovetail and since it appears to be more complete, lets see if we can scaffold it further and fill the gaps with pacbio reads.

```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling
ln -s /work/GIF/GIF3/archive3/Baum/091615_SCN/PacBio_Jobs/016/016450/data/reads_of_insert.fastq
ln -s /work/GIF/remkv6/Baum/04_Dovetail2Restart/01_GenomeDownload/nematode_sp_19Jul2018_IbtP1.fasta
cpan
install Perl4::CoreLibs
source ~/.bashrc
wget  https://www.baseclear.com/wp-content/uploads/BaseTools-License-agreement-v.-2014-03-31.pdf
wget https://www.baseclear.com/wp-content/uploads/SSPACE-longread-v.-1-1.tar.gz
tar -zxvf SSPACE-longread-v.-1-1.tar.gz
```

## scaffold assembly with ccs reads
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/01_ScaffoldCCS
echo "module load sspace-longread/1.1-5ae3ijg" >>gapfiller.sh
echo "perl SSPACE-LongRead.pl -c ../nematode_sp_19Jul2018_IbtP1.fasta -p ../reads_of_insert.fastq -t 16 -s ../ccs2genome.out" >>gapfiller.sh

RunTime=04:56:13 -- 1 node 16 processors

#assembly prior to scaffolding
#############################################################################

                                         Number of scaffolds        590
                                     Total size of scaffolds  166433103
                                            Longest scaffold   20440036
                                           Shortest scaffold        536
                                 Number of scaffolds > 1K nt        589  99.8%
                                Number of scaffolds > 10K nt        580  98.3%
                               Number of scaffolds > 100K nt        135  22.9%
                                 Number of scaffolds > 1M nt         16   2.7%
                                Number of scaffolds > 10M nt          6   1.0%
                                          Mean scaffold size     282090
                                        Median scaffold size      44252
                                         N50 scaffold length   10944860
                                          L50 scaffold count          6
                                                 scaffold %A      31.58
                                                 scaffold %C      18.39
                                                 scaffold %G      18.39
                                                 scaffold %T      31.56
                                                 scaffold %N       0.08
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      90.8%
              Percentage of assembly in unscaffolded contigs       9.2%
                      Average number of contigs per scaffold        3.3
Average length of break (>25 Ns) between contigs in scaffold        100

                                           Number of contigs       1951
                              Number of contigs in scaffolds       1577
                          Number of contigs not in scaffolds        374
                                       Total size of contigs  166297003
                                              Longest contig    1155511
                                             Shortest contig        173
                                   Number of contigs > 1K nt       1937  99.3%
                                  Number of contigs > 10K nt       1896  97.2%
                                 Number of contigs > 100K nt        449  23.0%
                                   Number of contigs > 1M nt          3   0.2%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size      85237
                                          Median contig size      50218
                                           N50 contig length     137872
                                            L50 contig count        302
                                                   contig %A      31.60
                                                   contig %C      18.40
                                                   contig %G      18.41
                                                   contig %T      31.58
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0
#############################################################################

#assembly post scaffolding
#############################################################################
Number of scaffolds        527
Total size of scaffolds  166493275
   Longest scaffold   31384897
  Shortest scaffold        536
Number of scaffolds > 1K nt        526  99.8%
Number of scaffolds > 10K nt        518  98.3%
Number of scaffolds > 100K nt        135  25.6%
Number of scaffolds > 1M nt         15   2.8%
Number of scaffolds > 10M nt          5   0.9%
 Mean scaffold size     315927
Median scaffold size      49377
N50 scaffold length   11380747
 L50 scaffold count          5
        scaffold %A      31.57
        scaffold %C      18.38
        scaffold %G      18.39
        scaffold %T      31.54
        scaffold %N       0.12
scaffold %non-ACGTN       0.00
Number of scaffold non-ACGTN nt          0

Percentage of assembly in scaffolded contigs      92.4%
Percentage of assembly in unscaffolded contigs       7.6%
Average number of contigs per scaffold        3.6
Average length of break (>25 Ns) between contigs in scaffold        141

  Number of contigs       1917
Number of contigs in scaffolds       1601
Number of contigs not in scaffolds        316
Total size of contigs  166297037
     Longest contig    1155511
    Shortest contig        173
Number of contigs > 1K nt       1904  99.3%
Number of contigs > 10K nt       1863  97.2%
Number of contigs > 100K nt        459  23.9%
Number of contigs > 1M nt          3   0.2%
Number of contigs > 10M nt          0   0.0%
   Mean contig size      86749
 Median contig size      51246
  N50 contig length     143023
   L50 contig count        298
          contig %A      31.61
          contig %C      18.40
          contig %G      18.41
          contig %T      31.58
          contig %N       0.00
  contig %non-ACGTN       0.00
Number of contig non-ACGTN nt          0
#############################################################################

This knocked out 63 contigs
```

## Lets try scaffolding again with the subreads

```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/02_ScaffoldSubreads

#subreads
ln -s /work/GIF/severin/Baum/09_Falcon_by_transcript/SCN.all.subreads.sl.7k.fasta
ln -s ../nematode_sp_19Jul2018_IbtP1.fasta

module load sspace-longread/1.1-5ae3ijg
perl SSPACE-LongRead.pl -c nematode_sp_19Jul2018_IbtP1.fasta -p SCN.all.subreads.sl.7k.fasta -t 16


#############################################################################
Number of scaffolds        396
Total size of scaffolds  166863111
   Longest scaffold   31384913
  Shortest scaffold        536
Number of scaffolds > 1K nt        395  99.7%
Number of scaffolds > 10K nt        390  98.5%
Number of scaffolds > 100K nt        138  34.8%
Number of scaffolds > 1M nt         15   3.8%
Number of scaffolds > 10M nt          5   1.3%
 Mean scaffold size     421371
Median scaffold size      64938
N50 scaffold length   11380747
 L50 scaffold count          5
        scaffold %A      31.50
        scaffold %C      18.35
        scaffold %G      18.34
        scaffold %T      31.47
        scaffold %N       0.34
scaffold %non-ACGTN       0.00
Number of scaffold non-ACGTN nt          0

Percentage of assembly in scaffolded contigs      95.6%
Percentage of assembly in unscaffolded contigs       4.4%
Average number of contigs per scaffold        4.8
Average length of break (>25 Ns) between contigs in scaffold        375

  Number of contigs       1902
Number of contigs in scaffolds       1711
Number of contigs not in scaffolds        191
Total size of contigs  166297033
     Longest contig    1155511
    Shortest contig        173
Number of contigs > 1K nt       1890  99.4%
Number of contigs > 10K nt       1853  97.4%
Number of contigs > 100K nt        459  24.1%
Number of contigs > 1M nt          3   0.2%
Number of contigs > 10M nt          0   0.0%
   Mean contig size      87433
 Median contig size      52046
  N50 contig length     143481
   L50 contig count        298
          contig %A      31.61
          contig %C      18.41
          contig %G      18.40
          contig %T      31.58
          contig %N       0.00
  contig %non-ACGTN       0.00
Number of contig non-ACGTN nt          0
#############################################################################

```

### Compare top 15 largest scaffolds
```
paste <(bioawk -c fastx '{print $name,length($seq)}' ../nematode_sp_19Jul2018_IbtP1.fasta |sort -k2,2nr) <(bioawk -c fastx '{print $name,length($seq)}' ../../01_ScaffoldCCS/PacBio_scaffolder_results/scaffolds.fasta |sort -k2,2nr)  <(bioawk -c fastx '{print $name,length($seq)}' scaffolds.fasta |sort -k2,2nr ) <(bioawk -c fastx '{print $name,length($seq)}' /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/01_DovetailOutput/download/serial_chicago_hic/nematode_sp._22Aug2017_DZkUC.fasta |sort -k2,2nr ) |awk 'NR<16' |cut -f 2,4,6,8 |cat <(printf "dovetail2unscaffolded\tdovetail2scaffoldccs\tdovetail2scaffoldsubread\tdovetail1738genome\n") - |less -S
###############################################################################
dovetail2unscaffolded   dovetail2scaffoldccs    dovetail2scaffoldsubread        dovetail1738genome
20440036        31384897        31384913        17452550
16328168        16328168        16355503        15749149
14072926        14072926        14177380        14051119
13519465        13743806        13779824        12619283
11380747        11380747        11380747        11216375
10944860        9630264 9630264 10460252
9630264 7744779 7744779 10159891
7744779 6172682 6172682 9664129
6172682 3481826 3576988 7104694
3335668 3105541 3154034 2044697
3105541 3044928 3105541 2003123
3044928 2723444 2723444 1961110
2723444 2528721 2551054 1667135
2448955 1823045 2010315 439826
1823045 1097846 1097846 397866
###############################################################################

Looks like scaffolding with subreads is the best.  
```

## Fill the gaps in the scaffolded assembly
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/03_GapFilling
ln -s ../reads_of_insert.fastq
ln -s ../02_ScaffoldSubreads/PacBio_scaffolder_results/scaffolds.fasta

# find the breakpoinks, snps, and indels -- this took 10 mins on interactive node
/opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/mindthegap-2.0.2-7jl4b5gf24sy56r7b2bsmtkczifs7hyd/bin/MindTheGap find -ref scaffolds.fasta -in reads_of_insert.fastq

# fill the gaps
module load mindthegap/2.0.2-7jl4b5g -- this may take 10-12 hrs according to the estimate from the program.
/opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/mindthegap-2.0.2-7jl4b5gf24sy56r7b2bsmtkczifs7hyd/bin/MindTheGap fill -bkpt MindTheGap_Expe-2018-10-11.10\:30.breakpoints -in reads_of_insert.fastq -out MTGFilledScaffolds -max-memory 120000

This just provided some structural variants, instead of filling gaps.  This is outputted in vcf format and may be useful.
```


##  Really fill the gaps in the scaffolded assembly

```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/04_RedundansGapFilling

#this should fill gaps, unless gaps are not closed with long reads.
module load redundans
redundans.py -o redundansOut --noreduction -l reads_of_insert.fastq -f scaffolds.fasta -t 16


#This collapsed too many scaffolds, try gmcloser
```

## GMcloser to fill gaps with pacbio reads
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/05_gmcloserGapfilling
ln -s ../reads_of_insert.fasta
ln -s ../02_ScaffoldSubreads/PacBio_scaffolder_results/scaffolds.fasta
ln -s ../../07_TN10IlluminaReads/27_27799_CAGATC_L002_R1_001.fastq
ln -s ../../07_TN10IlluminaReads/27_27799_CAGATC_L002_R2_001.fastq
#!/bin/bash
#SBATCH -A its-hpc-condo-las-free
#SBATCH -N 1
#SBATCH -p freecompute
#SBATCH --ntasks-per-node=12
#SBATCH -t 168:00:00
#SBATCH -J gmcloser_0
#SBATCH -o gmcloser_0.o%j
#SBATCH -e gmcloser_0.e%j
#SBATCH --mail-user=remkv6@istate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
module use /work/GIF/software/modules
module load GIF/gmcloser
module load GIF/perl
gmcloser  -t scaffolds.fasta -q reads_of_insert.fasta -ht -c -n 12 -lr -it 3 -mq 2 -p out4 -l 100 -i 500 -d 50 -c
scontrol show job $SLURM_JOB_ID



Results
################################################################################
#beginning assembly
new_Assemblathon.pl scaffolds.fasta
Number of scaffolds        396
Total size of scaffolds  166863111
    Longest scaffold   31384913
   Shortest scaffold        536
Number of scaffolds > 1K nt        395  99.7%
Number of scaffolds > 10K nt        390  98.5%
Number of scaffolds > 100K nt        138  34.8%
Number of scaffolds > 1M nt         15   3.8%
Number of scaffolds > 10M nt          5   1.3%
  Mean scaffold size     421371
Median scaffold size      64938
 N50 scaffold length   11380747
  L50 scaffold count          5
         scaffold %A      31.50
         scaffold %C      18.35
         scaffold %G      18.34
         scaffold %T      31.47
         scaffold %N       0.34
 scaffold %non-ACGTN       0.00

#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/05_gmcloserGapfilling/iteration-3
new_Assemblathon.pl out4.gapclosed.fa
Number of scaffolds        396
                                    Total size of scaffolds  166141519
                                           Longest scaffold   31274080
                                          Shortest scaffold        536
                                Number of scaffolds > 1K nt        395  99.7%
                               Number of scaffolds > 10K nt        390  98.5%
                              Number of scaffolds > 100K nt        138  34.8%
                                Number of scaffolds > 1M nt         15   3.8%
                               Number of scaffolds > 10M nt          5   1.3%
                                         Mean scaffold size     419549
                                       Median scaffold size      64938
                                        N50 scaffold length   11339714
                                         L50 scaffold count          5
                                                scaffold %A      31.50
                                                scaffold %C      18.35
                                                scaffold %G      18.34
                                                scaffold %T      31.47
                                                scaffold %N       0.34
                                        scaffold %non-ACGTN       0.00
                            Number of scaffold non-ACGTN nt          0
################################################################################

GMcloser didnt do a very good job.  using Pilon for gap filling now
```

## Try to fill in more of those N's using Pilon
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/07_TN10IlluminaReads

#Need to align the TN10 reads Matt gave me to the sspace scaffolded genome

#Create hisat database
module load hisat2
hisat2-build scaffolds.fasta sspaceScaffolded


runHisat2.sh
#################################################################################
#!/bin/bash

module load hisat2
module load samtools
DBDIR="/work/GIF/remkv6/Baum/04_Dovetail2Restart/07_TN10IlluminaReads"
GENOME="sspaceScaffolded"

p=12
R1_FQ="$1"
R2_FQ="$2"

OUTPUT=$(basename ${R1_FQ} |cut -f 1 -d "_");

hisat2 \
  -p ${p} \
  -x ${DBDIR}/${GENOME} \
  -1 ${R1_FQ} \
  -2 ${R2_FQ}  \
  -S  ${OUTPUT}.sam &> ${OUTPUT}.log
samtools view --threads 12 -b -o ${OUTPUT}.bam ${OUTPUT}.sam
samtools sort -m 5G -o ${OUTPUT}_sorted.bam -T ${OUTPUT}_temp --threads 12 ${OUTPUT}.bam
######################################################################################

sh runHISAT2.sh 27_27799_CAGATC_L002_R1_001.fastq 27_27799_CAGATC_L002_R2_001.fastq
```


## Now run Pilon with the bam file
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/07_Split1Pilon

#this runs a bit faster if I split the genome across two nodes


#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/07_Split1Pilon
ln -s ../04_GapFilling/02_ScaffoldSubreads/PacBio_scaffolder_results/scaffolds.fasta Gmcloser.Scaffolds.fa
module load cdbfasta
cdbfasta Gmcloser.Scaffolds.fa
awk '{print $1}' Gmcloser.Scaffolds1.bed |cdbyank Gmcloser.Scaffolds.fa.cidx > Gmcloser.Scaffolds1.fa

module load samtools
samtools view -@ 16 -b -L Gmcloser.Scaffolds1.bed ../06_pilonPolishGapFill/27_sorted.bam >Gmcloser.Scaffolds1.bam
samtools index Gmcloser.Scaffolds1.bam
module use /work/GIF/software/modules
module load pilon/1.22-s7zrot6
/opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/jdk-8u172-b11-rnauqmrxgpic4aopyat6bomafogtjavz/bin/java -Xms120000M -Xmx120000M -Xss2000M -jar /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/pilon-1.22-s7zrot6o5yqjh6oxpdxsxcdiswpjioyy/bin/pilon-1.22.jar  --genome Gmcloser.Scaffolds1.fa  --frags Gmcloser.Scaffolds1.bam --output Pilon1 --outdir /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/07_Split1Pilon/01_misassemblyTest --changes --vcf --tracks --fix all --threads 16 --mingap 0 --chunksize 1000000

#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/08_Split2Pilon
ln -s ../04_GapFilling/02_ScaffoldSubreads/PacBio_scaffolder_results/scaffolds.fasta Gmcloser.Scaffolds.fa
module load cdbfasta
cdbfasta Gmcloser.Scaffolds.fa
awk '{print $1}' Gmcloser.Scaffolds2.bed |cdbyank Gmcloser.Scaffolds.fa.cidx > Gmcloser.Scaffolds2.fa

module load samtools
samtools view -@ 16 -b -L Gmcloser.Scaffolds1.bed ../06_pilonPolishGapFill/27_sorted.bam >Gmcloser.Scaffolds1.bam
samtools index Gmcloser.Scaffolds1.bam

module use /work/GIF/software/modules
module load pilon/1.22-s7zrot6
/opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/jdk-8u172-b11-rnauqmrxgpic4aopyat6bomafogtjavz/bin/java -Xms120000M -Xmx120000M -Xss2000M -jar /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/pilon-1.22-s7zrot6o5yqjh6oxpdxsxcdiswpjioyy/bin/pilon-1.22.jar  --genome Gmcloser.Scaffolds2.fa  --frags Gmcloser.Scaffolds2.bam --output Pilon2 --outdir /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/06_pilonPolishGapFill --changes --vcf --tracks --fix all --threads 16 --mingap 0 --chunksize 300000
```
