# get featurecounts for SCN reads


/work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate

### get files
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/39_DifferentialExpression/02_mRNA

ln -s /work/GIF/remkv6/Baum/04_Dovetail2Restart/49_RenameChromosomes/01_Transfer2Box/OrderedSCNGenePredictions.gff3
ln -s ../../49_RenameChromosomes/01_Transfer2Box/SCNgenome.fasta

#grab all the trimmed reads from the previous run on a differently named genome
for f  in ../../11_AlignRNA/*val_1.fq.gz; do ln -s $f;done
for f  in ../../11_AlignRNA/*val_2.fq.gz; do ln -s $f;done

#One of the samples had odd trimming and had to be left untrimmed because it would not map.
ln -s /work/GIF/remkv6/Baum/03_GlandRepeat738/00_TrimReads/1703FL-02-01_S1_L001_R2_001.fastq.gz
ln -s /work/GIF/remkv6/Baum/03_GlandRepeat738/00_TrimReads/1703FL-02-01_S1_L001_R1_001.fastq.gz
```
### figure out stranding Information
```
ml miniconda3; source activate RSeQC

#need to create a bed12 annotation file


#I already have bams aligned unstranded, so I used those.
 for f in ../../11_AlignRNA/*bam; do echo $f; infer_experiment.py -i $f -r mikado.loci.ancestral.bed12  >${f}.Strand  ;done

 ../../11_AlignRNA/1703FL-02-01_S1_L001_R1_001_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.0962
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.8428
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.0610
 ../../11_AlignRNA/1703FL-02-02_S2_L001_R1_001_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.1704
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.7922
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.0373
 ../../11_AlignRNA/1703FL-02-03_S3_L001_R1_001_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.1324
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.8398
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.0278
 ../../11_AlignRNA/1703FL-02-04_S4_L001_R1_001_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.1149
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.8497
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.0354
 ../../11_AlignRNA/1703-TM101_S0_L003_R1_001_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.0993
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.8351
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.0656
 ../../11_AlignRNA/1703-TM102_S0_L003_R1_001_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.1150
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.8491
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.0359
 ../../11_AlignRNA/AllRNASEQ_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.2571
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.5031
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.2398
 ../../11_AlignRNA/SCNgenome.consensus_isoforms_sorted.bam


 This is SingleEnd Data
 Fraction of reads failed to determine: 0.1303
 Fraction of reads explained by "++,--": 0.8118
 Fraction of reads explained by "+-,-+": 0.0580
 ../../11_AlignRNA/SCNgenome.H.glycinesEST_sorted.bam


 This is SingleEnd Data
 Fraction of reads failed to determine: 0.0928
 Fraction of reads explained by "++,--": 0.7269
 Fraction of reads explained by "+-,-+": 0.1803
 ../../11_AlignRNA/SRR6230579_1_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.1536
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.8245
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.0220
 ../../11_AlignRNA/SRR6230580_1_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.0796
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.0229
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.8974
 ../../11_AlignRNA/SRR6230581_1_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.1618
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.0187
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.8195
 ../../11_AlignRNA/SRR6230582_1_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.1507
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.8252
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.0240
 ../../11_AlignRNA/SRR6230583_1_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.1450
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.8334
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.0216
 ../../11_AlignRNA/SRR6230584_1_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.1424
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.8339
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.0238
 ../../11_AlignRNA/SRR6230585_1_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.1442
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.0205
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.8353
 ../../11_AlignRNA/SRR6230586_1_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.1259
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.8562
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.0179
 ../../11_AlignRNA/SRR6230587_1_val_1.fq_sorted.bam


 This is PairEnd Data
 Fraction of reads failed to determine: 0.1234
 Fraction of reads explained by "1++,1--,2+-,2-+": 0.8559
 Fraction of reads explained by "1+-,1-+,2++,2--": 0.0206

```

### align to featurecounts
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/39_DifferentialExpression/02_mRNA
paste <(ls -1 *l_1.fq.gz)  <(ls -1 *l_2.fq.gz) |while read line; do echo "sh runFeatureCountsNBigwigStrandRFmRNA.sh "$line" /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3";done >>featurecounts.sh

#Ran the above script for each of the required analyses, deleted improper script/read pairs for those that required different strandedness settings from the rest of the samples

runFeatureCountsNBigwigStrandFRgeneMult.sh  runFeatureCountsNBigwigStrandFRmRNA.sh      runFeatureCountsNBigwigStrandRFgene.sh
runFeatureCountsNBigwigStrandFRgene.sh      runFeatureCountsNBigwigStrandRFgeneMult.sh  runFeatureCountsNBigwigStrandRFmRNA.sh


#runFeatureCountsNBigwigStrandFRgene.sh
################################################################################################################################


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
hisat2 -p ${PROC} --rna-strandness FR -x ${GENOME%.*} -1 ${R1_FQ} -2 ${R2_FQ}  -S ${R1_FQ%.*}Gene.sam

module load samtools
samtools view --threads ${PROC} -b -o ${R1_FQ%.*}Gene.bam ${R1_FQ%.*}Gene.sam
samtools sort  -o ${R1_FQ%.*}Gene_sorted.bam -T ${R1_FQ%.*}Gene_temp --threads ${PROC} ${R1_FQ%.*}Gene.bam

module load GIF2/java/1.8.0_25-b17
module load picard/2.17.0-ft5qztz
java -jar /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/picard-2.17.0-ft5qztzntoymuxiqt3b6yi6uqcmgzmds/bin/picard.jar CollectAlignmentSummaryMetrics  REFERENCE_SEQUENCE=${DBDIR}/${GENOME} INPUT=${R1_FQ%.*}Gene_sorted.bam OUTPUT=${R1_FQ%.*}Gene.bam_alignment.stats

module load subread
featureCounts -T ${PROC} -M   -t gene -s 1 -g ID -a ${GFF} -o ${R1_FQ%.*}Gene_counts_genes.txt ${R1_FQ%.*}Gene_sorted.bam

module load bedtools2
bedtools genomecov -ibam ${R1_FQ%.*}Gene_sorted.bam -bga -g ${GENOME} >${R1_FQ%.*}Gene_sorted.bam.bdg
#module load bioawk
#bioawk -c fastx '{print $name,length($seq)}'  ${GENOME} >Chr.sizes
/home/remkv6/common_scripts/kentUtils/bin/linux.x86_64/bedGraphToBigWig ${R1_FQ%.*}Gene_sorted.bam.bdg Chr.sizes ${R1_FQ%.*}Gene_sorted.bw
#######################################################################################################################

#the resulting scripts that were Ran
#featurecounts.sh
#########################################################################################################################
sh runFeatureCountsNBigwigStrandFRgeneMult.sh 1703FL-02-01_S1_L001_R1_001_val_1.fq.gz 1703FL-02-01_S1_L001_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgeneMult.sh 1703FL-02-02_S2_L001_R1_001_val_1.fq.gz 1703FL-02-02_S2_L001_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgeneMult.sh 1703FL-02-03_S3_L001_R1_001_val_1.fq.gz 1703FL-02-03_S3_L001_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgeneMult.sh 1703FL-02-04_S4_L001_R1_001_val_1.fq.gz 1703FL-02-04_S4_L001_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgeneMult.sh 1703-TM101_S0_L003_R1_001_val_1.fq.gz 1703-TM101_S0_L003_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgeneMult.sh 1703-TM102_S0_L003_R1_001_val_1.fq.gz 1703-TM102_S0_L003_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgeneMult.sh SRR6230579_1_val_1.fq.gz SRR6230579_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgeneMult.sh SRR6230582_1_val_1.fq.gz SRR6230582_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgeneMult.sh SRR6230583_1_val_1.fq.gz SRR6230583_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgeneMult.sh SRR6230584_1_val_1.fq.gz SRR6230584_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgeneMult.sh SRR6230586_1_val_1.fq.gz SRR6230586_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgeneMult.sh SRR6230587_1_val_1.fq.gz SRR6230587_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandRFgeneMult.sh SRR6230580_1_val_1.fq.gz SRR6230580_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandRFgeneMult.sh SRR6230581_1_val_1.fq.gz SRR6230581_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandRFgeneMult.sh SRR6230585_1_val_1.fq.gz SRR6230585_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgene.sh 1703FL-02-01_S1_L001_R1_001_val_1.fq.gz 1703FL-02-01_S1_L001_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgene.sh 1703FL-02-02_S2_L001_R1_001_val_1.fq.gz 1703FL-02-02_S2_L001_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgene.sh 1703FL-02-03_S3_L001_R1_001_val_1.fq.gz 1703FL-02-03_S3_L001_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgene.sh 1703FL-02-04_S4_L001_R1_001_val_1.fq.gz 1703FL-02-04_S4_L001_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgene.sh 1703-TM101_S0_L003_R1_001_val_1.fq.gz 1703-TM101_S0_L003_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgene.sh 1703-TM102_S0_L003_R1_001_val_1.fq.gz 1703-TM102_S0_L003_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgene.sh SRR6230579_1_val_1.fq.gz SRR6230579_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgene.sh SRR6230582_1_val_1.fq.gz SRR6230582_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgene.sh SRR6230583_1_val_1.fq.gz SRR6230583_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgene.sh SRR6230584_1_val_1.fq.gz SRR6230584_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgene.sh SRR6230586_1_val_1.fq.gz SRR6230586_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRgene.sh SRR6230587_1_val_1.fq.gz SRR6230587_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandRFgene.sh SRR6230580_1_val_1.fq.gz SRR6230580_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandRFgene.sh SRR6230581_1_val_1.fq.gz SRR6230581_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandRFgene.sh SRR6230585_1_val_1.fq.gz SRR6230585_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRmRNA.sh 1703FL-02-01_S1_L001_R1_001_val_1.fq.gz 1703FL-02-01_S1_L001_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRmRNA.sh 1703FL-02-02_S2_L001_R1_001_val_1.fq.gz 1703FL-02-02_S2_L001_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRmRNA.sh 1703FL-02-03_S3_L001_R1_001_val_1.fq.gz 1703FL-02-03_S3_L001_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRmRNA.sh 1703FL-02-04_S4_L001_R1_001_val_1.fq.gz 1703FL-02-04_S4_L001_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRmRNA.sh 1703-TM101_S0_L003_R1_001_val_1.fq.gz 1703-TM101_S0_L003_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRmRNA.sh 1703-TM102_S0_L003_R1_001_val_1.fq.gz 1703-TM102_S0_L003_R2_001_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRmRNA.sh SRR6230579_1_val_1.fq.gz SRR6230579_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRmRNA.sh SRR6230582_1_val_1.fq.gz SRR6230582_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRmRNA.sh SRR6230583_1_val_1.fq.gz SRR6230583_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRmRNA.sh SRR6230584_1_val_1.fq.gz SRR6230584_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRmRNA.sh SRR6230586_1_val_1.fq.gz SRR6230586_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandFRmRNA.sh SRR6230587_1_val_1.fq.gz SRR6230587_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandRFmRNA.sh SRR6230580_1_val_1.fq.gz SRR6230580_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandRFmRNA.sh SRR6230581_1_val_1.fq.gz SRR6230581_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
sh runFeatureCountsNBigwigStrandRFmRNA.sh SRR6230585_1_val_1.fq.gz SRR6230585_2_val_2.fq.gz /work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate SCNgenome.fasta OrderedSCNGenePredictions.gff3
#########################################################################################################################
```


### Mapping rates
```
paste <(ls -1 *stats) <(cat *stats |awk '$1=="PAIR" {print $7}' ) |less

1703FL-02-01_S1_L001_R1_001.fastqGene.bam_alignment.stats       0.756435
1703FL-02-01_S1_L001_R1_001.fastqmRNA.bam_alignment.stats       0.756435
1703FL-02-01_S1_L001_R1_001.fastqMult.bam_alignment.stats       0.756435
1703FL-02-02_S2_L001_R1_001_val_1.fqGene.bam_alignment.stats    0.893617
1703FL-02-02_S2_L001_R1_001_val_1.fqmRNA.bam_alignment.stats    0.893617
1703FL-02-02_S2_L001_R1_001_val_1.fqMult.bam_alignment.stats    0.893617
1703FL-02-03_S3_L001_R1_001_val_1.fqGene.bam_alignment.stats    0.9458
1703FL-02-03_S3_L001_R1_001_val_1.fqmRNA.bam_alignment.stats    0.9458
1703FL-02-03_S3_L001_R1_001_val_1.fqMult.bam_alignment.stats    0.9458
1703FL-02-04_S4_L001_R1_001_val_1.fqGene.bam_alignment.stats    0.95839
1703FL-02-04_S4_L001_R1_001_val_1.fqmRNA.bam_alignment.stats    0.95839
1703FL-02-04_S4_L001_R1_001_val_1.fqMult.bam_alignment.stats    0.95839
1703-TM101_S0_L003_R1_001_val_1.fqGene.bam_alignment.stats      0.952391
1703-TM101_S0_L003_R1_001_val_1.fqmRNA.bam_alignment.stats      0.952391
1703-TM101_S0_L003_R1_001_val_1.fqMult.bam_alignment.stats      0.952391
1703-TM102_S0_L003_R1_001_val_1.fqGene.bam_alignment.stats      0.924849
1703-TM102_S0_L003_R1_001_val_1.fqmRNA.bam_alignment.stats      0.924849
1703-TM102_S0_L003_R1_001_val_1.fqMult.bam_alignment.stats      0.924849
SRR6230579_1_val_1.fqGene.bam_alignment.stats   0.933199
SRR6230579_1_val_1.fqmRNA.bam_alignment.stats   0.933199
SRR6230579_1_val_1.fqMult.bam_alignment.stats   0.933199
SRR6230580_1_val_1.fqGene.bam_alignment.stats   0.936017
SRR6230580_1_val_1.fqmRNA.bam_alignment.stats   0.936017
SRR6230580_1_val_1.fqMult.bam_alignment.stats   0.936017
SRR6230581_1_val_1.fqGene.bam_alignment.stats   0.949976
SRR6230581_1_val_1.fqmRNA.bam_alignment.stats   0.949976
SRR6230581_1_val_1.fqMult.bam_alignment.stats   0.949976
SRR6230582_1_val_1.fqGene.bam_alignment.stats   0.951359
SRR6230582_1_val_1.fqmRNA.bam_alignment.stats   0.951359
SRR6230582_1_val_1.fqMult.bam_alignment.stats   0.951359
SRR6230583_1_val_1.fqGene.bam_alignment.stats   0.955742
SRR6230583_1_val_1.fqmRNA.bam_alignment.stats   0.955742
SRR6230583_1_val_1.fqMult.bam_alignment.stats   0.955742
SRR6230584_1_val_1.fqGene.bam_alignment.stats   0.95413
SRR6230584_1_val_1.fqmRNA.bam_alignment.stats   0.95413
SRR6230584_1_val_1.fqMult.bam_alignment.stats   0.95413
SRR6230585_1_val_1.fqGene.bam_alignment.stats   0.949704
SRR6230585_1_val_1.fqmRNA.bam_alignment.stats   0.949704
SRR6230585_1_val_1.fqMult.bam_alignment.stats   0.949704
SRR6230586_1_val_1.fqGene.bam_alignment.stats   0.941641
SRR6230586_1_val_1.fqmRNA.bam_alignment.stats   0.941641
SRR6230586_1_val_1.fqMult.bam_alignment.stats   0.941641
SRR6230587_1_val_1.fqGene.bam_alignment.stats   0.956628
SRR6230587_1_val_1.fqmRNA.bam_alignment.stats   0.956628
SRR6230587_1_val_1.fqMult.bam_alignment.stats   0.956628
```
