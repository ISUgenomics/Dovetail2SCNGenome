# Generate differential expression for gland and whole worm RNASEQ


# generate bigwig and count files
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate

ln -s ../../25_AnnotateGenes/SCNgenomeFunctionalGeneAn
ln -s ../../10_RepeatModeler/SCNgenome.fasta
for f in  ../../11_AlignRNA/*sorted.bam; do ln -s $f;done



for f in *bam; do echo "ml samtools; samtools index "$f"; ml subread; featureCounts -T 16 -p -t gene -s 2 -g ID -a SCNgenomeFunctionalGeneAnnotationsHighConfidenceAllFeatures.gff3  -o "${f%.*}"_GeneCounts.txt " $f"; ml bedtools2; bedtools genomecov -ibam "$f" -bga -g SCNgenome.fasta > "${f%.*}"_sorted.bam.bdg; ml bioawk; bioawk -c fastx '{print \$name,length(\$seq)}'  SCNgenome.fasta >Chr.sizes ; /home/remkv6/common_scripts/kentUtils/bin/linux.x86_64/bedGraphToBigWig "${f%.*}"_sorted.bam.bdg Chr.sizes "${f%.*}"_sorted.bw"; done
```


### Deseq setup for All whole worm vs All Gland
```
#creates the command below
ls -1 *txt |while read line; do echo "<(awk '{print \$1,\$7}' "$line")"; done |tr "\n" " " |sed 's/^/paste /g' |sed 's/$/ |less/g' |less

#gets all the samples in the correct order
paste <(awk '{print $1,$7}' 1703FL-02-01_S1_L001_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703FL-02-02_S2_L001_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703FL-02-03_S3_L001_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703FL-02-04_S4_L001_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703-TM101_S0_L003_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703-TM102_S0_L003_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230579_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230580_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230581_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230582_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230583_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230584_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230585_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230586_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230587_1_val_1.fq_sorted_GeneCounts.txt)  |less

#createst the header to the deseq table
paste <(awk '{print $1,$7}' 1703FL-02-01_S1_L001_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703FL-02-02_S2_L001_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703FL-02-03_S3_L001_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703FL-02-04_S4_L001_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703-TM101_S0_L003_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703-TM102_S0_L003_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230579_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230580_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230581_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230582_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230583_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230584_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230585_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230586_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230587_1_val_1.fq_sorted_GeneCounts.txt)   |awk 'NR==2 {print $1,$14,$16,$18,$20,$22,$24,$26,$28,$30,"X"$2,"X"$4,"X"$6,"X"$8,"X"$10,"X"$12}' |tr " " "\t" |sed 's/_val_1\.fq_sorted.bam//g' >header.txt

Creates all the counts for the deseq table
paste <(awk '{print $1,$7}' 1703FL-02-01_S1_L001_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703FL-02-02_S2_L001_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703FL-02-03_S3_L001_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703FL-02-04_S4_L001_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703-TM101_S0_L003_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' 1703-TM102_S0_L003_R1_001_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230579_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230580_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230581_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230582_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230583_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230584_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230585_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230586_1_val_1.fq_sorted_GeneCounts.txt) <(awk '{print $1,$7}' SRR6230587_1_val_1.fq_sorted_GeneCounts.txt)   |awk 'NR>2 {print $1,$14,$16,$18,$20,$22,$24,$26,$28,$30,$2,$4,$6,$8,$10,$12}' |tr " " "\t" >tailer.txt

cat header.txt tailer.txt >DeseqTable.txt

 tr "\t" "\n" <header.txt |awk '{if(NR==1){print $0} else if(NR<11) {print $0"\t0"} else {print $0"\t1"}}'>ConditionTable.txt




ml r-deseq2/1.20.0-py2-r3.5-openmpi3-zhebatp


 library("DESeq2")
 dat<-read.table("DeseqTable.txt",header = T,quote = "",row.names = 1)
 dat <- as.matrix(dat)
 condition <- factor(c(rep("WT",9),rep("Mut",6)))
 condition=relevel(condition,ref = "WT")
 coldata <-read.table("ConditionTable.txt",header = T,row.names = 1)
 coldata <- data.frame(row.names=colnames(dat), condition)
 dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata,design=~ condition)
 dds <- DESeq(dds)
 res <- results(dds)
 table(res$padj<0.05)
 res <- res[order(res$padj), ]
 resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
 names(resdata)[1] <- "Gene"
 write.csv(resdata, file="AllGlandvsAllWholeWorm",quote = FALSE,row.names = F)

 estimating size factors
 estimating dispersions
 gene-wise dispersion estimates
 mean-dispersion relationship
 final dispersion estimates
 fitting model and testing
 -- replacing outliers and refitting for 1065 genes
 -- DESeq argument 'minReplicatesForReplace' = 7
 -- original counts are preserved in counts(dds)
 estimating dispersions
 fitting model and testing


 FALSE  TRUE
  8012  8514




  rld <- rlogTransformation(dds)
  library(ggplot2)
  p <- plotPCA(rld)
  p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
  print(p)
  q()
```

### Compare MM10 GLAND AND PA3 Gland
```
less DeseqTable.txt |cut -f 1,11,12,13,14,15,16 |awk '{print $1,$2,$6,$7,$3,$4,$5}' |tr " " "\t"  >PA3vsMM10GlandDeseqTable.txt


 awk 'NR==1' PA3vsMM10GlandDeseqTable.txt |tr "\t" "\n"|awk '{if(NR==1){print $0} else if(NR<5) {print $0"\t0"} else {print $0"\t1"}}' >PA3vsMM10GlandConditionTable.txt

 library("DESeq2")
 dat<-read.table("PA3vsMM10GlandDeseqTable.txt",header = T,quote = "",row.names = 1)
 dat <- as.matrix(dat)
 condition <- factor(c(rep("WT",3),rep("Mut",3)))
 condition=relevel(condition,ref = "WT")
 coldata <-read.table("PA3vsMM10GlandConditionTable.txt",header = T,row.names = 1)
 coldata <- data.frame(row.names=colnames(dat), condition)
 dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata,design=~ condition)
 dds <- DESeq(dds)
 res <- results(dds)
 table(res$padj<0.05)
 res <- res[order(res$padj), ]
 resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
 names(resdata)[1] <- "Gene"
 write.csv(resdata, file="PA3GlandvsMM10Gland",quote = FALSE,row.names = F)

 estimating size factors
 estimating dispersions
 gene-wise dispersion estimates
 mean-dispersion relationship
 final dispersion estimates
 fitting model and testing


 FALSE  TRUE
  1687     5


 rld <- rlogTransformation(dds)
 library(ggplot2)
 p <- plotPCA(rld)
 p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
 print(p)
 q()
```

### Compare PA3 Gland to PA3 parasitic compatible whole worm
```
less DeseqTable.txt |cut -f 1,6,7,10,11,15,16 >pPA3CvsPA3GlandDeseqTable.txt
awk 'NR==1' pPA3CvsPA3GlandDeseqTable.txt |tr "\t" "\n"|awk '{if(NR==1){print $0} else if(NR<5) {print $0"\t0"} else {print $0"\t1"}}' >pPA3CvsPA3GlandConditionTable.txt

library("DESeq2")
dat<-read.table("pPA3CvsPA3GlandDeseqTable.txt",header = T,quote = "",row.names = 1)
dat <- as.matrix(dat)
condition <- factor(c(rep("WT",3),rep("Mut",3)))
condition=relevel(condition,ref = "WT")
coldata <-read.table("pPA3CvsPA3GlandConditionTable.txt",header = T,row.names = 1)
coldata <- data.frame(row.names=colnames(dat), condition)
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata,design=~ condition)
dds <- DESeq(dds)
res <- results(dds)
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
write.csv(resdata, file="pPA3CvsPA3Gland",quote = FALSE,row.names = F)

estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing


FALSE  TRUE
 5212   912

```

### Compare PA3 preparasitic J2 to parasitic J2 compatible
```
less DeseqTable.txt |awk  '{print $1,$2,$9,$5,$6,$7,$10}' |tr " " "\t" >ppPA3vspPA3CDeseqTable.txt

less DeseqTable.txt |awk  '{print $1,$2,$4,$5,$6,$7,$10}' |tr " " "\t" >ppPA3vspPA3CDeseqTable.txt
 awk 'NR==1' ppPA3vspPA3CDeseqTable.txt |tr "\t" "\n"|awk '{if(NR==1){print $0} else if(NR<5) {print $0"\t0"} else {print $0"\t1"}}' >ppPA3vspPA3CConditionTable.txt

library("DESeq2")
dat<-read.table("ppPA3vspPA3CDeseqTable.txt",header = T,quote = "",row.names = 1)
dat <- as.matrix(dat)
condition <- factor(c(rep("WT",3),rep("Mut",3)))
condition=relevel(condition,ref = "WT")
coldata <-read.table("ppPA3vspPA3CConditionTable.txt",header = T,row.names = 1)
coldata <- data.frame(row.names=colnames(dat), condition)
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata,design=~ condition)
dds <- DESeq(dds)
res <- results(dds)
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
write.csv(resdata, file="ppPA3vspPA3C",quote = FALSE,row.names = F)

estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing


FALSE  TRUE
 2268    16


```

### Compare PA3 preparasitic J2 to parasitic J2 incompatible
```
less DeseqTable.txt |awk '{print  $1,$2,$9,$5,$3,$8,$4}' |tr " " "\t" >ppPA3vspPA3ICDeseqTable.txt
awk 'NR==1' ppPA3vspPA3ICDeseqTable.txt |tr "\t" "\n"|awk '{if(NR==1){print $0} else if(NR<5) {print $0"\t0"} else {print $0"\t1"}}' >ppPA3vspPA3ICConditionTable.txt

library("DESeq2")
dat<-read.table("ppPA3vspPA3ICDeseqTable.txt",header = T,quote = "",row.names = 1)
dat <- as.matrix(dat)
condition <- factor(c(rep("WT",3),rep("Mut",3)))
condition=relevel(condition,ref = "WT")
coldata <-read.table("ppPA3vspPA3ICConditionTable.txt",header = T,row.names = 1)
coldata <- data.frame(row.names=colnames(dat), condition)
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata,design=~ condition)
dds <- DESeq(dds)
res <- results(dds)
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
write.csv(resdata, file="ppPA3vspPA3IC",quote = FALSE,row.names = F)

estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing

FALSE  TRUE
 6971 10720


```

### Compare PA3 compatible parasitic to PA3 incompatible parasitic
```
less DeseqTable.txt |awk '{print $1,$6,$7,$10,$3,$8,$4}' |tr " " "\t" >pPA3CvspPA3ICDeseqTable.txt
awk 'NR==1' pPA3CvspPA3ICDeseqTable.txt |tr "\t" "\n"|awk '{if(NR==1){print $0} else if(NR<5) {print $0"\t0"} else {print $0"\t1"}}' >pPA3CvspPA3ICConditionTable.txt

library("DESeq2")
dat<-read.table("pPA3CvspPA3ICDeseqTable.txt",header = T,quote = "",row.names = 1)
dat <- as.matrix(dat)
condition <- factor(c(rep("WT",3),rep("Mut",3)))
condition=relevel(condition,ref = "WT")
coldata <-read.table("pPA3CvspPA3ICConditionTable.txt",header = T,row.names = 1)
coldata <- data.frame(row.names=colnames(dat), condition)
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata,design=~ condition)
dds <- DESeq(dds)
res <- results(dds)
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
write.csv(resdata, file="pPA3CvspPA3IC",quote = FALSE,row.names = F)


estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing

FALSE  TRUE
 6159 10353


```

### Compare 2 PA3 rep to 3 MM10 rep
```
less DeseqTable.txt |cut -f 1,11,12,13,14,15,16 |awk '{print $1,$6,$7,$3,$4,$5}' |tr " " "\t"  >2PA3vs3MM10GlandDeseqTable.txt


 awk 'NR==1' 2PA3vs3MM10GlandDeseqTable.txt |tr "\t" "\n"|awk '{if(NR==1){print $0} else if(NR<4) {print $0"\t0"} else {print $0"\t1"}}' >2PA3vs3MM10GlandConditionTable.txt

 library("DESeq2")
 dat<-read.table("2PA3vs3MM10GlandDeseqTable.txt",header = T,quote = "",row.names = 1)
 dat <- as.matrix(dat)
 condition <- factor(c(rep("WT",2),rep("Mut",3)))
 condition=relevel(condition,ref = "WT")
 coldata <-read.table("2PA3vs3MM10GlandConditionTable.txt",header = T,row.names = 1)
 coldata <- data.frame(row.names=colnames(dat), condition)
 dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata,design=~ condition)
 dds <- DESeq(dds)
 res <- results(dds)
 table(res$padj<0.05)
 res <- res[order(res$padj), ]
 resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
 names(resdata)[1] <- "Gene"
 write.csv(resdata, file="2PA3vs3MM10GlandGland",quote = FALSE,row.names = F)


 FALSE  TRUE
  2157    20

 rld <- rlogTransformation(dds)
 library(ggplot2)
 p <- plotPCA(rld)
 p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
 print(p)
 q()




 less ../39_DifferentialExpression/2PA3vs3MM10GlandGland |sed 's/,/\t/g' |awk 'NR<7284 && NR>1 && $3<0{print $1}' >2PA3vs3MM10GlandListDown.tab
 less ../39_DifferentialExpression/2PA3vs3MM10GlandGland |sed 's/,/\t/g' |awk 'NR<7284 && NR>1 && $3>0{print $1}' >2PA3vs3MM10GlandListUp.tab

```

### Diff Expression 2 MM10's and 2 PA3's
```
less DeseqTable.txt |cut -f 1,11,12,13,14,15,16 |awk '{print $1,$6,$7,$4,$5}' |tr " " "\t" >2PA3vs2MM10GlandDeseqTable.txt



 awk 'NR==1' 2PA3vs2MM10GlandDeseqTable.txt |tr "\t" "\n"|awk '{if(NR==1){print $0} else if(NR<4) {print $0"\t0"} else {print $0"\t1"}}' >2PA3vs2MM10GlandConditionTable.txt

 library("DESeq2")
 dat<-read.table("2PA3vs2MM10GlandDeseqTable.txt",header = T,quote = "",row.names = 1)
 dat <- as.matrix(dat)
 condition <- factor(c(rep("WT",2),rep("Mut",2)))
 condition=relevel(condition,ref = "WT")
 coldata <-read.table("2PA3vs2MM10GlandConditionTable.txt",header = T,row.names = 1)
 coldata <- data.frame(row.names=colnames(dat), condition)
 dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata,design=~ condition)
 dds <- DESeq(dds)
 res <- results(dds)
 table(res$padj<0.05)
 res <- res[order(res$padj), ]
 resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
 names(resdata)[1] <- "Gene"
 write.csv(resdata, file="2PA3vs2MM10GlandGland",quote = FALSE,row.names = F)

 FALSE  TRUE
  2669    79



 less ../39_DifferentialExpression/2PA3vs2MM10GlandGland |sed 's/,/\t/g' |awk 'NR<51 && NR>1 && $3<0{print $1}' >2PA3vs2MM10GlandListDown.tab
 less ../39_DifferentialExpression/2PA3vs2MM10GlandGland |sed 's/,/\t/g' |awk 'NR<51 && NR>1 && $3>0{print $1}' >2PA3vs2MM10GlandListUp.tab

```
