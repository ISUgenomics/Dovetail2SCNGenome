
### Investigate genes, How many are clearly transposable elements?
```

#How many mRNA's in the gene prediction overlap with repetitive elements?

#60% overlap
(mikado) [remkv6@nova023 05_FinalGenePrediction]$ bedtools intersect -f .6 -wo -a <(awk '$3=="mRNA"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
   5150    5150   49819
#50% overlap
(mikado) [remkv6@nova023 05_FinalGenePrediction]$ bedtools intersect -f .5 -wo -a <(awk '$3=="mRNA"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
   6133    6133   59294
#30% overlap
(mikado) [remkv6@nova023 05_FinalGenePrediction]$ bedtools intersect -f .3 -wo -a <(awk '$3=="mRNA"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
   8700    8700   84152
#20% overlap   
(mikado) [remkv6@nova023 05_FinalGenePrediction]$ bedtools intersect -f .2 -wo -a <(awk '$3=="mRNA"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
  10813   10813  104663

#Any intersect
(mikado) [remkv6@nova023 05_FinalGenePrediction]$ bedtools intersect  -wo -a <(awk '$3=="mRNA"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
  25929   25929  251739


#What if we consider only CDS overlaps with repeats
#Any overlap between any cds and any repeat
bedtools intersect -f .5 -wo -a <(awk '$3=="CDS"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/Parent=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
14159   14159  137011

# 50% of the CDS length, only one CDS required
bedtools intersect -f .5 -wo -a <(awk '$3=="CDS"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/Parent=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
  11791   11791  114040

# requires at least 2 CDS to overlap any part of repeat
bedtools intersect -f .5 -wo -a <(awk '$3=="CDS"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/Parent=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq -c |awk '$1>1' |wc
  9775   19550  172686
#requires at least 30% overlap of cds by the repeat, and two CDS's
bedtools intersect -f .3 -wo -a <(awk '$3=="CDS"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/Parent=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq -c |awk '$1>1' |wc
   8313   16626  146776

#Probably the most conservative measure, same as above
bedtools intersect -f .3 -wo -a <(awk '$3=="CDS"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/Parent=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq -c |awk '$1>1 {print $2}' >RepetitiveBy30PercOverlap2CDS.list

#How many of the mRNAs have an annotation that is transposon, helitron, or transposase?
less SCNgenomeFunctionalGeneAndMrnaAnnotations_sorted.gff |awk '$3=="mRNA"' |cut -f 9 |grep -i -e  "transposon" -e "transposase" -e  "helitron"  |wc
643   13971  513328

less SCNgenomeFunctionalGeneAndMrnaAnnotations_sorted.gff |awk '$3=="mRNA"' |cut -f 9 |grep -i -e  "transposon" -e "transposase" -e  "helitron" >GenesWithTransposonAnnotations.list
```

# Are all of the effectors there?
```

#identify by name only "gland" and "effector"
bedtools intersect -wo -f .5 -r -a <(awk '$3=="mRNA"'  SCNGenePredictions.gff3)  -b <(awk '$3=="mRNA"'  ../OldPredictionSCNgenomeFunctionalGeneAnnotations.gff3 ) |grep -i  "effector" |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort |uniq >effector.list
bedtools intersect -wo -f .5 -r -a <(awk '$3=="mRNA"'  SCNGenePredictions.gff3)  -b <(awk '$3=="mRNA"'  ../OldPredictionSCNgenomeFunctionalGeneAnnotations.gff3 ) |grep -i  "gland" |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort |uniq >>glandGenes.list
cat glandGenes.list effector.list |sort|uniq >EffectorGenes.list
wc EffectorGenes.list
 396  396 3780 EffectorGenes.list


#how many of these effectors are repetitive as called by repeatmasker
cat <(sort EffectorGenes.list |uniq) <(sort RepetitiveBy30PercOverlap2CDS.list |uniq) |sort|uniq -c |awk '$1==2' |wc
    176     352    3031

```

### Examine the expression of genes, leaving repetitive elements
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate

#zero expression
less DeseqTable.txt |awk '{print $1,$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16}' |awk '$2=="0"' |wc
  10379   20758  121089
#one read expression  
less DeseqTable.txt |awk '{print $1,$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16}' |awk '$2=="1"' |wc
   1981    3962   23137

# make a list of all mrna and gene names, grep to get them and get new gene count
awk '$3=="mRNA" || $3=="transcript"' OrderedSCNGenePredictions.gff3 |cut -f 9 |sed 's/ID=//g' |sed 's/Parent=//g' |sed 's/;/\t/g' |cut -f 1,2 >mrnaGeneAll.list

#make gene list for mikado to remove from the gff
less DeseqTable.txt |awk '{print $1,$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16}' |awk '$2=="0" {print $1}' |grep -w -f - mrnaGeneAll.list >NoExpressedMrnaGene.list
wc NoExpressedMrnaGene.list
 11268  22536 218395 NoExpressedMrnaGene.list

#remove these from the annotation

```


# List of genes to keep for counts.  i.e. not expressed or repeats.  
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn/05_FinalGenePrediction

awk '$3=="mRNA" ' SCNGenePredictions.gff3 |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 >AllmRNAs.list

#Remove the overlap with repeats, and annotated transposon proteins
#Figure in non expressed transcripts later when featurecounts is available
cat AllmRNAs.list GenesWithTransposonAnnotations.list RepetitiveBy30PercOverlap2CDS.list NoExpressedMrnaGene.list |awk '{print $1}' |sort|uniq -c |awk '$1==1{print $2}' |cat - EffectorGenes.list |sort|uniq -c |awk '$1==1{print $2}' >NoRepeatExpressedFinalGene.List

wc NoRepeatExpressedFinalGene.List
22921  22921 223262 NoRepeatExpressedFinalGene.List

grep -w -f NoRepeatExpressedFinalGene.List mrnaGeneAll.list >NoRepeatExpressedFinalGeneAndmRNA.List

#This places all exons cds and mrna next to the parent so no errors in the mikado grep
gt gff3 -checkids OrderedSCNGenePredictions.gff3 >test.gff

mikado util grep NoRepeatExpressedFinalGeneAndmRNADone.list test.gff  FilteredSCNGene
Predictions.gff3

#extract the functional annotations, noticed that this moves Name= before Note= and moves them to lower case
gt -q yes gff3 -checkids -tidy SCNgenomeFunctionalGeneAndMrnaAnnotations_sorted.gff >test.gff
mikado util grep NoRepeatExpressedFinalGeneAndmRNADone.list test.gff FunctionalFilteredSCNgenomeAnnotations.gff

```


### Merge good alignmetns with cufflinks gffread -M ### testing still 4/14/20
 While this does substantially lower the number of genes that are generated, it does combine genes with mrnas on two strands, which is unlikely.  Having them all displayed as a potential mRNAs may be helpful, but also display lots of spurious transcripts toois also
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn/03_MethionineProteinsOnly

#rename all gth alignments from each transcript assembler, so that gffread merge can tell the difference in the transcripts
ls */*aln |while read line; do awk '{print substr($1,1,5),$0}' $f;done |sed 's/AllRN/CLASS2/1' |sed 's/Brake/BRAKER/1' |sed 's/NCBIE/NCBIEST/1' |sed 's/SCNge/MAKER/1' |sed 's/Spade/SPADES/1' |sed 's/Strin/STRINGTIE/1' |sed 's/Trini/TRINITY/1' |while read line; do echo "sed -i 's/ID=/ID= "$line" /g' "$line" ; sed -i 's/Parent=/Parent=" $line "/g' ";done |awk '{print $1,$2,$3$4$6,$8,$9,$10,$11,$12$13$15,$14}' >RenameAllGenomethreaderAlignments.sh
sh RenameAllGenomethreaderAlignments.sh

#/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn/05_FinalGenePrediction

#merge with all genomethreader alignments of Met start encoded proteins and to All H.glycines EST alignments
gffread -E -K -M -Q -VHEJ -g ../SCNgenome.fasta  ../AllRNASEQClass2.gtf ../NonRrnaRNASEQ_stringtie.gtf ../SCNgenome.transcripts.gff3 ../SCNgenome.Trinity-GG.gff3 ../SCNgenome.DovetailSCNMaker4.all.maker.transcripts.gff3 ../brakermasked.gff ../tidiedBrakerUnmasked.gff3 ../mikado.loci.gff3  -o AllProtsMetESTMergeOrderedSCNGenePredictions.gff3 -x  MergingTestSCNVHEJ_transcripts.fasta -y  MergingTestSCNVHEJ_proteins.fasta

#How many genes are there if we merge with all transcripts starting with M from the genomethreader alignments of all genome-guided transcriptomes?
bedtools intersect -wo -f .3 -a <(awk '$3=="locus"' AllProtsMetMergeESTmergeMergeOrderedSCNGenePredictions.gff3) -b FixedaugustusreformatName.gff3 |awk '$12=="gene"' |cut -f 1-9  |wc
 12518  112662 1999026

#make a gff to look at in jbrowse, genometools throws a fit because the parent ID combos are not unique via the parallel gth alignments.
awk '{if($3=="locus") {print $1,$2,"gene",$4,$5,$6,$7,$8,$9} else {print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' AllProtsMetESTMergeOrderedSCNGenePredictions.gff3 |tr " " "\t" >MergingTestSCN.gff3


ml tabix/2013-12-16-py2-talztjl
bgzip MergingTestSCN.gff3
tabix -p gff MergingTestSCN.gff3.gz


# Continued 05/14/20, as the previous prediction did poorly on busco3 protein mode.
Lets convert to proteins and see how it does with busco protein mode

ml cufflinks
gffread -VHEJ MergingTestSCN.gff3  -g ../SCNgenome.fasta  -t mRNA -x  MergingTestSCNVHEJ_transcripts.fasta -y  MergingTestSCNVHEJ_proteins.fasta


grep -c ">" MergingTestSCNVHEJ_proteins.fasta
98024
58% busco complete, not good enough



```
