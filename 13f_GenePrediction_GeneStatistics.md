## Compare output with existing annotations
### Schachtii gene call stats
```
#gene length
less ../../../12_SchachtiiSynteny/H_sch_gene_calls_v1_CP.gff |awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  78,025,824
Count:  26,739
Mean:   2,918
Median: 2,129
Min:    252
Max:    56,948

#CDS length
(mikado) [remkv6@condofree032 02_Round02]$ less ../../../12_SchachtiiSynteny/H_sch_gene_calls_v1_CP.gff |awk '$3=="CDS"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  36,075,582
Count:  222,897
Mean:   161
Median: 125
Min:    1
Max:    13,080

#transcript length
less ../../../12_SchachtiiSynteny/H_sch_gene_calls_v1_CP.gff |awk '$3=="transcript"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  101,553,992
Count:  32,624
Mean:   3,112
Median: 2,275
Min:    252
Max:    56,948

#exons per transcript
less ../../../12_SchachtiiSynteny/H_sch_gene_calls_v1_CP.gff |awk '$3=="CDS" ' |cut -f 9 |awk '{print $1}' |sort |uniq -c |awk '{print $1}' |summary.sh
Total:  222,897
Count:  32,624
Mean:   6
Median: 5
Min:    1
Max:    113


```

### Filtered gene predictions
```
#gene length
less FilteredSCNGenePredictions.gff3 |awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  57,286,114
Count:  19,653
Mean:   2,914
Median: 1,859
Min:    200
Max:    324,453


#transcript length
less FilteredSCNGenePredictions.gff3 |awk '$3=="mRNA"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
less FilteredSCNGenePredictions.gff3 |awk '$3=="mRNA"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  66,883,744
Count:  22,533
Mean:   2,968
Median: 1,915
Min:    200
Max:    324,453


#exons per transcript
less FilteredSCNGenePredictions.gff3 |awk '$3=="CDS" ' |cut -f 9 |sed 's/\./\t/3' |awk '{print $1}' |sort|uniq -c |awk '{print $1}' |summary.sh
Total:  150,908
Count:  22,533
Mean:   6
Median: 5
Min:    1
Max:    209

less FilteredSCNGenePredictions.gff3 |awk '$3=="CDS"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  23,713,115
Count:  150,908
Mean:   157
Median: 125
Min:    0
Max:    10,527





```

### mikado.loci.gff3
```
#exons per transcript
less FilteredSCNGenePredictions.gff3 |awk '$3=="CDS"' |bedtools intersect -v -wo -a - -b ../../10_RepeatModeler/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/\./\t/3' |awk '{print $1}' |sort|uniq -c |awk '{print $1}' |summary.sh
Total:  228,567
Count:  31,566
Mean:   7
Median: 5
Min:    1
Max:    194



```
### gene stats of braker masked annotation
```
less augustus.hints.gff |awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  67,011,405
Count:  22,408
Mean:   2,990
Median: 2,015
Min:    200
Max:    52,212

```

### gene stats of braker unmasked annotation
```
less augustus.hints.gff |awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  83,144,932
Count:  35,514
Mean:   2,341
Median: 1,546
Min:    200
Max:    43,263
```

### gene stats of maker annotation on old dovetail genome
```
/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding

less 12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff|awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  67,168,685
Count:  22,856
Mean:   2,938
Median: 1,811
Min:    5
Max:    65,204

```



### gene calls from 738 genome
```

less ../CamTechGenomeComparison/58_Renamatorium/1_genomeNgff/fixed.augustus.gff3|awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  71,832,760
Count:  29,769
Mean:   2,413
Median: 1,603
Min:    79
Max:    65,717

```

### Functional annotation stats
```
#interproscan
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="gene"' |grep ";ipr:" |wc
  15896  333598 7927368
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="mRNA"' |grep ";ipr:" |wc
  17858  379630 9452960


#proteins to uniprot
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="gene"' |grep "|UPBlastp:" |wc
2669   95182 2433753
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="mRNA"' |grep "|UPBlastp:" |wc
3116  111134 2933213

#proteins to NR
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="gene"' |grep "NRBlastp:" |wc
2239   65191 1546985
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="mRNA"' |grep "NRBlastp:" |wc
2557   75387 1854103

#transcripts to uniprot
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="gene"' |grep "UPBlastx:" |wc
   8475  262490 6617663
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="mRNA"' |grep "UPBlastx:" |wc
   9714  301477 7862747

#transcripts to NT
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="gene"' |grep "NTBlastn:" |wc
   1917   60266 1471066
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="mRNA"' |grep "NTBlastn:" |wc
2163   69204 1778873





#All databases together
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="gene"' |grep -e  "note=" -e "ipr:" |wc
  16487  344745 8046194
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="mRNA"' |grep -e  "note=" -e "ipr:" |wc
  18459  390976 9584008

#Genes with PFAM annotations
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="gene"' |grep -i  "pfam" |wc
  8413  234661 6271921

#Genes with KEGG annotations
less FunctionalFilteredSCNgenomeAnnotations.gff |awk '$3=="gene"' |grep -i  "kegg" |wc
    707   19566  673143


#Genes with GO annotations

Hmm. interproscan didnt seem to add GO terms

````
