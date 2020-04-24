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

### Mikado.loci.gff3 round2, including repeat genes
```
#gene length
less mikado.loci.gff3 |awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  120,922,982
Count:  38,646
Mean:   3,128
Median: 821
Min:    199
Max:    465,216

#transcript length
less mikado.loci.gff3 |awk '$3=="mRNA"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  142,902,933
Count:  45,286
Mean:   3,155
Median: 1,034
Min:    199
Max:    465,216

#exons per transcript
less mikado.loci.gff3 |awk '$3=="CDS" ' |cut -f 9 |sed 's/\./\t/3' |awk '{print $1}' |sort|uniq -c |awk '{print $1}' |summary.sh
Total:  258,802
Count:  45,286
Mean:   5
Median: 3
Min:    1
Max:    194

```

### mikado.loci.gff3
```
#exons per transcript
less mikado.loci.gff3 |awk '$3=="CDS"' |bedtools intersect -v -wo -a - -b ../../10_RepeatModeler/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/\./\t/3' |awk '{print $1}' |sort|uniq -c |awk '{print $1}' |summary.sh
Total:  228,567
Count:  31,566
Mean:   7
Median: 5
Min:    1
Max:    194

#how many genes are not repetitive
less mikado.loci.gff3 |awk '$3=="CDS"' |bedtools intersect -v -wo -a - -b ../../10_RepeatModeler/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/\./\t/2' |awk '{print $1}' |sort|uniq|wc
25180   25180  650224

#how many are repetitive then?
By subtraction that is: 13,466 genes

#gene sizes that are not repetitive
less NonRepetitiveGenes.gff3 |awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  104,488,696
Count:  25,180
Mean:   4,149
Median: 1,726
Min:    199
Max:    465,216

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


## Round2 of mikado using all ESTs from tylenchida
```
less mikado.loci.gff3| awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  131,514,741
Count:  39,516
Mean:   3,328
Median: 920
Min:    199
Max:    465,216
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
less 01_Interpro/interproAnnot.tsv |awk '{print $1}' |sort|uniq|wc
  25779   25779  639842

#proteins to uniprot
less 04_ProtsUniprot/mikado_proteins.vs.uniprot_sprot.cul5.1e5.blastp.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1 }' |sort|uniq|wc
   1604    1604   39727

#transcripts to uniref
less 05_TransUniprot/mikado_transcripts.vs.uniprot_sprot.cul5.1e5.blastx.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1 }' |sort|uniq|wc
  12796   12796  317987

#prots to nr
less 02_Prots2Nr/mikado_proteinsFixed.vs.
nr.cul5.1e5.blastp.out |grep -v "hypothetical" |grep -v "uncharacterized" |aw
k '{print $1}' |sort|uniq|wc
   3056    3056   75814

#transcripts to nt
ess 03_Transcrips2Nt/mikado_transcripts.vs.nt.cul5.1e5.blastn.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1}' |sort|uniq|wc
   2266    2266   56298



#All databases together
cat <(less 01_Interpro/interproAnnot.tsv |awk '{print $1}' |sort|uniq) <(less 02_Prots2Nr/mikado_proteinsFixed.vs.nr.cul5.1e5.blastp.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1}' |sort|uniq) <( less 03_Transcrips2Nt/mikado_transcripts.vs.nt.cul5.1e5.blastn.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1}' |sort|uniq) <(less 04_ProtsUniprot/mikado_proteins.vs.uniprot_sprot.cul5.1e5.blastp.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1 }' |sort|uniq) <(less 05_TransUniprot/mikado_transcripts.vs.uniprot_sprot.cul5.1e5.blastx.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1 }' |sort|uniq) |sort|uniq|wc
  26951   26951  668844

````
