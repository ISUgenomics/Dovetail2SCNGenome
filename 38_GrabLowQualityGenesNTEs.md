## Create jbrowse track and downloads for low quality genes

###  Repeat Overlap of 30% of 2 CDS from a single transcript
```
/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn/05_FinalGenePrediction

#could have used the created file, but just used the same command instead.
bedtools intersect -f .3 -wo -a <(awk '$3=="CDS"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/Parent=//g' |sed 's/;/\ -f 1 |sort|uniq -c |awk '$1>1{print $2}' |grep -w -f - <(sed 's/ID=/ID=\t/g' SCNGenePredictions.gff3|sed 's/Parent=/Parent=\t/g'|sed 's/;/\t/1') >RepeatsByOverlap.gff3

#need changed back to normal gff
```

### TE's by annotation
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn/05_FinalGenePrediction

grep -w -f GenesWithTransposonAnnotations.list <(sed 's/ID=/ID=\t/g' SCNGenePredictions.gff3|sed 's/Parent=/Parent=\t/g'|sed 's/;/\t/1') >GenesWithTransposonAnnotations.gff3 &

```

### Genes lacking expression
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn/05_FinalGenePrediction

grep -w -f <(awk '{print $1}' NoExpressedMrnaGene.list) <(sed 's/ID=/ID=\t/g' SCNGenePredictions.gff3|sed 's/Parent=/Parent=\t/g'|sed 's/;/\t/1') >NoExpressedMrnaGene.gff3 &

```


### Move to renaming folder
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes

cp ../01_mikadoRerurn/05_FinalGenePrediction/RepeatsByOverlap.gff3 .
cp ../01_mikadoRerurn/05_FinalGenePrediction/GenesWithTransposonAnnotations.gff3 .
cp ../01_mikadoRerurn/05_FinalGenePrediction/NoExpressedMrnaGene.gff3 .

ml maker
 for f in  RepeatsByOverlap.gff3 GenesWithTransposonAnnotations.gff3 NoExpressedMrnaGene.gff3 ; do map_data_ids chromosome.map $f;done
