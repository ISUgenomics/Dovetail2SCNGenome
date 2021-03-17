# Make functional annotations display in jbrowse



###  Create a list with annotated and unannoated genes
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/01_Transfer2Box

#Interpro is the most informative, place these annotations to display first
less FunctionalAnnotations.tab |sed 's/InterPro:/\tInterPro:/1' |sed 's/Note=/Note=\t/g' |awk  -F "\t" '{print $1,$2$4$3}' >ReformatFunctionalAnnotations.tab

less OrderedSCNGenePredictions.gff3 |awk '$3=="mRNA"' |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |cat ReformatFunctionalAnnotations.tab - |sort -k1,1 -u >AllGenesandFunctionalAnnotations.tab


#did this do as expected?
less OrderedSCNGenePredictions.gff3 |awk '$3=="mRNA"' |wc
  23933  215397 2908981
wc FunctionalAnnotations.tab
     20322   514856 14101030 FunctionalAnnotations.tab
less OrderedSCNGenePredictions.gff3 |awk '$3=="mRNA"' |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |cat FunctionalAnnotations.tab - |sort -k1,1 -u |awk 'NF>1' |wc
      20322  514856 14101030

Yes, all are there
```

### Attach this list to the GFF
```
/work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/01_Transfer2Box

#adds annotation to the mRNA
paste <(awk '$3=="mRNA"' OrderedSCNGenePredictions.gff3 |sort -k9,9V) <(sort -k1,1V AllGenesandFunctionalAnnotations.tab|cut -f 2 )  |sed 's/\t/;/9' >mrnsFunctional.gff3
paste <(awk '$3=="gene"' OrderedSCNGenePredictions.gff3 |sort -k9,9V) <(sed 's/\./\t/1' AllGenesandFunctionalAnnotations.tab|sort -k1,1 -u  |sed 's/ /\t/1' |cut -f 3)|sed 's/\t/;/9'  >geneFunctional.gff3


cat <(awk '$3!="mRNA" && $3!="gene"' OrderedSCNGenePredictions.gff3)  mrnsFunctionalSorted.gff3  geneFunctional.gff3 |sed 's/; /;/g' >CombineTest.gff3


perl ../gff3sort/gff3sort.pl --precise --chr_order natural CombineTest.gff3  >CombineTestFunctionalSorted.gff3


```
