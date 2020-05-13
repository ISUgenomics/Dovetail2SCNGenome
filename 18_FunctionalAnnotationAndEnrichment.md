## Compile all annotations for gff


### NR Blast compilation
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/25_AnnotateGenes/06_Combine
cat ../02_Prots2Nr/SuperSplitter/*blastp.out |sort -k1,1V -u |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tNRBlastp:/g' >NRBlastp.tab
```
### NT Blastn annotation
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/25_AnnotateGenes/06_Combine
cat ../03_Transcrips2Nt/OrderedSCNGenePredictions_VHEJtranscripts.part*blastn.out |sort -k1,1V -u |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tNTBlastn:/g' >NTBlastn.tab
```
### Uniprot Blastp
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/25_AnnotateGenes/06_Combine
cat ../04_ProtsUniprot/OrderedSCNGenePredictionsVHEJ_proteins.part*out |sort -u -k1,1V |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tUPBlastp:/g'> UPBlastp.tab
```
### Uniprot Blastx
```
cat ../05_TransUniprot/Ordered*blastx.out |sort -k1,1 -u |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tUPBlastx:/g' >UPBlastx.tab
```

### Interproscan annotation
```
#had to rerun this on the new annotation 4/9/20

less  ../01_Interpro/interproAnnot_1.gff3 |awk 'NR<282848 && $3!="polypeptide"' |cut -f 1,2,9- |sed 's/;/\t/3' |cut -f 1,2,4 |awk 'substr($1,1,1)!="#"' |sed 's/\t/:/2' |uniq |awk '{arr[$1]=arr[$1] "\t" $2}END{for(i in arr)print i,arr[i]}' |sed 's/\t/#/1' |sed 's/\t/;/g' |sed 's/#/\t/1' |sed 's/\t/\tIPR:/1'  >interproAnnot.tab1    
```
### Combine all annotations
```
#Fix the oddball annotation of transcript in $3
 awk -F"\t" '{if($3=="transcript"){print $1"\t"$2"\tmRNA\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9} else {print $0}}' OrderedSCNGenePredictions.gff3 >OddTranscriptFixOrderedSCNGenePredictions.gff3


#get all the blasts together
awk -F"\t" '{arr[$1]=arr[$1] "\t" $0}END{for(i in arr)print i,arr[i]}' *tab |cut -f 1,3,5,7,9,11  |sed 's/\t/#/1'|sed 's/\t/|/g' |sed 's/=//g'|sed 's/#/\tNote=/1' >CombineAnnot.tab1

#add in the interpro annotations
awk '{arr[$1]=arr[$1] "\t" $0}END{for(i in arr)print i,arr[i]}' *tab1 |sed 's/ ;/\t/1' |cut -f 2,3,5 |sed 's/\t/;/2' >CombeinAnnotIPRs.tab2


less ../07_NewGenes/OddTranscriptFixOrderedSCNGenePredictions.gff3 |awk '$3=="mRNA"' |sed 's/ID=/ID=\t/g' |sed 's/;/\t;/1' >OddTranscriptFixOrderedSCNGenePredictionsgff3.GREPMOD

#find the mrnas lacking annotations
awk '{print $1}' CombeinAnnotIPRs.tab2 |cat - <(cut -f 10 OddTranscriptFixOrderedSCNGenePredictionsgff3.GREPMOD) |sort|uniq -c |awk '$1==1 {print $2}' |cat - CombeinAnnotIPRs.tab2 >AllGenesWWOannot.tab3

#get the mrnas in the gff in the proper order and paste
paste <(sort -k10,10 OddTranscriptFixOrderedSCNGenePredictionsgff3.GREPMOD|sed 's/;Name=.*//g') <(sort -k1,1 AllGenesWWOannot.tab3) |awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9$10$11";"$13}' >AllGenesAnnotated.gff


#combine with exons,genes, utrs, cds, etc annotations
awk '$3!="mRNA"' ../07_NewGenes/OddTranscriptFixOrderedSCNGenePredictions.gff3 |sed 's/Name=.*//g'|grep -v "#" |cat - AllGenesAnnotated.gff >AllTypesWWOAnnotationsDisorganized.gff

#get the proper order for the gff, so it will show up in jbrowse
 perl gff3sort/gff3sort.pl --precise --chr_order natural AllTypesWWOAnnotationsDisorganized.gff > SCNgenomeFunctionalGeneAnnotations.gff3
```


#### Summary of mrnas that were annotated
```
#How many mRNAs were annotated?
less SCNgenomeFunctionalGeneAnnotations.gff3 |awk '$3=="mRNA"' |grep "Note" |sed 's/;/\t/1' |cut -f 9 |sed 's/\./\t/2' |cut -f 1 |sort|uniq|wc
  18942   18942  240910


#How many are there total again?
awk '$3=="mRNA"' SCNgenomeFunctionalGeneAnnotations.gff3|wc
  39481  730635 21161939


```

### Transfer annotations to gene for display in jbrowse
```
#  Will just get the top isoform of each mrna

#make grepmod of all gene annotations from
less SCNgenomeFunctionalGeneAnnotations.gff3 |awk '$3=="gene" ' |sed 's/ID=/ID=\t/g' |sed 's/;/\t;/1' >SCNgenomeFunctionalGeneAnnotationsgff3.GeneGREPMOD


#I will make a matrix of mRNA to gene names, so I can swap the mrna names out with gene names to make a tab file for genes to annotations.
 less SCNgenomeFunctionalGeneAnnotations.gff3 |awk '$3=="mRNA" {print $9}' |sed 's/Parent=//g' |sed 's/ID=//g' |sed 's/\;/\t/g' |cut -f 1,2 >mRNA2Gene.list

#Swap names in tab file, ditch ncrnas
paste <( sort -k1,1 mRNA2Gene.list) <(sort -k1,1 AllGenesWWOannot.tab3|awk -F"\t" '{if(NF>1){print $0} else{print $1"\t"}}' ) <( sort -k1,1 mRNA2Gene.list) |cut -f 2,4 |sort -k2,2Vr|sort -k1,1 -u  |awk 'substr($1,1,2)!="nc"' >GeneAnnotations.tabular


#Get the genes in gff format
paste <(sort -k10,10V SCNgenomeFunctionalGeneAnnotationsgff3.GeneGREPMOD|sed 's/Name=.*//g') <(less GeneAnnotations.tabular|sort -k1,1V) |cut -f 1-10,13 |awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9$10";"$11}' >GeneOnlyAnnot.gff3

awk '$3!="gene"' SCNgenomeFunctionalGeneAnnotations.gff3 |cat - GeneOnlyAnnot.gff3 >AllTypesWWOAnnotationsGenesMrnasDisorganized.gff

#get the proper order for the gff, so it will show up in jbrowse
 perl gff3sort/gff3sort.pl --precise --chr_order natural AllTypesWWOAnnotationsGenesMrnasDisorganized.gff > SCNgenomeFunctionalGeneAndMrnaAnnotations.gff3


#format for jbrowse #some of the genes look funny
sh ~/common_scripts/runTabix.sh SCNgenomeFunctionalGeneAndMrnaAnnotations.gff3








#How many genes are there total again?
awk '$3=="gene"' SCNgenomeFunctionalGeneAnnotations.gff3 |wc
  34454  310086 4984494

```
