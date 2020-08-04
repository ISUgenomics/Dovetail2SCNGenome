## Compile all annotations for gff


### NR Blast compilation
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/25_AnnotateGenes/02_Prots2Nr/SuperSplitter

cat *blastp.out |sort -k1,1V -u |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tNRBlastp:/g' >NRBlastp.tab
```
### NT Blastn annotation
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/25_AnnotateGenes/03_Transcrips2Nt
cat *blastn.out |sort -k1,1V -u |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tNTBlastn:/g' >NTBlastn.tab
```
### Uniprot Blastp
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/25_AnnotateGenes/04_ProtsUniprot
cat *blastp.out |sort -u -k1,1V |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tUPBlastp:/g'> UPBlastp.tab
```
### Uniprot Blastx
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/25_AnnotateGenes/05_TransUniprot
cat *blastx.out |sort -k1,1 -u |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tUPBlastx:/g' >UPBlastx.tab
```

### Interproscan annotation
```
```


### Interproscan annotation
```
#How many informative interproscan annotations are there?
cat *gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |wc
 157033 2259940 31882682
#How many mRNAs were annotated with informative annotations
cat *gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |sort -k1,1V -u |wc
19215  273421 3758676

#The fields are all different for each of these separations of the interproscan annotation, so 6 different files created
cat *gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") && NF==10) {print $1"\t"$9} else if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") && NF==11) {print $1"\t"$9","$11} else if(($2=="PANTHER" ||$2=="SMART" ||$2=="SUPERFAMILY") &&  NF==12) {print $1"\t"$8","$10","$12} else {next}}' >PantherSmartSuperfamily.tab

cat *gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if($2=="CDD" && NF==11) {print $1"\t"$9","$10} else if($2=="CDD" && NF==12) {print $1"\t"$9","$10","$13} if ($2=="CDD" && NF==13) {print $1"\t"$8","$10","$11","$13} else {next}}' >CDD.tab

cat *gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if(($2=="PIRSF" ||$2=="Gene3D" )&& NF==10) {print $1"\t"$9} else if (($2=="PIRSF"||$2=="Gene3D" ) && NF==11) {print $1"\t"$9","$11} else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==12) {print $1"\t"$10}else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==13) {print $1"\t"$9","$11} else if(($2=="Gene3D" ||$2=="PIRSF")  && NF==14) {print $1"\t"$9","$11} else {next}}' > Gene3dPirsf.tab

cat *gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if($2=="SFLD" || $2=="MobiDBLite") {print $1"\t"$10","$9} else if($2=="Coils") {print $1"\t"$9} else {next}}' >CoilsSfldMobidblite.tab

cat *gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if($2=="Pfam" &&NF==11) {print $1"\t"$9","$10} else if($2=="Pfam" &&NF==12) {print $1"\t"$9","$10","$12} else if($2=="Pfam" && NF==13) {print $1"\t"$8","$10","$12} else {next}}' >Pfam.tab

cat *gff3  |awk -F"\t" '$3!="polypeptide"' |awk -F"\t" '{if(substr($1,1,1)==">") {exit} else {print $0}}' |cut -f 1,2,4,5,7,9 |grep -v "#" |sed 's/;/\t/g' |awk -F"\t" '{if(($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap")&& NF==11 ) {print $1"\t"$9","$10} else if (($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap" )&& NF==12 ) {print $1"\t"$9","$10","$12} else if (($2=="PRINTS" ||$2=="ProDom" ||$2=="ProSiteProfiles" || $2=="TIGRFAM" ||$2=="ProSitePatterns" ||$2=="Hamap")&& NF==13 ) {print $1"\t"$8","$10","$11","$13} else {next}}' >PrintsProdomPrositepatternsPrositeprofilesTigrfamHamap.tab

#did I get all of the mRNA's?
cat *tab |sort -u -k1,1 |wc
  19215   55377 1437044
#yes

#did I get all of the annotations?
cat *tab |wc
157030  489415 14373210

#How many are exactly the same?
$ cat *tab |sort|uniq|wc
 105681  322508 9549698



```
### Combine all annotations
```

#get all the blasts together and interproscan together
awk -F"\t" '{arr[$1]=arr[$1] "\t" $2}END{for(i in arr)print i,arr[i]}' *tab  |sed 's/\t/#/1'|sed 's/\t/|/g' |sed 's/Name=//g'|sed 's/=//g' |sed 's/#/\tNote=/1'  >CombineAnnot.tab1

#How many mRNAs have an annotation?
sort -u -k1,1V CombineAnnot.tab1 |wc
  20322  514856 14300744

#How many genes have an annotation?
less CombineAnnot.tab1 |sed 's/\./\t/2' |cut -f 1,3 |sort -k1,1 -u |wc
 19136  482377 13406686

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
