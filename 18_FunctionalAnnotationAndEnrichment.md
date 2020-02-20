



## Compile all annotations for gff


### NR Blast compilation
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/25_AnnotateGenes/06_Combine
cat ../02_Prots2Nr/*blastp.out |sort -k1,1V -u |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tNRBlastp:/g' >NRBlastp.tab

```

### NT Blastn annotation
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/25_AnnotateGenes/06_Combine


less ../03_Transcrips2Nt/mikado_transcripts.vs.nt.cul5.1e5.blastn.out |sort -k1,1V -u |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tNTBlastn:/g' >NTBlastn.tab

```

### Uniprot Blastp
```
sort -u -k1,1V mikado_proteinsFixed.vs.uniprot_sprot.cul5.1e5.blastp.out |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tUPBlastp:/g'> ../06_Combine/UPBlastp.tab



```

### Uniprot Blastx
```
sort -k1,1 -u mikado_transcripts.vs.uniprot_sprot.cul5.1e5.blastx.out |cut -f 4,5,17 |sed 's/\t/:/2' |sed 's/\t/\tUPBlastx:/g' ../06_Combine/UPBlastx.tab




```

### Interproscan annotation
```
less  ../01_Interpro/interproAnnot.gff3 |awk 'NR<234320 && $3!="polypeptide"' |cut -f 1,2,9- |sed 's/;/\t/3' |cut -f 1,2,4 |awk 'substr($1,1,1)!="#"' |sed 's/\t/:/2' |uniq |awk '{arr[$1]=arr[$1] "\t" $2}END{for(i in arr)print i,arr[i]}' |sed 's/\t/#/1' |sed 's/\t/;/g' |sed 's/#/\t/1' |sed 's/\t/\tIPR:/1'  >interproAnnot.tab1
```
### Combine all annotations
```

#get all the blasts together
awk '{arr[$1]=arr[$1] ";" $2}END{for(i in arr)print i,arr[i]}' *tab |sed 's/ ;/\tNote=/1'   >CombineAnnot.tab1

#add in the interpro annotations
awk '{arr[$1]=arr[$1] ";" $2}END{for(i in arr)print i,arr[i]}' *tab1 |sed 's/ ;/\t/1' >CombeinAnnotIPRs.tab2

#find the mrnas lacking annotations
awk '{print $1}' CombeinAnnotIPRs.tab2 |cat - <(cut -f 10 Uniquemikado.grepmod.gff) |sort|uniq -c |awk '$1==1 {print $2 }' |cat - CombeinAnnotIPRs.tab2 >AllGenesWWOannot.tab3

#get the mrnas in the gff in the proper order for a paste below
less AllGenesWWOannot.tab3 |awk '{print $1}' |while read line; do grep -m 1 -w $line mikado.grepmod.gff >>AllGenesWWOannot.gffmod; done

#get all mrnas in gff format
less AllGenesWWOannot.gffmod |sed 's/\t//9' |sed 's/\t/;/9' |sed 's/$/;/g' |paste - <(awk '{print $2}'  AllGenesWWOannot.tab3) |sed 's/\t//9' >AllGenesAnnotated.gff

#combine with exons,genes, utrs, cds, etc annotations
awk '$3!="mRNA"' ../mikado.loci.gff3 |grep -v "#" |cat - AllGenesAnnotated.gff >AllTypesWWOAnnotationsDisorganized.gff

#get the proper order for the gff, so it will show up in jbrowse
 perl gff3sort/gff3sort.pl --precise --chr_order natural AllTypesWWOAnnotationsDisorganized.gff > SCNgenomeFunctionalGeneAnnotations.gff3

#How many genes were annotated?

 less SCNgenomeFunctionalGeneAnnotations.gff3 |awk '$3=="mRNA"' |grep "Note" |sed 's/;/\t/1' |cut -f 9 |sed 's/\./\t/2' |cut -f 1 |sort|uniq|wc
   14900   14900  384818


```
