# Signal peptide prediction

### SignalP 5.0
```
wget https://services.healthtech.dtu.dk/download/66d87146-6d25-45b7-9bd9-62f21520bcbe/signalp-5.0b.Linux.tar.gz

tar -zxvf signalp-5.0b.Linux.tar.gz
cd signalp-5.0b/bin/

ln -s /work/GIF/remkv6/Baum/04_Dovetail2Restart/25_AnnotateGenes/07_NewGenes/OrderedSCNGenePredictionsVHEJ_proteins.fasta

./signalp -fasta OrderedSCNGenePredictionsVHEJ_proteins.fasta -gff3
less OrderedSCNGenePredictionsVHEJ_proteins_summary.signalp5 |grep -v "#" |awk '{print $1,$3}' |tr " " "\t" >SignalP5.tab

### how many??
less OrderedSCNGenePredictionsVHEJ_proteins_summary.signalp5 |grep -v "#" |awk '{print $1,$3}' |tr " " "\t" |awk '$2>.6' |wc
   3437    6874   64432
[remkv6@condo042 bin]$ less OrderedSCNGenePredictionsVHEJ_proteins_summary.signalp5 |grep -v "#" |awk '{print $1,$3}' |tr " " "\t" |awk '$2>.5' |wc
   3714    7428   69635
[remkv6@condo042 bin]$ less OrderedSCNGenePredictionsVHEJ_proteins_summary.signalp5 |grep -v "#" |awk '{print $1,$3}' |tr " " "\t" |awk '$2>.8' |wc
   2841    5682   53265
[remkv6@condo042 bin]$ less OrderedSCNGenePredictionsVHEJ_proteins_summary.signalp5 |grep -v "#" |awk '{print $1,$3}' |tr " " "\t" |awk '$2>.9' |wc
   2353    4706   44102

```

#### SignalP 4.1

```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/32_SignalP
 ln -s ../47_MikadoFinalize/mikado.loci.ancestralVHEJ_proteins.fasta

ml signalp/4.1f-py2-binsfdx

fasta-splitter.pl --n-parts 100 mikado.loci.ancestralVHEJ_proteins.fasta
for f in *part* ; do signalp -f short $f > ${f}.signalP4.out; done

 cat *signalP4.out |grep -v "#" |awk '{print $1"\t"$9}' >mikado.loci.ancestralVHEJ_proteinsSignalP4.tab

 awk '$2>.6' mikado.loci.ancestralVHEJ_proteinsSignalP4.tab |wc
    2317    4634   71334
```

### Signalp 5.0 Rerun with all transcripts assembled.  
```


./signalp -fasta MergingTestSCNVHEJ_proteins.fasta -gff3


#so there are some duplicate mrna names I guess, so there may be some flaws to using this
less MergingTestSCNOnly_proteins_summary.signalp5|awk '$2!="OTHER"' |cut -f 1 |awk 'NR>2'|wc
   8898    8898  114428
[remkv6@condo045 bin]$ less MergingTestSCNOnly_proteins_summary.signalp5|awk '$2!="OTHER"' |cut -f 1 |awk 'NR>2'|sort|uniq|wc
   7964    7964  100750


less MergingTestSCNgffreadSecretedOnly.gff3|awk  '$2=="Secreted"' |awk '{print $4,"Secreted",$6,$7,$8,$9,$10,$11,"ID="$3 }' |tr " " "\t" >mRNAMergingTestSCNgffreadSecretedOnly.gff3

bedtools intersect -wo -f .7 -a <(awk '$3=="mRNA"' ../../../../25_AnnotateGenes/06_Combine/SCNgenomeFunctionalGeneAndMrnaAnnotations.gff3) -b mRNAMergingTestSCNgffreadSecretedOnly.gff3   | cut -f 1-9 |sort |uniq|sort -k1,1V -k4,5nr >AnnotatedSecretedmRNAMergingTestSCNgffread.gff3

#what are the number of secreted heterodera avenae gland expressed?
bedtools intersect -wo -f .7 -a <(awk '$3=="mRNA"' ../../../../25_AnnotateGenes/06_Combine/SCNgenomeFunctionalGeneAndMrnaAnnotations.gff3) -b mRNAMergingTestSCNgffread
SecretedOnly.gff3  |cut -f 1-9 |sort|uniq|grep -c -i "heterodera avenae "
69
## and the total number of annotations
awk '$3=="mRNA"' ../../../../25_AnnotateGenes/06_Combine/SCNgenomeFunctionalGeneAndMrnaAnnotations.gff3 |grep -c -i "heterodera avenae "
88

#what are the number of secreted gmapped effectors?
bedtools intersect -wo -f .7 -a <(awk '$3=="mRNA"' ../../../../25_AnnotateGenes/06_Combine/SCNgenomeFunctionalGeneAndMrnaAnnotations.gff3) -b ../../../../29_Effectors/SCNgenome.effector_sorted.gff |cut -f 10-18 |awk '$3=="mRNA"' |sort|uniq|wc
     91     819   25744

#and the total number of gmapped effectors?
awk '$3=="mRNA"' ../../../../29_Effectors/SCNgenome.effector_sorted.gff|sort|uniq|wc                                                                                       126    1134   35026


```


### SignalP 3.0

```

less mikado.loci.ancestralVHEJ_proteins.part-01.signalp3.out |grep ">" -A 4 |grep "cleavage" -B 4 |grep -v "\-\-" |tr "\n" "\t" |sed 's/>/\n/g' |awk '{print $1,$8,$20}'|less

```
