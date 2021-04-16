# Need a comparison between effectors among X12 and this TN10 assembly

### gather resources
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/58_X12_Effectors

ln -s ../01_CondoPorts/28_X12GenomeComparison/SCN_genome.fa
ln -s ../01_CondoPorts/28_X12GenomeComparison/pasa2.longest.filter.pep
wget https://scnbase.org/files/download/effector.fasta
````

### Run DNA alignment
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/58_X12_Effectors

less SCNgenome.effector.gff3|awk '$3=="gene"' |wc
    121    1089   17191

No paths found for lcl|flhggfha4D09|Pioneer
No paths found for lcl|flhggfha45D07|ChorismateMutase

#make the table
less SCNgenome.effector.gff3|awk '$3=="gene"' |cut -f 1,4,5,9 |sed 's/ID=//g' |sed 's/\t/:/1' |sed 's/\t/-/1'|less

#grab tn10's stats too
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/01_Transfer2Box
less SCNgenome.effector_sorted.gff|awk '$3=="gene"' |cut -f 1,4,5,9 |sed 's/ID=//g' |sed 's/\t/:/1' |sed 's/\t/-/1'|less
```


### Run Signalp 5.1 on these
```
 ln -s ../01_CondoPorts/28_X12GenomeComparison/pasa2.longest.filter.pep X12proteins.fa
 ~/common_scripts/fasta-splitter.pl --n-parts 8 X12proteins.fa

 for f in *part*fa;do echo "signalp -fasta "$f" -gff3"; done >signalp5.sh

 ml bioawk
bioawk -c fastx '{print $name,length($seq)}' X12proteins.fa >ProteinLengths.txt

ml samtools ; samtools faidx X12proteins.fa


cat *signalp5 |grep -v "#" |sort -k1,1nr |paste - <(sort -k1,1nr ProteinLengths.txt)  |awk '$3>.4999999' |awk '{print $1,substr($7,1,5),$12}' |sed 's/-/ /g' |awk '{print "samtools faidx X12proteins.fa  "$1":"$3"-"$4 " >>SignalPeptidesSubtracted.fasta"}' >ExtractSignalPContainingProts.sh
sh ExtractSignalPContainingProts.sh


#How many are secreted in X12?
grep -c ">" SignalPeptidesSubtracted.fasta
1306


 tmhmm SignalPeptidesSubtracted.fasta >SignalPeptidesSubtractedtmhmm.out

#remove genes with tm domains and get position
less *tab |awk '$2=="0" {print $1}'|

 ln -s ../01_CondoPorts/28_X12GenomeComparison/pasa2.longest.filter.gff3 X12annotation.gff3

  less X12annotation.gff3 |awk '$3=="mRNA"' |sed 's/;/\t/g' |sed 's/ID=//g' >X12annotationGrepMOD

ml bedtools2;bedtools intersect -wo -a X12EffectorPos.bed -b Secreted.bed |less  
```
