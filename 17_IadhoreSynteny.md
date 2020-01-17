# Run iadhore, as mcscanx does not seem reliable

### X12 synteny with TN10 pseudomolecule
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/01_X12/02_iadhoreCircos

#70% identity, 50% length

awk '{print $2}' ../pasa2mikado.blastout |cdbyank ../mikado_proteinsFixed.fasta.cidx |bioawk -c fastx '{print $name,length($seq)}' |paste ../pasa2mikado.blastout - |awk '($4/$14)>.5 && $3>70' |awk '{print $1,$2}' |sort -u -k1,1V >PairwiseOrthology.list



#Could not get synteny with the X12 genes, mapping our genes to x12 genome

#/work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/01_X12/03_mapOurGenes
ln -s ../../../25_AnnotateGenes/mikado_transcripts.fasta
ln -s ../../../28_X12GenomeComparison/SCN_genome.fa
sh runGmap.sh SCN_genome /work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/02_738Assembly/03_mapOurGenes/ SCN_genome.fa mikado_transcripts.fasta

#runGmap.sh
###################################################################
#!/bin/bash

#Makes a database and searches your sequences.
#sh runGmap.sh <database name> <folder of database file ending with a "/"> <Fasta file> <query file>

#examples
#sh run_gmap.sh red_abalone_02Jun2017_5fUJu /work/GIF/remkv6/Serb/03_DavideGMAP/ red_abalone_02Jun2017_5fUJu.fasta DavideQuerydna.fasta
#sh run_gmap.sh  m.yessoensisGenome /work/GIF/remkv6/Serb/03_DavideGMAP DavideQuerydna.fasta
#sh run_gmap.sh Crassostreagigasgenome /work/GIF/remkv6/Serb/03_DavideGMAP Crassostreagigasgenome.fa DavideQuerydna.fasta


module load gmap-gsnap/2018-03-25-qa3kh3t
dbname=$1
dbloc=$2
dbfasta=$3
query=$4
gmap_build -d $dbname  -D $dbloc $dbfasta
gmap -D $dbloc -d $dbname -B 5 -t 16  --input-buffer-size=1000000 --output-buffer-size=1000000 -f gff3_gene  $query >${dbname%.*}.${query%.*}.gff3
###################################################################

#get the ortholog list based on mapping.  2400 genes did not map to the X12 assembly
less GMAP_0.e772363 | grep "No paths" |awk '{print $5}' |grep -v -f - <(grep ">" mikado_transcripts.fasta ) |awk '{print $1}' |sed 's/>//g' |awk '{print $1,$1}' |tr " " "\t" > OrthologousGenes.list


# create iadhore files
/work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/01_X12/02_iadhoreCircos/subject
 less ../../03_mapOurGenes/Path1GenesOnlyMikado2X12.gff |sed 's/ID=//g' |sed 's/;/\t/g' |sed 's/\.path1//g' |awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'

 paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini

```

### make circos plot
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/01_X12/02_iadhoreCircos/01_Circos
 awk '$3=="gene" ' ../../03_mapOurGenes/Path1GenesOnlyMikado2X12.gff |sed 's/ID=//g' |sed 's/;/\t/g' |sed 's/\.path1//g' >X12grepmod.gff
 awk '$3=="gene" ' ../../Gene.mikado.loci.gff3 |sed 's/ID=//g' |sed 's/;/\t/g'  >TN10grepmod.gff

ln -s ../output/segments.txt

 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="TN10genome") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line TN10grepmod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="TN10genome") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line TN10grepmod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="TN10genome") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -m 1 -w $line X12grepmod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="TN10genome") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -m 1 -w $line X12grepmod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8

 #This code uses the list from the first step to paste the scaffold names and scaffold positions in the proper orientation as in the example output below (space separated)
 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '{if($5=="TN10genome") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' |awk 'NR>1' |paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf

ln -s ../../SCNgenome.fasta
ln -s ../../../../28_X12GenomeComparison/SCN_genome.fa X12genome.fa


  bioawk -c fastx '{print $name,length($seq)}' SCNgenome.fasta |awk '{print "chr","-",$1,$1,"0",$2,"blue"}' >TN10Karyotype.conf
  bioawk -c fastx '{print $name,length($seq)}'  X12genome.fa |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >X12Karyotype.conf

awk '{print $1}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' TN10Karyotype.conf >>tmpKaryotype.conf1";done >TN10Karyotype.sh
sh TN10Karyotype.sh
awk '{print $4}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' X12Karyotype.conf >>tmpKaryotype.conf2";done >X12Karyotype.sh
sh X12Karyotype.sh


cat <(sort tmpKaryotype.conf1 |uniq) <(sort tmpKaryotype.conf2 |uniq) >karyotype.conf


circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - Scaffold_1,Scaffold_2,Scaffold_3,Scaffold_4,Scaffold_5,Scaffold_6,Scaffold_7,Scaffold_8,Scaffold_9 -static_rx Scaffold_1,Scaffold_2,Scaffold_3,Scaffold_4,Scaffold_5,Scaffold_6,Scaffold_7,Scaffold_8,Scaffold_9
calculating round 0
report round 0 minimize init 251751 final 165894 change 34.10%
calculating round 1
report round 1 minimize init 165894 final 59291 change 64.26%
calculating round 2
report round 2 minimize init 59291 final 58615 change 1.14%
calculating round 3
report round 3 minimize init 58615 final 58608 change 0.01%
scorereport init 251751 final 58608 change 76.72%
chromosomes_order = chr1,Scaffold_4,Scaffold_8,Scaffold_6,Scaffold_3,Scaffold_5,Scaffold_2,Scaffold_9,chr8,Scaffold_1,chr2,chr_26_pilon,chr9,chr7,chr_21_pilon,chr4,chr_18_pilon,chr_540_pilon,chr_274_pilon,chr3,chr5,chr_288_pilon,chr_292_pilon,chr_303_pilon,chr_725_pilon,Scaffold_7,chr_644_pilon,chr_576_pilon,chr_566_pilon,chr_407_pilon,chr_236_pilon,chr_664_pilon,chr_638_pilon,chr6,chr_528_pilon,chr_651_pilon,chr_443_pilon,chr_130_pilon


less RevisedSyntenicRibbons.conf |awk '{print $3-$2}' |summary.sh
Total:  31,682,741
Count:  500
Mean:   63,365
Median: 40,171
Min:    987
Max:    502,222
[remkv6@condo075 01_Circos]$ less SyntenicRibbons.conf |awk '{print $3-$2}' |summary.sh
Total:  63,601,468
Count:  898
Mean:   70,825
Median: 47,429
Min:    60
Max:    509,173

```
