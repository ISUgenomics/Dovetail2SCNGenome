# Run iadhore, as mcscanx does not seem reliable

### X12 putative orthologs with with TN10 pseudomolecule genome annotations
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/01_X12/02_iadhoreCircos
#Now /work/gif/remkv6/Baum/04_DovetailSCNGenome/01_CondoPorts/31_Synteny/01_X12

ln -s ../../../49_RenameChromosomes/mikado.loci.ancestralVHEJ_transcripts.fasta
ln -s ../../../49_RenameChromosomes/mikado.loci.ancestralVHEJ_proteins.fasta

ln -s ../../../49_RenameChromosomes/01_Transfer2Box/OrderedSCNGenePredictions.gff3
ln -s ../../../49_RenameChromosomes/01_Transfer2Box/SCNgenome.fasta TN10genome.fa
ml blast-plus/2.11.0-py3-4pqzweg
makeblastdb -in mikado.loci.ancestralVHEJ_proteins.fasta -dbtype prot -out mikado.loci.ancestralVHEJ_proteins
blastp -db mikado.loci.ancestralVHEJ_proteins -query pasa2.longest.filter.pep -outfmt 6 -num_threads 35 -out  X12pep2TN10pep.blastout



#70% identity, 50% length
ml cdbfasta; ml bioawk
cdbfasta mikado.loci.ancestralVHEJ_proteins.fasta

awk '$3>70 {print $2} ' X12pep2TN10pep.blastout |while read line; do cdbyank mikado.loci.ancestralVHEJ_proteins.fasta.cidx -a $line; done |bioawk -c fastx '{print $name,length($seq)}' |paste <(awk '$3>70' X12pep2TN10pep.blastout) - |awk 'substr($1,1,1)!="#" && ($4/$14)>.5'|awk '{print $1,$2}' |sort -u -k1,1V  >PairwiseOrthology.list &
```


#### Synteny for X12 and TN10 pseudomolecule genomes with tn10 gmapped genes
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/01_X12/03_mapOurGenes
ln -s ../../../../49_RenameChromosomes/mikado.loci.ancestralVHEJ_transcripts.fasta
ln -s ../../../28_X12GenomeComparison/SCN_genome.fa


sh runGmap.sh SCN_genome /work/gif/remkv6/Baum/04_DovetailSCNGenome/01_CondoPorts/31_Synteny/01_X12/03_mapOurGenes/ SCN_genome.fa mikado.loci.ancestralVHEJ_transcripts.fasta

#runGmap.sh
###################################################################
#!/bin/bash

#Makes a database and searches your sequences.
#sh runGmap.sh <database name> <folder of database file ending with a "/"> <Fasta file> <query file>

#examples
#sh run_gmap.sh red_abalone_02Jun2017_5fUJu /work/GIF/remkv6/Serb/03_DavideGMAP/ red_abalone_02Jun2017_5fUJu.fasta DavideQuerydna.fasta
#sh run_gmap.sh  m.yessoensisGenome /work/GIF/remkv6/Serb/03_DavideGMAP DavideQuerydna.fasta
#sh run_gmap.sh Crassostreagigasgenome /work/GIF/remkv6/Serb/03_DavideGMAP Crassostreagigasgenome.fa DavideQuerydna.fasta


module load gmap-gsnap/2019-05-12-zjqshxf
dbname=$1
dbloc=$2
dbfasta=$3
query=$4
gmap_build -d $dbname  -D $dbloc $dbfasta
gmap -D $dbloc -d $dbname -B 5 -t 16  --input-buffer-size=1000000 --output-buffer-size=1000000 -f gff3_gene  $query >${dbname%.*}.${query%.*}.gff3
###################################################################

#How many genes did not map? 22,173 mapped.
less GMAP_0.e3149854|grep "path" |wc
   1760    8800   59840

#get the ortholog list based on mapping.
less GMAP_0.e3149854 | grep "No paths" |awk '{print $5}' |grep -v -f - <(grep ">" mikado.loci.ancestralVHEJ_transcripts.fasta ) |awk '{print $1}' |sed 's/>//g' |awk '{print $1,$1}' |tr " " "\t" > OrthologousGenes.list


# create iadhore files
/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_CondoPorts/31_Synteny/01_X12/02_iadhoreCircos/subject
#gene because this was gmapped with mRNAs which would show as genes.
 awk '$3=="gene"'  ../../03_mapOurGenes/SCN_genome.mikado.loci.ancestralVHEJ_transcripts.gff3 |grep "\.path1;" |sed 's/ID=//g' |sed 's/;/\t/g' |grep "\.t1"|sed 's/\.path1//g' |awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'
 sed -i 's/ .*//g' *.lst
 sed -i 's/Hetgly/XHetgly/g' *lst
 ls *lst >input.txt
 paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini

#/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_CondoPorts/31_Synteny/01_X12/02_iadhoreCircos/query
awk '$3=="mRNA"' ../../OrderedSCNGenePredictions.gff3 |sed 's/ID=//g' |grep "\.t1" |sed 's/;/\t/g' |awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini

cat query/query.ini subject/subject.ini |tr "\t" " " >iadhore.ini

#get the names different so iadhore isnt confused add X to second column.
awk '{print $1"\tX"$2}' NewOrthologousGenes.list |grep "\.t1" >ModNewOrthologousGenes.list

blast_table=Orthologues.list
prob_cutoff=0.05
anchor_points=3
number_of_threads=16
visualizeAlignment=false
output_path= output
alignment_method=gg2
gap_size=5
cluster_gap=15
level_2_only=true
q_value=.05



```

### make circos plot X12 to TN10 pseudo
```
##/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_CondoPorts/31_Synteny/01_X12/02_iadhoreCircos/01_Circos

awk '$3=="gene"' ../../03_mapOurGenes/SCN_genome.mikado.loci.ancestralVHEJ_transcripts.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |grep "path1" |sed 's/\.path1//g' |awk '{print $1,$2,$3,$4,$5,$6,$7,$8,"X"$9}' |tr " " "\t" >X12grepmod.gff

 awk '$3=="mRNA"' ../../OrderedSCNGenePredictions.gff3 |sed 's/ID=//g' |sed 's/;/\t/g'  |grep "\.t1" >TN10grepmod.gff


ln -s ../output/segments.txt

 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="TN10genome") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line TN10grepmod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="TN10genome") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line TN10grepmod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="TN10genome") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -m 1 -w $line X12grepmod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="TN10genome") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -m 1 -w $line X12grepmod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8

 #This code uses the list from the first step to paste the scaffold names and scaffold positions in the proper orientation as in the example output below (space separated)
 less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '{if($5=="TN10genome") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' |awk 'NR>1' |paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf



#Create the karyotypes
ln -s ../../TN10genome.fa
ln -s ../../../../28_X12GenomeComparison/SCN_genome.fa X12genome.fa


bioawk -c fastx '{print $name,length($seq)}' TN10genome.fa |awk '{print "chr","-",$1,$1,"0",$2,"blue"}' >TN10Karyotype.conf
bioawk -c fastx '{print $name,length($seq)}'  X12genome.fa |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >X12Karyotype.conf

awk '{print $1}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' TN10Karyotype.conf >>tmpKaryotype.conf1";done >TN10Karyotype.sh
sh TN10Karyotype.sh
awk '{print $4}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' X12Karyotype.conf >>tmpKaryotype.conf2";done >X12Karyotype.sh
sh X12Karyotype.sh

cat <(sort tmpKaryotype.conf1 |uniq) <(sort tmpKaryotype.conf2 |uniq) >karyotype.conf


#Order the synteny bands
ml circos/0.69-6-fejhzgj
circos-tools-0.22/tools/orderchr/bin/orderchr -links RevisedSyntenicRibbons.conf -karyotype karyotype.conf  -static  Chromosome_1,Chromosome_2,Chromosome_3,Chromosome_4,Chromosome_5,Chromosome_6,Chromosome_7,Chromosome_8,Chromosome_9 -static_rx Chromosome_1,Chromosome_2,Chromosome_3,Chromosome_4,Chromosome_5,Chromosome_6,Chromosome_7,Chromosome_8,Chromosome_9 -init_order  Chromosome_1,Chromosome_2,Chromosome_3,Chromosome_4,Chromosome_5,Chromosome_6,Chromosome_7,Chromosome_8,Chromosome_9


chromosomes_order = Chromosome_9,Chromosome_8,Chromosome_7,Chromosome_2,Chromosome_3,Chromosome_4,Chromosome_6,Chromosome_1,Chromosome_5,chr_508_pilon,chr_303_pilon,chr_725_pilon,chr6,chr_528_pilon,chr_651_pilon,chr_443_pilon,chr_402_pilon,chr_21_pilon,chr4,chr_18_pilon,chr5,chr1,chr_274_pilon,chr3,chr2,chr_287_pilon,chr_26_pilon,chr9,chr_644_pilon,chr_576_pilon,chr_592_pilon,chr_407_pilon,chr_566_pilon,chr_236_pilon,chr_638_pilon,chr7,chr8


#some false alignments here,i.e. super small to super large regions a

less SyntenicRibbons.conf |sort|uniq|awk '{print $3-$2}' |summary.sh
Total:  93,442,625
Count:  738
Mean:   126,616
Median: 67,936
Min:    680
Max:    1,276,810




#remove synteny where it is 5x greater in length in one genome vs the other -- too unlikely
less SyntenicRibbons.conf  |awk '{if($6>$5) {print $3-$2, $6-$5,$0}else {print $3-$2, $5-$6,$0}}'|awk '{if($1> (5*$2) || $2> ($1*5)) {next; } else {print $3,$4,$5,$6,$7,$8}}' |sort|uniq   >RevisedSyntenicRibbons.conf


less RevisedSyntenicRibbons.conf |sort|uniq|awk '{print $3-$2}' |summary.sh
Total:  80,613,520
Count:  608
Mean:   132,588
Median: 74,145
Min:    2,387
Max:    1,276,810



#rename chromosomes to tn10 and x12

sed -i 's/chr/X12_/2'  karyotype.conf |
sed -i 's/Chromosome_/TN10_/g' karyotype.conf
sed -i 's/_pilon//g' karyotype.conf
 sed -i 's/__/_/g' karyotype.conf

sed -i 's/Chromosome_/TN10_/g' RevisedSyntenicRibbons.conf
sed  -i 's/chr/X12_/g' RevisedSyntenicRibbons.conf
sed -i 's/_pilon//g' RevisedSyntenicRibbons.conf
sed -i 's/__/_/g' RevisedSyntenicRibbons.conf

sed -i 's/__/_/g' circos.conf
sed -i 's/_pilon//g' circos.conf
sed -i 's/,chr/,X12_/g' circos.conf
sed -i 's/Chromosome_/TN10_/g' circos.conf


circos -conf circos.conf
```




### Redo 738 to pseudo synteny using mikado genes mapped to 738
##### map pseudo genes to 738
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_CondoPorts/31_Synteny/02_738Assembly


wget https://scnbase.org/files/download/genome738sl.polished.mitofixed.fasta
ln -s ../../../49_RenameChromosomes/mikado.loci.ancestralVHEJ_transcripts.fasta
ln -s ../../../49_RenameChromosomes/01_Transfer2Box/OrderedSCNGenePredictions.gff3
ln -s ../../../49_RenameChromosomes/01_Transfer2Box/SCNgenome.fasta TN10genome.fa



echo "sh runGmap.sh genome738sl.polished.mitofixed /work/gif/remkv6/Baum/04_DovetailSCNGenome/01_CondoPorts/31_Synteny/02_738Assembly/  genome738sl.polished.mitofixed.fa mikado.loci.ancestralVHEJ_transcripts.fasta" >gmap.sh

#How many did not map?
grep "No paths" gmap_0.e3150021 |wc
    712    3560   24208

******


less genome738sl.polished.mikado.loci.ancestralVHEJ_transcripts.gff3 |awk '$3=="gene"' |grep "path1" |cut -f 9 |awk '{print $1,"X"$1}' |sed 's/\.path/\t/2' |sed 's/;/\t/1' |sed 's/ID=//g' |tr " " "\t" |cut -f 1,3 |sed 's/\.path1//g' |grep "\.t1" >Orthologues.list
```

##### run iadhore on 738 and pseudo tn10
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_CondoPorts/31_Synteny/02_738Assembly/02_iadhore

# Make the query files for iadhore
cd query
cp ../../../01_X12/02_iadhoreCircos/query/* .
cd ../subject/
awk '$3=="gene"' ../../genome738sl.polished.mikado.loci.ancestralVHEJ_transcripts.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print $9$7,$1}' |grep "path1" |grep "\.t1" |sed 's/\.path1//g' |awk '{print >> $2 ".lst"; close ($2)}'
 sed -i 's/ .*//g' *.lst
 sed -i 's/^/X/g' *lst
 ls *lst >input.txt
 paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
cd ../; cat query/query.ini subject/subject.ini |tr "\t" " " >iadhore.ini


#modified the iadhore.ini file parameters to match the x12 comparison
blast_table=FixedGmapOrthologues.list
prob_cutoff=0.05
anchor_points=3
number_of_threads=16
visualizeAlignment=false
output_path= output
alignment_method=gg2
gap_size=5
cluster_gap=15
level_2_only=true
q_value=.05

 i-adhore iadhore.ini

#create grepmods
awk '$3=="gene"' ../../genome738sl.polished.mikado.loci.ancestralVHEJ_transcripts.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |grep "path1" |sed 's/\.path1//g' |awk '{print $1,$2,$3,$4,$5,$6,$7,$8,"X"$9}' |tr " " "\t" |grep "\t1" >738grepmod.gff
awk '$3=="mRNA"' ../../OrderedSCNGenePredictions.gff3 |sed 's/ID=//g' |sed 's/;/\t/g'  |grep "\.t1" >TN10grepmod.gff
#Get syntenic ribbons

less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line TN10grepmod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line TN10grepmod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -m 1 -w $line 738grepmod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -m 1 -w $line 738grepmod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '$5!=$1' |awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' | paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf

less SyntenicRibbons.conf |awk '{print $3-$2}' |sort|uniq|summary.sh
Total:  74,028,909
Count:  1,252
Mean:   59,128
Median: 36,023
Min:    183
Max:    871,367

#remove those that are five fold larger, as they are likely artifacts
less SyntenicRibbons.conf  |awk '{if($6>$5) {print $3-$2, $6-$5,$0}else {print $3-$2, $5-$6,$0}}'|awk '{if($1> (5*$2) || $2> ($1*5)) {next; } else {print $3,$4,$5,$6,$7,$8}}' |sort|uniq   >RevisedSyntenicRibbons.conf

less RevisedSyntenicRibbons.conf |awk '{print $3-$2}' |sort|uniq|summary.sh
Total:  65,626,657
Count:  1,130
Mean:   58,076
Median: 33,807
Min:    1,398
Max:    871,367

ln -s ../../../01_X12/TN10genome.fa
ln -s ../../genome738sl.polished.mitofixed.fasta
bioawk -c fastx '{print $name,length($seq)}' TN10genome.fa |awk '{print "chr","-",$1,$1,"0",$2,"blue"}' >TN10Karyotype.conf
bioawk -c fastx '{print $name,length($seq)}' genome738sl.polished.mitofixed.fasta  |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >738Karyotype.conf

awk '{print $1}' RevisedSyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' TN10Karyotype.conf >>tmpKaryotype.conf1";done >TN10Karyotype.sh
sh TN10Karyotype.sh
awk '{print $4}' RevisedSyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' 738Karyotype.conf >>tmpKaryotype.conf2";done >X12Karyotype.sh
sh X12Karyotype.sh

cat <(sort tmpKaryotype.conf1 |uniq) <(sort tmpKaryotype.conf2 |uniq) >karyotype.conf

#Order the synteny bands
ml circos/0.69-6-fejhzgj
circos-tools-0.22/tools/orderchr/bin/orderchr -links RevisedSyntenicRibbons.conf -karyotype karyotype.conf  -static  Chromosome_1,Chromosome_2,Chromosome_3,Chromosome_4,Chromosome_5,Chromosome_6,Chromosome_7,Chromosome_8,Chromosome_9 -static_rx Chromosome_1,Chromosome_2,Chromosome_3,Chromosome_4,Chromosome_5,Chromosome_6,Chromosome_7,Chromosome_8,Chromosome_9 -init_order  Chromosome_1,Chromosome_2,Chromosome_3,Chromosome_4,Chromosome_5,Chromosome_6,Chromosome_7,Chromosome_8,Chromosome_9

#note that this never came out of the above command correctly.  It would rearrange names.  Not sure of the problem, too many perhaps.  cleaned this up manually
chromosomes_order = Chromosome_8,Chromosome_9,Chromosome_6,Chromosome_7,Chromosome_2,Chromosome_4,Chromosome_5,Chromosome_3,Chromosome_1,scaffold_713,scaffold_351,scaffold_230,scaffold_382,scaffold_179,scaffold_505,scaffold_31,scaffold_484,scaffold_332,scaffold_326,scaffold_412,scaffold_223,scaffold_411,scaffold_315,scaffold_57,scaffold_30,scaffold_453,scaffold_423,scaffold_407,scaffold_659,scaffold_164,scaffold_219,scaffold_295,scaffold_383,scaffold_126,scaffold_446,scaffold_322,scaffold_638,scaffold_105,scaffold_343,scaffold_39,scaffold_24,scaffold_101,scaffold_243,scaffold_84,scaffold_170,scaffold_248,scaffold_500,scaffold_66,scaffold_291,scaffold_222,scaffold_607,scaffold_52,scaffold_11,scaffold_76,scaffold_532,scaffold_147,scaffold_129,scaffold_309,scaffold_227,scaffold_267,scaffold_137,scaffold_528,scaffold_452,scaffold_118,scaffold_183,scaffold_121,scaffold_373,scaffold_301,scaffold_139,scaffold_399,scaffold_255,scaffold_292,scaffold_402,scaffold_203,scaffold_398,scaffold_396,scaffold_308,scaffold_241,scaffold_388,scaffold_674,scaffold_174,scaffold_406,scaffold_35,scaffold_311,scaffold_134,scaffold_268,scaffold_169,scaffold_545,scaffold_258,scaffold_729,scaffold_621,scaffold_98,scaffold_346,scaffold_234,scaffold_40,scaffold_91,scaffold_333,scaffold_159,scaffold_114,scaffold_642,scaffold_526,scaffold_130,scaffold_350,scaffold_145,scaffold_184,scaffold_279,scaffold_46,scaffold_9,scaffold_304,scaffold_307,scaffold_29,scaffold_671,scaffold_497,scaffold_294,scaffold_471,scaffold_275,scaffold_588,scaffold_65,scaffold_355,scaffold_609,scaffold_481,scaffold_342,scaffold_477,scaffold_358,scaffold_50,scaffold_72,scaffold_92,scaffold_162,scaffold_193,scaffold_208,scaffold_168,scaffold_393,scaffold_156,scaffold_151,scaffold_216,scaffold_185,scaffold_154,scaffold_459,scaffold_397,scaffold_177,scaffold_289,scaffold_23,scaffold_249,scaffold_345,scaffold_265,scaffold_125,scaffold_87,scaffold_440,scaffold_425,scaffold_176,scaffold_88,scaffold_353,scaffold_38,scaffold_149,scaffold_215,scaffold_401,scaffold_64,scaffold_313,scaffold_288,scaffold_676,scaffold_119,scaffold_400,scaffold_42,scaffold_10,scaffold_205,scaffold_106,scaffold_61,scaffold_504,scaffold_86,scaffold_666,scaffold_20,scaffold_298,scaffold_36,scaffold_664,scaffold_33,scaffold_135,scaffold_122,scaffold_60,scaffold_233,scaffold_344,scaffold_68,scaffold_319,scaffold_662,scaffold_611,scaffold_15,scaffold_206,scaffold_45,scaffold_80,scaffold_585,scaffold_246,scaffold_272,scaffold_82,scaffold_49,scaffold_224,scaffold_7,scaffold_466,scaffold_94,scaffold_262,scaffold_270,scaffold_658,scaffold_190,scaffold_2,scaffold_356,scaffold_316,scaffold_643,scaffold_34,scaffold_132,scaffold_376,scaffold_591,scaffold_181,scaffold_493,scaffold_44,scaffold_652,scaffold_571,scaffold_305,scaffold_67,scaffold_16,scaffold_437,scaffold_59,scaffold_253,scaffold_606,scaffold_109,scaffold_12,scaffold_527,scaffold_240,scaffold_210,scaffold_73,scaffold_97,scaffold_422,scaffold_117,scaffold_150,scaffold_451,scaffold_191,scaffold_328,scaffold_217,scaffold_324,scaffold_323,scaffold_79,scaffold_680,scaffold_214,scaffold_75,scaffold_325,scaffold_431,scaffold_48,scaffold_256,scaffold_14,scaffold_413,scaffold_449,scaffold_17,scaffold_63,scaffold_599,scaffold_336,scaffold_475,scaffold_348,scaffold_257,scaffold_18,scaffold_510,scaffold_96,scaffold_417,scaffold_562,scaffold_394,scaffold_589,scaffold_225,scaffold_192,scaffold_550,scaffold_314,scaffold_300,scaffold_1,scaffold_352,scaffold_442,scaffold_278,scaffold_273,scaffold_327,scaffold_245,scaffold_282,scaffold_62,scaffold_99,scaffold_4,scaffold_276,scaffold_136,scaffold_363,scaffold_594,scaffold_283,scaffold_171,scaffold_182,scaffold_153,scaffold_58,scaffold_188,scaffold_74,scaffold_120,scaffold_28,scaffold_429,scaffold_95,scaffold_143,scaffold_239,scaffold_116,scaffold_104,scaffold_53,scaffold_112,scaffold_242,scaffold_107,scaffold_261,scaffold_103,scaffold_85,scaffold_472,scaffold_287,scaffold_146,scaffold_377,scaffold_25,scaffold_250,scaffold_251,scaffold_155,scaffold_349,scaffold_281,scaffold_144,scaffold_285,scaffold_340,scaffold_8,scaffold_380,scaffold_102,scaffold_111,scaffold_90,scaffold_228,scaffold_54,scaffold_142,scaffold_163,scaffold_108,scaffold_204,scaffold_89,scaffold_264,scaffold_274,scaffold_419,scaffold_209,scaffold_577,scaffold_293,scaffold_613,scaffold_19,scaffold_69,scaffold_220,scaffold_434,scaffold_157,scaffold_496,scaffold_448,scaffold_523,scaffold_543,scaffold_197,scaffold_201,scaffold_312,scaffold_637,scaffold_244,scaffold_83,scaffold_266,scaffold_247,scaffold_166,scaffold_133,scaffold_51,scaffold_260,scaffold_148,scaffold_468,scaffold_160,scaffold_26,scaffold_13,scaffold_424,scaffold_232,scaffold_408,scaffold_77,scaffold_578,scaffold_221,scaffold_93,scaffold_347,scaffold_207,scaffold_172,scaffold_329,scaffold_341,scaffold_280,scaffold_22,scaffold_100,scaffold_213,scaffold_37,scaffold_180,scaffold_439,scaffold_113,scaffold_438,scaffold_378,scaffold_360,scaffold_123,scaffold_187,scaffold_21,scaffold_56,scaffold_81,scaffold_444,scaffold_284,scaffold_6,scaffold_271,scaffold_173,scaffold_457,scaffold_43,scaffold_110,scaffold_533,scaffold_140,scaffold_506,scaffold_462,scaffold_512,scaffold_416,scaffold_3,scaffold_303,scaffold_41,scaffold_5,scaffold_127,scaffold_701,scaffold_385,scaffold_226,scaffold_78,scaffold_483,scaffold_229,scaffold_71,scaffold_317,scaffold_366,scaffold_70,scaffold_124,scaffold_331,scaffold_339,scaffold_719,scaffold_189,scaffold_252,scaffold_47,scaffold_138,scaffold_445,scaffold_55
```


### Make circos histograms representing coverage --
```
/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_CondoPorts/31_Synteny/12_DisplayHistograms
ln -s ../01_X12/TN10genome.fa

#make histogram windows
ml bioawk ;bioawk -c fastx '{print $name"\t0\t"length($seq)}' TN10genome.fa  >chrom.sizes
bedtools makewindows -w 10000 -b chrom.sizes >SCNwindows.bed

ln -s /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/01_Transfer2Box/EDTAtransposons.gff3
ln -s /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/01_Transfer2Box/OrderedSCNGenePredictions.gff3
gunzip  /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/01_Transfer2Box/RepeatMaskerFormatted.gff3.gz
gunzip  /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/01_Transfer2Box/RiboCoords.bed.gz
gunzip  /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/01_Transfer2Box/SortedTRF.bed.gz

ln -s   /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/01_Transfer2Box/RepeatMaskerFormatted.gff3
ln -s   /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/01_Transfer2Box/RiboCoords.bed
ln -s  /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/01_Transfer2Box/SortedTRF.bed


#create histograms
bedtools annotate -both -i SCNwindows.bed -files EDTAtransposons.gff3 |awk '{print $1,$2,$3,$5}' >EDTA.histogram
bedtools annotate -both -i SCNwindows.bed -files RepeatMaskerFormatted.gff3 |awk '{print $1,$2,$3,$5}' >Repeatmasker.histogram
cat  <(awk '{print $1,$4,$5}' EDTAtransposons.gff3 |tr " " "\t")  <(awk '{print $1,$4,$5}' RepeatMaskerFormatted.gff3|tr " " "\t") |sort -k1,1V -k2,3n|bedtools merge -i - >AllRepeats.bed


bedtools annotate -both -i SCNwindows.bed -files <(awk '$3=="exon"'  OrderedSCNGenePredictions.gff3) |awk '{print $1,$2,$3,$5}' >Exons.histogram
bedtools annotate -both -i SCNwindows.bed -files <(awk '$3=="gene"'  OrderedSCNGenePredictions.gff3) |awk '{print $1,$2,$3,$5}' >Genes.histogram

bedtools annotate -both -i SCNwindows.bed -files SortedTRF.bed |awk '{print $1,$2,$3,$5}' >TRF.histogram

bedtools annotate -both -i SCNwindows.bed -files RiboCoords.bed |awk '{print $1,$2,$3,$5}' >Ribosomal.histogram


#circos.conf
###################################################################
karyotype = ./karyotype.conf
chromosomes_units = 1000000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <plots>
 <plot>
   type = histogram
   fill_color =black
   file =Genes.histogram
   r1 = 0.99r
   r0 = 0.94r
   orientation = out
   extend_bin  = no
   <backgrounds>
     <background>
     spacing   = 0.05r
   color = yellow
    thickness = 3u
    </background>
  </backgrounds>
</plot>
 <plot>
   type = histogram
   fill_color =black
   file =Exons.histogram
   r1 = 0.94r
   r0 = 0.89r
   orientation = out
   extend_bin  = no
   <backgrounds>
     <background>
     spacing   = 0.05r
   color = purple
    thickness = 3u
    </background>
  </backgrounds>
</plot>
 <plot>
   type = histogram
   fill_color =black
   file =Ribosomal.histogram
   r1 = 0.89r
   r0 = 0.84r
   orientation = out
   extend_bin  = no
   <backgrounds>
     <background>
     spacing   = 0.05r
   color = blue
    thickness = 3u
    </background>
  </backgrounds>
</plot>

 <plot>
   type = histogram
   fill_color =black
   file =AllRepeats.histogram
   r1 = 0.84r
   r0 = 0.79r
   orientation = out
   extend_bin  = no
   <backgrounds>
     <background>
     spacing   = 0.05r
   color = red
    thickness = 3u
    </background>
  </backgrounds>
</plot>
 <plot>
   type = histogram
   fill_color =black
   file =TRF.histogram
   r1 = 0.79r
   r0 = 0.74r
   orientation = out
   extend_bin  = no
   <backgrounds>
     <background>
     spacing   = 0.05r
   color = green
    thickness = 3u
    </background>
    </backgrounds>
</plot>
</plots>
<<include /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/circos-0.69-6-fejhzgj2bzy7nkrgzto3hr5kcvajbvyz/lib/circos/etc/housekeeping.conf>>

<image>
 <<include /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/circos-0.69-6-fejhzgj2bzy7nkrgzto3hr5kcvajbvyz/lib/circos/etc/image.conf>>

</image>

<<include /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/circos-0.69-6-fejhzgj2bzy7nkrgzto3hr5kcvajbvyz/lib/circos/etc/colors_fonts_patterns.conf>>

###################################################################
