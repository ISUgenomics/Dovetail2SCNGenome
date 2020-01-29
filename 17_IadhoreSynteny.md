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

#some false alignments here,i.e. super small to super large regions
[remkv6@condo075 01_Circos]$ less SyntenicRibbons.conf |awk '{print $3-$2}' |summary.sh
Total:  63,601,468
Count:  898
Mean:   70,825
Median: 47,429
Min:    60
Max:    509,173

```

### iadhore synteny with 738 genome
```
ln -s ../fixed.augustus.gff3
ln -s ../../01_X12/mikado.loci.gff3



#Primary protein to primary protein isoform blastp

awk '{print $2}' 7382mikado.blastout |cdbyank mikado_proteinsFixed.fasta.cidx |bioawk -c fastx '{print $name,length($seq)}' |paste 7382mikado.blastout - |awk '($4/$14)>.5 && $3>70' |awk '{print $1,$2,"Family"substr($1,9,length($1))}' |sed 's/\.t/\t/2' |awk '{print $1,$3"\n"$2,$3}' |tr " " "\t" >Orthologues.list

less Orthologues.list |sed 's/\.t/\t/g' |sed 's/\./\t/2' |awk '{print $1"\t"$3}' >FixedOrthologues.list


grep ">" ../augustus.aa |sed 's/>//g' |sed 's/\./\t/2' |awk '{print $1}' |sort|uniq>738Genes
grep ">" mikado_proteinsFixed.fasta |sed 's/>//g' |sed 's/\./\t/2' |awk '{print $1}' |sort|uniq>MikadoGenes
cat GenesInFixedOrthologues.list 738Genes|sort|uniq -c |awk '$1==1 {print $2}'|awk '{print $1"\tFamily4727"NR}' >Missing738Genes.list
cat GenesInFixedOrthologues.list MikadoGenes |sort|uniq -c |awk '$1==1 {print $2}'  |awk '{print $1"\tFamily4726"NR}' >MissingPseudomoleculeGenes.list
cat FixedOrthologues.list MissingPseudomoleculeGenes.list Missing738Genes.list >FixedOrthologues.list2


#create non-family based orthologues, to match the synteny calls for x12
awk '{print $2}' 7382mikado.blastout |cdbyank ../mikado_proteinsFixed.fasta.cidx |bioawk -c fastx '{print $name,length($seq)}' |paste 7382mikado.blastout - |awk '($4/$14)>.5 && $3>70' |awk '{print $1,$2}' |sort -u -k1,1V >PairwiseOrthology.list

mkdir query
cd query
awk '$3=="gene"' ../mikado.loci.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'

sed -i 's/ .*//g' *.lst
ls *lst >input.txt
#This can vary also if you have periods "." in your gene names.
 paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini

cd ../
mkdir subject
cd subject/
awk '$3=="gene"' ../fixed.augustus.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close
($2)}'
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini

cat query/query.ini subject/subject.ini |tr "\t" " " >iadhore.ini
vi iadhore.ini

ml GIF2/iAdHoRe/3.0.01
i-adhore iadhore.ini
```

### 738 to pseudo circos plot
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/02_738Assembly/02_iadhore/01_circos

#circos version 0.69.2

mkdir 01_circos
cd 01_circos/
ln -s ../mikado.loci.gff3
ln -s ../fixed.augustus.gff3
ln -s ../../genome738sl.polished.mitoFixed.fa
ln -s ../../../01_X12/SCNgenome.fasta
ln -s ../output/segments.txt

awk '$3=="gene"' mikado.loci.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' >mikadoGrepMod
awk '$3=="gene"' fixed.augustus.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' >738GrepMod

less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line mikadoGrepMod; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line mikadoGrepMod; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -m 1 -w $line 738GrepMod; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -m 1 -w $line 738GrepMod; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '$5!=$1' |awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' | paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf





bioawk -c fastx '{print $name,length($seq)}' SCNgenome.fasta |awk '{print "chr","-",$1,$1,"0",$2,"blue"}' >TN10Karyotype.conf
bioawk -c fastx '{print $name,length($seq)}'  genome738sl.polished.mitoFixed.fa |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >
738Karyotype.conf

awk '{print $1}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' TN10Karyotype.conf >>tmpKaryotype.conf1";done >TN10Karyotype.sh
sh TN10Karyotype.sh
awk '{print $4}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' 738Karyotype.conf >>tmpKaryotype.conf2";done >X12Karyotype.sh
 sh X12Karyotype.sh
cat <(sort tmpKaryotype.conf1 |uniq) <(sort tmpKaryotype.conf2 |uniq) >karyotype.conf
circos -conf circos.conf


vi circos.conf
vi ticks.conf
vi bands.conf
vi ideogram.conf
cp /work/GIF/software/programs/circos/0.69-4/etc/housekeeping.conf .
vi housekeeping.conf

sort tmpKaryotype.conf1 |uniq|awk '{print $3}' |tr "\n" "," |sed 's/.$//' |awk '{print "circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - "$0" -static_rx "$0 }' |less

wget http://circos.ca/distribution/circos-tools-0.22.tgz
tar -zxvf circos-tools-0.22.tgz
circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - Scaffold_1,Scaffold_2,Sca
ffold_3,Scaffold_4,Scaffold_5,Scaffold_6,Scaffold_7,Scaffold_8,Scaffold_9 -static_rx Scaffold_1,Scaffold_2,Scaffold_3,Scaffold_4,Scaffold_5,Scaffold_6,Scaffo
ld_7,Scaffold_8,Scaffold_9

ml GIF2/circos
circos -conf circos.conf



less SyntenicRibbons.conf |awk '{print $3-$2}' |summary.sh
Total:  12,465,597
Count:  247
Mean:   50,468
Median: 36,444
Min:    4,386
Max:    289,224


```
### Make circos histograms representing coverage
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/02_738Assembly/02_iadhore/01_circos

ln -s ../../../../10_RepeatModeler/SCNgenome.fasta.out.gff

less TN10Karyotype.conf |awk '{print $3,$5,$6}' |tr " " "\t" >chrom.sizes

bedtools annotate -both -i SCNwindows.bed -files SCNgenome.fasta.out.gff |awk '{print $1,$2,$3,$5}' >Repeats.histogram

bedtools makewindows -w 10000 -b chrom.sizes >SCNwindows.bed
 bedtools annotate -both -i SCNwindows.bed -files <(awk '$3=="exon"'  mikado.loci.gff3) |awk '{print $1,$2,$3,$5}' >Genes.histogram

less -S 01_TRF/trf.out |awk -v x=0 '{if (substr($1,1,1)=="@") {x=$1} else {print x,$1,$2}}' |sort|uniq|tr " " "\t" |sed 's/@//g'   >TRF.bed
bedtools annotate -both -i SCNwindows.bed -files TRF.bed |awk '{print $1,$2,$3,$5}' >TRF.histogram

bedtools annotate -both -i SCNwindows.bed -files RiboCoords.bed |awk '{print $1,$2,$3,$5}' >Ribosomal.histogram

 bedtools annotate -both -i SCNwindows.bed -files <(awk '$3=="exon"' SCNgenome.effector.gff3) |awk '{print $1,$2,$3,$5}' >Effector.histogram


 less ../../../../29_Effectors/01_Diamond/diamond.out |awk '{print $2}' |grep -f - mikado.loci.gff3 |awk '$3=="exon"' |tr " " "\t">diamondEffector.bed
bedtools annotate -both -i SCNwindows.bed -files diamondEffector.bed |awk '{print $1,$2,$3,$5}' >DiamondEffector.histogram

#circos.conf
###################################################################
karyotype = ./karyotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=SyntenicRibbons.conf
    radius = 0.69r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>
<plots>
 <plot>
   type = histogram
   fill_color =black
   file =Repeats.histogram
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
   file =TRF.histogram
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
   file =Genes.histogram
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
   file =Effector.histogram
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

 <plot>
   type = histogram
   fill_color =black
   file =DiamondEffector.histogram
   r1 = 0.74r
   r0 = 0.69r
   orientation = out
   extend_bin  = no
   <backgrounds>
     <background>
     spacing   = 0.05r
   color = orange
    thickness = 3u
    </background>
  </backgrounds>
</plot>


</plots>
<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 chromosomes_order = scaffold_323,scaffold_106,scaffold_384,scaffold_228,scaffold_218,scaffold_177,scaffold_86,scaffold_36,scaffold_33,scaffold_135,scaffold_425,scaffold_401,scaffold_289,scaffold_68,scaffold_60,scaffold_190,scaffold_122,scaffold_158,scaffold_80,scaffold_165,scaffold_585,scaffold_49,scaffold_2,Scaffold_8,scaffold_64,Scaffold_5,scaffold_326,scaffold_57,scaffold_295,scaffold_223,scaffold_170,scaffold_446,scaffold_39,scaffold_105,scaffold_15,scaffold_84,scaffold_248,scaffold_456,scaffold_334,scaffold_11,scaffold_129,scaffold_227,scaffold_308,scaffold_498,scaffold_47,scaffold_118,scaffold_432,scaffold_372,scaffold_408,Scaffold_3,scaffold_238,scaffold_114,scaffold_307,scaffold_46,scaffold_294,scaffold_40,scaffold_482,scaffold_302,scaffold_471,scaffold_255,scaffold_184,scaffold_50,scaffold_92,scaffold_193,scaffold_101,scaffold_216,scaffold_333,scaffold_23,scaffold_93,Scaffold_2,scaffold_41,scaffold_303,scaffold_5,scaffold_127,scaffold_229,scaffold_438,scaffold_592,scaffold_70,scaffold_462,scaffold_465,scaffold_222,Scaffold_6,scaffold_58,scaffold_142,scaffold_496,scaffold_197,scaffold_83,scaffold_266,scaffold_1,scaffold_51,scaffold_451,scaffold_148,scaffold_424,scaffold_22,scaffold_527,scaffold_418,scaffold_87,Scaffold_7,scaffold_53,scaffold_261,scaffold_287,scaffold_146,scaffold_377,scaffold_155,scaffold_19,scaffold_281,scaffold_285,scaffold_380,scaffold_251,scaffold_111,scaffold_89,scaffold_63,scaffold_95,scaffold_231,scaffold_359,Scaffold_1,scaffold_257,scaffold_79,scaffold_192,scaffold_291,scaffold_282,scaffold_62,scaffold_99,scaffold_30,scaffold_442,scaffold_413,scaffold_356,scaffold_171,Scaffold_4,scaffold_504,scaffold_125,scaffold_45,scaffold_239,scaffold_67,scaffold_16,scaffold_59,scaffold_109,scaffold_12,scaffold_150,scaffold_117,scaffold_73,scaffold_684,Scaffold_9,scaffold_563,scaffold_615,scaffold_427,scaffold_521,scaffold_180,scaffold_113,scaffold_378,scaffold_360,scaffold_123,scaffold_21,scaffold_56,scaffold_444,scaffold_43,scaffold_416

###################################################################
```

### Redo 738 to pseudo synteny using mikado genes mapped to 738

```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/02_738Assembly/02_iadhore/gmapMikadoGenes

ln -s ../../../01_X12/03_mapOurGenes/mikado_transcripts.fasta
ln -s ../01_circos/genome738sl.polished.mitoFixed.fa

echo "sh runGmap.sh genome738sl.polished.mitoFixed /work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/02_738Assembly/02_iadhore/gmapMikadoGenes/ genome738sl.polished.mitoFixed.fa mikado_transcripts.fasta" >gmap.sh

less genome738sl.polished.mikado_transcripts.gff3 |awk '$3=="gene"' |grep "path1" |cut -f 9 |awk '{print $1,$1}' |sed 's/\.path/\t/2' |sed 's/;/\t/1' |sed 's/ID=//g' |tr " " "\t" |cut -f 1,3 >Orthologues.list

cd ..
# Make the query files for iadhore

cd subject/

awk '$3=="gene"' ../gmapMikadoGenes/genome738sl.polished.mikado_transcripts.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print $9$7,$1}' |grep "\.1\.path1" |awk '{print >> $2 ".lst"; close ($2)}'
 sed -i 's/ .*//g' *.lst
 ls *lst >input.txt
 paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
 #removed a few scaffolds from the iadhore.ini that did not have genes mapped to them

 cat gmapMikadoGenes/Orthologues.list |grep "\.1\.path1" >GmapOrthologues.list
 sed 's/\.1/\t/2'  GmapOrthologues.list |awk '{print $2,$1}' |tr " " "\t"  >FixedGmapOrthologues.list

#modified the iadhore.ini file parameters to match the x12 comparison
blast_table=FixedGmapOrthologues.list
prob_cutoff=0.001
anchor_points=3
number_of_threads=16
visualizeAlignment=false
output_path= output
alignment_method=gg2
gap_size=5
cluster_gap=15
level_2_only=true
q_value=.01

 i-adhore iadhore.ini


#Get syntenic ribbons
 awk '$3=="gene"' ../gmapMikadoGenes/genome738sl.polished.mikado_transcripts.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |grep "\.1\.path1" >738GrepMod

less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line mikadoGrepMod; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line mikadoGrepMod; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -m 1 -w $line 738GrepMod; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '$5!=$1'|awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -m 1 -w $line 738GrepMod; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8
less segments.txt |awk 'NR>1' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '$5!=$1' |awk '{if($5=="Pseudomolecule") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' | paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf

less SyntenicRibbons.conf |awk '{print $3-$2}' |summary.sh
Total:  67,013,343
Count:  1,547
Mean:   43,318
Median: 26,042
Min:    11
Max:    533,651

```
