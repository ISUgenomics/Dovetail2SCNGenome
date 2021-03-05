# Synteny for all related Tylenchida genomes

### G. pallida
```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/01_pallida
#grabbing only primary isoform protein fastas
ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/globodera_pallida.PRJEB123.WBPS10.protein.fa
ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/DovetailSCNMaker4.all.maker.proteins.Isoform1Only.fasta

#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/
module load orthofinder
orthofinder -f 01_pallida

#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/01_pallida/01_iadhore/
less ../Results_Jun21/Orthogroups.txt |tr " " "\n" |awk -v Family=0 '{if(substr($1,1,2)=="OG") {Family=$1} else {print $0,Family}}' |sed 's/://g' |cat - <(dos2unix ../Results_Jun21/Orthogroups_UnassignedGenes.csv|awk 'NR>1 {print $2,$1}') |tr " " "\t" >Orthologues.list

ln -s ../../../12_MakerGenesOrthofinder/globodera_pallida.PRJEB123.WBPS10.annotations.gff3
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff

mkdir query
less ../DovetailSCNMaker4.all.NOFASTA.gff| awk '$3=="mRNA"' |sed 's/;/\t/g' |sed 's/ID=//g' |grep "RA"|awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
  paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini

mkdir subject
#make sure to get only the primary transcript
less ../globodera_pallida.PRJEB123.WBPS10.annotations.gff3 |awk '$3=="gene"' |sed 's/;/\t/g' |sed 's/ID=gene://g' |awk '{print $9$7,$1}' |sed 's/pathogens_Gpal_//g' |awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
cd ..
cat query/query.ini subject/subject.ini |tr "\t" " " >iadhore.ini

blast_table=Orthologues.list
table_type=family
prob_cutoff=0.001
anchor_points=3
number_of_threads=16
visualizeAlignment=true
output_path= output
alignment_method=gg2
gap_size=15
cluster_gap=20
level_2_only=true
q_value=.05
```
Circos setup
```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/01_pallida/02_circos
ln -s ../01_iadhore/output/segments.txt
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.genome.fasta
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff
ln -s ../../../12_MakerGenesOrthofinder/globodera_pallida.PRJEB123.WBPS10.annotations.gff3
ln -s ../../../12_MakerGenesOrthofinder/globodera_pallida.PRJEB123.WBPS10.genomic.fa

sed 's/;/\t/g' DovetailSCNMaker4.all.NOFASTA.gff |sed 's/ID=//g' |awk '$3=="mRNA"' |grep "RA" >SCNGrepMod.gff

 sed 's/;/\t/g' globodera_pallida.PRJEB123.WBPS10.annotations.gff3 |sed 's/ID=gene://g' |awk '$3=="gene"' |sed 's/pathogens_Gpal_//g' >GPALGrepMod.gff

less segments.txt |awk 'NR>1'  |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="G.pallida") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line GPALGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
 less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="G.pallida") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line GPALGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
 less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="G.pallida") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
 less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'|awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="G.pallida") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8

less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '{if($5=="G.pallida") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' |awk 'NR>1' |paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf



bioawk -c fastx '{print $name,length($seq)}' globodera_pallida.PRJEB123.WBPS10.genomic.fa |awk '{print "chr","-",$1,$1,"0",$2,"blue"}' |sed 's/pathogens_Gpal_//g' >G.pallidaKaryotype.conf
#had to modify the G.pallida karyotype file, as some duplicately named scaffold was present (scaffold_214.1 and scaffold_214.2).  I deleted scaffold_214.1 because it was only 15kb.
bioawk -c fastx '{print $name,length($seq)}' DovetailSCNMaker4.genome.fasta |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >H.glycinesKaryotype.conf


awk '{print $1}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' G.pallidaKaryotype.conf >>tmpKaryotype.conf1";done >G.pallidaKaryotype.sh
 sh G.pallidaKaryotype.sh

awk '{print $4}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' H.glycinesKaryotype.conf >>tmpKaryotype.conf2";done >H.glycinesKaryotype.sh
sh H.glycinesKaryotype.sh
cat <(sort tmpKaryotype.conf1 |uniq) <(sort tmpKaryotype.conf2 |uniq) >karyotype.conf

sort tmpKaryotype.conf2 |uniq|awk '{print $3}' |tr "\n" "," |sed 's/.$//' |awk '{print "circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - "$0" -static_rx "$0 }' |less
circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_51,Scaffold_60,Scaffold_70,Scaffold_73,Scaffold_80,Scaffold_90 -static_rx Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_51,Scaffold_60,Scaffold_70,Scaffold_73,Scaffold_80,Scaffold_90
calculating round 0
report round 0 minimize init 22109 final 2873 change 87.01%
calculating round 1
report round 1 minimize init 2873 final 2565 change 10.72%
calculating round 2
report round 2 minimize init 2565 final 2082 change 18.83%
calculating round 3
report round 3 minimize init 2082 final 1906 change 8.45%
scorereport init 22109 final 1906 change 91.38%
chromosomes_order = Scaffold_3,scaffold_85,scaffold_25,scaffold_110,scaffold_102,scaffold_222,scaffold_539,scaffold_72,scaffold_53,scaffold_87,scaffold_325,scaffold_444,scaffold_8,scaffold_43,scaffold_278,scaffold_623,scaffold_447,scaffold_412,scaffold_111,scaffold_163,scaffold_1578,scaffold_2159,scaffold_139,scaffold_41,scaffold_930,scaffold_146,scaffold_74,scaffold_619,scaffold_596,scaffold_91,scaffold_390,scaffold_1503,scaffold_233,scaffold_564,scaffold_483,scaffold_296,Scaffold_73,scaffold_826,scaffold_92,scaffold_34,scaffold_311,scaffold_224,scaffold_526,scaffold_196,scaffold_410,scaffold_30,scaffold_145,scaffold_499,scaffold_240,scaffold_47,scaffold_67,scaffold_15,scaffold_231,scaffold_293,scaffold_125,scaffold_257,scaffold_327,scaffold_27,scaffold_271,scaffold_446,scaffold_561,scaffold_226,scaffold_5,scaffold_42,scaffold_1020,scaffold_23,Scaffold_60,scaffold_147,scaffold_90,scaffold_1181,scaffold_155,scaffold_100,scaffold_6,scaffold_38,scaffold_14,scaffold_13,scaffold_20,scaffold_704,scaffold_11,scaffold_131,scaffold_171,scaffold_180,scaffold_50,scaffold_223,scaffold_103,scaffold_71,scaffold_194,scaffold_616,scaffold_573,Scaffold_25,scaffold_22,scaffold_381,scaffold_588,scaffold_3,scaffold_285,scaffold_785,scaffold_323,scaffold_540,scaffold_369,scaffold_193,scaffold_18,scaffold_228,scaffold_170,scaffold_148,scaffold_1364,scaffold_402,scaffold_1325,scaffold_373,scaffold_48,scaffold_31,scaffold_300,scaffold_114,scaffold_84,scaffold_557,scaffold_414,scaffold_107,scaffold_287,scaffold_136,Scaffold_44,scaffold_2,scaffold_443,scaffold_404,scaffold_109,scaffold_159,scaffold_527,scaffold_377,scaffold_260,scaffold_117,scaffold_112,scaffold_12,scaffold_81,scaffold_79,scaffold_1389,scaffold_40,scaffold_214,scaffold_309,scaffold_1446,scaffold_73,scaffold_77,Scaffold_19,scaffold_19,scaffold_328,scaffold_253,scaffold_683,scaffold_119,scaffold_207,scaffold_910,scaffold_335,scaffold_552,scaffold_168,scaffold_54,scaffold_56,scaffold_123,scaffold_441,scaffold_59,scaffold_153,scaffold_374,scaffold_246,scaffold_618,scaffold_143,scaffold_209,Scaffold_103,scaffold_129,scaffold_217,scaffold_44,scaffold_69,scaffold_382,scaffold_211,scaffold_329,scaffold_1166,scaffold_653,scaffold_399,scaffold_9,scaffold_284,scaffold_45,scaffold_36,scaffold_63,scaffold_545,Scaffold_80,scaffold_178,scaffold_505,scaffold_611,scaffold_160,scaffold_21,scaffold_68,scaffold_529,scaffold_530,scaffold_236,scaffold_1244,Scaffold_104,scaffold_94,scaffold_52,scaffold_225,scaffold_1394,scaffold_210,scaffold_29,scaffold_58,scaffold_386,Scaffold_13,scaffold_1,scaffold_35,Scaffold_70,Scaffold_51,scaffold_126,Scaffold_90




#circos.conf
#############################################################################
karyotype = ./karyotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
chromosomes_order = Scaffold_3,scaffold_85,scaffold_25,scaffold_110,scaffold_102,scaffold_222,scaffold_539,scaffold_72,scaffold_53,scaffold_87,scaffold_325,scaffold_444,scaffold_8,scaffold_43,scaffold_278,scaffold_623,scaffold_447,scaffold_412,scaffold_111,scaffold_163,scaffold_1578,scaffold_2159,scaffold_139,scaffold_41,scaffold_930,scaffold_146,scaffold_74,scaffold_619,scaffold_596,scaffold_91,scaffold_390,scaffold_1503,scaffold_233,scaffold_564,scaffold_483,scaffold_296,Scaffold_73,scaffold_826,scaffold_92,scaffold_34,scaffold_311,scaffold_224,scaffold_526,scaffold_196,scaffold_410,scaffold_30,scaffold_145,scaffold_499,scaffold_240,scaffold_47,scaffold_67,scaffold_15,scaffold_231,scaffold_293,scaffold_125,scaffold_257,scaffold_327,scaffold_27,scaffold_271,scaffold_446,scaffold_561,scaffold_226,scaffold_5,scaffold_42,scaffold_1020,scaffold_23,Scaffold_60,scaffold_147,scaffold_90,scaffold_1181,scaffold_155,scaffold_100,scaffold_6,scaffold_38,scaffold_14,scaffold_13,scaffold_20,scaffold_704,scaffold_11,scaffold_131,scaffold_171,scaffold_180,scaffold_50,scaffold_223,scaffold_103,scaffold_71,scaffold_194,scaffold_616,scaffold_573,Scaffold_25,scaffold_22,scaffold_381,scaffold_588,scaffold_3,scaffold_285,scaffold_785,scaffold_323,scaffold_540,scaffold_369,scaffold_193,scaffold_18,scaffold_228,scaffold_170,scaffold_148,scaffold_1364,scaffold_402,scaffold_1325,scaffold_373,scaffold_48,scaffold_31,scaffold_300,scaffold_114,scaffold_84,scaffold_557,scaffold_414,scaffold_107,scaffold_287,scaffold_136,Scaffold_44,scaffold_2,scaffold_443,scaffold_404,scaffold_109,scaffold_159,scaffold_527,scaffold_377,scaffold_260,scaffold_117,scaffold_112,scaffold_12,scaffold_81,scaffold_79,scaffold_1389,scaffold_40,scaffold_214,scaffold_309,scaffold_1446,scaffold_73,scaffold_77,Scaffold_19,scaffold_19,scaffold_328,scaffold_253,scaffold_683,scaffold_119,scaffold_207,scaffold_910,scaffold_335,scaffold_552,scaffold_168,scaffold_54,scaffold_56,scaffold_123,scaffold_441,scaffold_59,scaffold_153,scaffold_374,scaffold_246,scaffold_618,scaffold_143,scaffold_209,Scaffold_103,scaffold_129,scaffold_217,scaffold_44,scaffold_69,scaffold_382,scaffold_211,scaffold_329,scaffold_1166,scaffold_653,scaffold_399,scaffold_9,scaffold_284,scaffold_45,scaffold_36,scaffold_63,scaffold_545,Scaffold_80,scaffold_178,scaffold_505,scaffold_611,scaffold_160,scaffold_21,scaffold_68,scaffold_529,scaffold_530,scaffold_236,scaffold_1244,Scaffold_104,scaffold_94,scaffold_52,scaffold_225,scaffold_1394,scaffold_210,scaffold_29,scaffold_58,scaffold_386,Scaffold_13,scaffold_1,scaffold_35,Scaffold_70,Scaffold_51,scaffold_126,Scaffold_90

#############################################################################

ticks.conf
###############################################################################
show_ticks          = yes
show_tick_labels    = yes
<ticks>
    radius           = 1r
    color            = black
    thickness        = 10p
    multiplier       = 1e-5
    format           = %d
 <tick>
    spacing        = 5u
    size           = 25p
    show_label     = yes
    label_size     = 25p
    label_offset   = 10p
    format         = %d
  </tick>

</ticks>
###############################################################################

bands.conf
###############################################################################
<bands>
   show_bands = yes
   fill_bands = yes
   band_transparency = 4
</bands>
###############################################################################

cp /work/GIF/software/programs/circos/0.69-4/etc/housekeeping.conf .


 ```

### G. pallida
```
 #/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/02_ellingtonae
 #grabbing only primary isoform protein fastas
 ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/DovetailSCNMaker4.all.maker.proteins.Isoform1Only.fasta
 ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/G.ellingtonae_protein.Isoform1Only.fasta

 #/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/
 module load orthofinder
 orthofinder -f 02_ellingtonae

#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/02_ellingtonae/01_iadhore
less ../Results_Jun21/Orthogroups.txt |tr " " "\n" |awk -v Family=0 '{if(substr($1,1,2)=="OG") {Family=$1} else {print $0,Family}}' |sed 's/://g' |cat - <(dos2unix ../Results_Jun21/Orthogroups_UnassignedGenes.csv|awk 'NR>1 {print $2,$1}') |tr " " "\t" >Orthologues.list
#had to modify the gene names for ellingtonae
sed -i 's/Gell//g' Orthologues.list
less Orthologues.list |sed 's/\.t1//g' >fixed.Orthologues.list

ln -s ../../../12_MakerGenesOrthofinder/G.ellingtonae.gff3
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff
mkdir query
less ../DovetailSCNMaker4.all.NOFASTA.gff| awk '$3=="mRNA"' |sed 's/;/\t/g' |sed 's/ID=//g' |grep "RA"|awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini

 mkdir subject
 #make sure to get only the primary transcript
 less ../G.ellingtonae.gff3 |awk '$3=="gene"' |sed 's/;/\t/g' |sed 's/ID=//g' |awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'
 sed -i 's/ .*//g' *.lst
 ls *lst >input.txt
 paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
 cd ..
 cat query/query.ini subject/subject.ini |tr "\t" " " >iadhore.ini

 blast_table=Orthologues.list
 table_type=family
 prob_cutoff=0.001
 anchor_points=3
 number_of_threads=16
 visualizeAlignment=true
 output_path= output
 alignment_method=gg2
 gap_size=15
 cluster_gap=20
 level_2_only=true
 q_value=.05
 ```
 Circos setup
 ```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/02_ellingtonae/02_circos
ln -s ../01_iadhore/output/segments.txt
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.genome.fasta
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff
ln -s ../../../12_MakerGenesOrthofinder/G.ellingtonae.gff3
ln -s ../../../12_MakerGenesOrthofinder/G.ellingtonae.genomic.fa

sed 's/;/\t/g' DovetailSCNMaker4.all.NOFASTA.gff |sed 's/ID=//g' |awk '$3=="mRNA"' |grep "RA" >SCNGrepMod.gff

sed 's/;/\t/g' G.ellingtonae.gff3 |sed 's/ID=//g' |awk '$3=="gene"' >GELLGrepMod.gff

 less segments.txt |awk 'NR>1'  |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="G.ellingtonae") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line GELLGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
  less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="G.ellingtonae") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line GELLGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
  less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="G.ellingtonae") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
  less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="G.ellingtonae") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8

 less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '{if($5=="G.ellingtonae") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' |awk 'NR>1' |paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf



 bioawk -c fastx '{print $name,length($seq)}' G.ellingtonae.genomic.fa |awk '{print "chr","-",$1,$1,"0",$2,"blue"}'  >G.ellingtonaeKaryotype.conf

 bioawk -c fastx '{print $name,length($seq)}' DovetailSCNMaker4.genome.fasta |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >H.glycinesKaryotype.conf


 awk '{print $1}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' G.ellingtonaeKaryotype.conf >>tmpKaryotype.conf1";done >G.ellingtonaeKaryotype.sh
  sh G.ellingtonaeKaryotype.sh

 awk '{print $4}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' H.glycinesKaryotype.conf >>tmpKaryotype.conf2";done >H.glycinesKaryotype.sh
 sh H.glycinesKaryotype.sh
 cat <(sort tmpKaryotype.conf1 |uniq) <(sort tmpKaryotype.conf2 |uniq) >karyotype.conf
[remkv6@condo011 02_circos]$ circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_15,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_4,Scaffold_58,Scaffold_60,Scaffold_61,Scaffold_62,Scaffold_63,Scaffold_66,Scaffold_69,Scaffold_70,Scaffold_71,Scaffold_72,Scaffold_73,Scaffold_80,Scaffold_90 -static_rx Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_15,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_4,Scaffold_58,Scaffold_60,Scaffold_61,Scaffold_62,Scaffold_63,Scaffold_66,Scaffold_69,Scaffold_70,Scaffold_71,Scaffold_72,Scaffold_73,Scaffold_80,Scaffold_90
calculating round 0
report round 0 minimize init 56411 final 6983 change 87.62%
calculating round 1
report round 1 minimize init 6983 final 4426 change 36.62%
calculating round 2
report round 2 minimize init 4426 final 2363 change 46.61%
calculating round 3
report round 3 minimize init 2363 final 1947 change 17.60%
scorereport init 56411 final 1947 change 96.55%
chromosomes_order = Scaffold_73,MEIZ01000131,MEIZ01000005,MEIZ01000002,MEIZ01000069,MEIZ01000067,MEIZ01000095,MEIZ01000062,MEIZ01000027,MEIZ01000008,Scaffold_3,MEIZ01000051,MEIZ01000063,MEIZ01000256,MEIZ01000031,MEIZ01000137,MEIZ01000012,MEIZ01000026,MEIZ01000146,MEIZ01000252,Scaffold_103,MEIZ01000080,MEIZ01000126,MEIZ01000325,MEIZ01000127,MEIZ01000221,MEIZ01000052,MEIZ01000076,MEIZ01000017,MEIZ01000090,MEIZ01000009,MEIZ01000158,MEIZ01000019,MEIZ01000033,MEIZ01000113,MEIZ01000138,MEIZ01000187,Scaffold_60,MEIZ01000024,MEIZ01000125,MEIZ01000010,MEIZ01000038,MEIZ01000313,MEIZ01000241,MEIZ01000194,MEIZ01000189,MEIZ01000007,MEIZ01000021,MEIZ01000098,MEIZ01000330,Scaffold_44,MEIZ01000016,MEIZ01000037,MEIZ01000060,MEIZ01000032,MEIZ01000014,MEIZ01000050,Scaffold_25,MEIZ01000084,MEIZ01000040,MEIZ01000004,MEIZ01000025,MEIZ01000094,MEIZ01000030,Scaffold_80,MEIZ01000116,MEIZ01000023,MEIZ01000182,MEIZ01000217,MEIZ01000018,MEIZ01000160,MEIZ01000022,MEIZ01000329,MEIZ01000292,MEIZ01000177,MEIZ01000467,MEIZ01000391,MEIZ01000013,MEIZ01000039,MEIZ01000059,Scaffold_19,MEIZ01000003,MEIZ01000107,MEIZ01000155,MEIZ01000164,MEIZ01000045,MEIZ01000091,MEIZ01000015,MEIZ01000036,MEIZ01000035,MEIZ01000001,Scaffold_104,Scaffold_70,Scaffold_61,Scaffold_72,Scaffold_13,Scaffold_15,MEIZ01000078,MEIZ01000061,MEIZ01000077,Scaffold_62,Scaffold_66,Scaffold_63,MEIZ01000006,MEIZ01000057,Scaffold_4,Scaffold_69,Scaffold_71,MEIZ01000142,Scaffold_58,MEIZ01000361,Scaffold_90





 #circos.conf
 #############################################################################
 karyotype = ./karyotype.conf
 chromosomes_units = 100000
   <<include ideogram.conf>>
   <<include ticks.conf>>
   <<include bands.conf>>

   <links>
   <link>
     file=SyntenicRibbons.conf
     radius = 0.94r
     bezier_radius = 0.1r
     thickness = 1
     ribbon = yes
   </link>
   </links>



 <image>
   <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
 angle_offset* = -46
 </image>
 <<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
  <<include ./housekeeping.conf>>
 chromosomes_order = Scaffold_3,scaffold_85,scaffold_25,scaffold_110,scaffold_102,scaffold_222,scaffold_539,scaffold_72,scaffold_53,scaffold_87,scaffold_325,scaffold_444,scaffold_8,scaffold_43,scaffold_278,scaffold_623,scaffold_447,scaffold_412,scaffold_111,scaffold_163,scaffold_1578,scaffold_2159,scaffold_139,scaffold_41,scaffold_930,scaffold_146,scaffold_74,scaffold_619,scaffold_596,scaffold_91,scaffold_390,scaffold_1503,scaffold_233,scaffold_564,scaffold_483,scaffold_296,Scaffold_73,scaffold_826,scaffold_92,scaffold_34,scaffold_311,scaffold_224,scaffold_526,scaffold_196,scaffold_410,scaffold_30,scaffold_145,scaffold_499,scaffold_240,scaffold_47,scaffold_67,scaffold_15,scaffold_231,scaffold_293,scaffold_125,scaffold_257,scaffold_327,scaffold_27,scaffold_271,scaffold_446,scaffold_561,scaffold_226,scaffold_5,scaffold_42,scaffold_1020,scaffold_23,Scaffold_60,scaffold_147,scaffold_90,scaffold_1181,scaffold_155,scaffold_100,scaffold_6,scaffold_38,scaffold_14,scaffold_13,scaffold_20,scaffold_704,scaffold_11,scaffold_131,scaffold_171,scaffold_180,scaffold_50,scaffold_223,scaffold_103,scaffold_71,scaffold_194,scaffold_616,scaffold_573,Scaffold_25,scaffold_22,scaffold_381,scaffold_588,scaffold_3,scaffold_285,scaffold_785,scaffold_323,scaffold_540,scaffold_369,scaffold_193,scaffold_18,scaffold_228,scaffold_170,scaffold_148,scaffold_1364,scaffold_402,scaffold_1325,scaffold_373,scaffold_48,scaffold_31,scaffold_300,scaffold_114,scaffold_84,scaffold_557,scaffold_414,scaffold_107,scaffold_287,scaffold_136,Scaffold_44,scaffold_2,scaffold_443,scaffold_404,scaffold_109,scaffold_159,scaffold_527,scaffold_377,scaffold_260,scaffold_117,scaffold_112,scaffold_12,scaffold_81,scaffold_79,scaffold_1389,scaffold_40,scaffold_214,scaffold_309,scaffold_1446,scaffold_73,scaffold_77,Scaffold_19,scaffold_19,scaffold_328,scaffold_253,scaffold_683,scaffold_119,scaffold_207,scaffold_910,scaffold_335,scaffold_552,scaffold_168,scaffold_54,scaffold_56,scaffold_123,scaffold_441,scaffold_59,scaffold_153,scaffold_374,scaffold_246,scaffold_618,scaffold_143,scaffold_209,Scaffold_103,scaffold_129,scaffold_217,scaffold_44,scaffold_69,scaffold_382,scaffold_211,scaffold_329,scaffold_1166,scaffold_653,scaffold_399,scaffold_9,scaffold_284,scaffold_45,scaffold_36,scaffold_63,scaffold_545,Scaffold_80,scaffold_178,scaffold_505,scaffold_611,scaffold_160,scaffold_21,scaffold_68,scaffold_529,scaffold_530,scaffold_236,scaffold_1244,Scaffold_104,scaffold_94,scaffold_52,scaffold_225,scaffold_1394,scaffold_210,scaffold_29,scaffold_58,scaffold_386,Scaffold_13,scaffold_1,scaffold_35,Scaffold_70,Scaffold_51,scaffold_126,Scaffold_90

 #############################################################################

 ticks.conf
 ###############################################################################
 show_ticks          = yes
show_tick_labels    = yes
<ticks>
    radius           = 1r
    color            = black
    thickness        = 10p
    multiplier       = 1e-5
    format           = %d
 <tick>
    spacing        = 10u
    size           = 25p
    show_label     = yes
    label_size     = 25p
    label_offset   = 10p
    format         = %d
  </tick>

</ticks>

 ###############################################################################

 bands.conf
 ###############################################################################
 <bands>
    show_bands = yes
    fill_bands = yes
    band_transparency = 4
 </bands>
 ###############################################################################

ideogram.conf
###############################################################################
<ideogram>
  <spacing>
    default = 0.006r
    break   = 30u
    axis_break_at_edge = yes
    axis_break         = yes
    axis_break_style   = 2
    <break_style 1>
          stroke_color     = black
          thickness        = 0.45r
          stroke_thickness = 2p
    </break>
    <break_style 2>
          stroke_color     = black
          stroke_thickness = 5p
          thickness        = 4r
    </break>
  </spacing>
  radius           = 0.84r
  thickness        = 80p
  fill             = yes
  stroke_color     = white
  stroke_thickness = 4p
  fill_color       = black
  show_label       = yes
  label_font       = bold
  label_size       = 16
  label_parallel   = no

  label_radius = dims(ideogram,radius_outer) + 0.06r
</ideogram>
###############################################################################

 cp /work/GIF/software/programs/circos/0.69-4/etc/housekeeping.conf .

```


### G. rostochiensis
```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/03_rostochiensis
 #grabbing only primary isoform protein fastas
 ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/DovetailSCNMaker4.all.maker.proteins.Isoform1Only.fasta
 ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/globodera_rostochiensis.PRJEB13504.WBPS10.protein.Isoform1Only.fasta

 #/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/
 module load orthofinder
 orthofinder -f 01_rostochiensis

/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/03_rostochiensis/01_iadhore
less ../Results_Jun21/Orthogroups.txt |tr " " "\n" |awk -v Family=0 '{if(substr($1,1,2)=="OG") {Family=$1} else {print $0,Family}}' |sed 's/://g' |cat - <(dos2unix ../Results_Jun21/Orthogroups_UnassignedGenes.csv|awk 'NR>1 {print $2,$1}') |tr " " "\t" >Orthologues.list
#had to modify the gene names for ellingtonae
sed -i 's/Gell//g' Orthologues.list
less Orthologues.list |sed 's/\.t1//g' >fixed.Orthologues.list

ln -s ../../../12_MakerGenesOrthofinder/globodera_rostochiensis.PRJEB13504.WBPS10.annotations.gff3
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff
mkdir query
less ../DovetailSCNMaker4.all.NOFASTA.gff| awk '$3=="mRNA"' |sed 's/;/\t/g' |sed 's/ID=//g' |grep "RA"|awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini

 mkdir subject
 #make sure to get only the primary transcript
 less ../globodera_rostochiensis.PRJEB13504.WBPS10.annotations.gff3 |awk '$3=="gene"' |sed 's/;/\t/g' |sed 's/ID=//g' |awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'
 sed -i 's/ .*//g' *.lst
 ls *lst >input.txt
 paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
 cd ..
 cat query/query.ini subject/subject.ini |tr "\t" " " >iadhore.ini

 blast_table=Orthologues.list
 table_type=family
 prob_cutoff=0.001
 anchor_points=3
 number_of_threads=16
 visualizeAlignment=true
 output_path= output
 alignment_method=gg2
 gap_size=15
 cluster_gap=20
 level_2_only=true
 q_value=.05
 ```
 Circos setup
 ```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/03_rostochiensis/02_circos
ln -s ../01_iadhore/output/segments.txt
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.genome.fasta
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff
ln -s ../../../12_MakerGenesOrthofinder/globodera_rostochiensis.PRJEB13504.WBPS10.annotations.gff3
ln -s ../../../12_MakerGenesOrthofinder/globodera_rostochiensis.PRJEB13504.WBPS10.genomic.fa

sed 's/;/\t/g' DovetailSCNMaker4.all.NOFASTA.gff |sed 's/ID=//g' |awk '$3=="mRNA"' |grep "RA" >SCNGrepMod.gff

less globodera_rostochiensis.PRJEB13504.WBPS10.annotations.gff3 |awk '$3=="gene"' |sed 's/;/\t/g' |sed 's/ID=gene://g' >GrosGrepMod.gff


 less segments.txt |awk 'NR>1'  |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="G.rostochiensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line GrosGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
  less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="G.rostochiensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line GrosGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
  less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="G.rostochiensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
  less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="G.rostochiensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8

 less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '{if($5=="G.rostochiensis") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' |awk 'NR>1' |paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf



 bioawk -c fastx '{print $name,length($seq)}' globodera_rostochiensis.PRJEB13504.WBPS10.genomic.fa |awk '{print "chr","-",$1,$1,"0",$2,"blue"}'  >G.rostochiensisKaryotype.conf

 bioawk -c fastx '{print $name,length($seq)}' DovetailSCNMaker4.genome.fasta |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >H.glycinesKaryotype.conf


 awk '{print $1}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' G.rostochiensisKaryotype.conf >>tmpKaryotype.conf1";done >G.rostochiensisKaryotype.sh
  sh G.rostochiensisKaryotype.sh

 awk '{print $4}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' H.glycinesKaryotype.conf >>tmpKaryotype.conf2";done >H.glycinesKaryotype.sh
 sh H.glycinesKaryotype.sh
 cat <(sort tmpKaryotype.conf1 |uniq) <(sort tmpKaryotype.conf2 |uniq) >karyotype.conf
 circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_15,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_4,Scaffold_58,Scaffold_60,Scaffold_61,Scaffold_62,Scaffold_63,Scaffold_66,Scaffold_69,Scaffold_70,Scaffold_71,Scaffold_72,Scaffold_73,Scaffold_80,Scaffold_90 -static_rx Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_15,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_4,Scaffold_58,Scaffold_60,Scaffold_61,Scaffold_62,Scaffold_63,Scaffold_66,Scaffold_69,Scaffold_70,Scaffold_71,Scaffold_72,Scaffold_73,Scaffold_80,Scaffold_90
 calculating round 0
 report round 0 minimize init 58927 final 8622 change 85.37%
 calculating round 1
 report round 1 minimize init 8622 final 6986 change 18.97%
 calculating round 2
 report round 2 minimize init 6986 final 5499 change 21.29%
 calculating round 3
 report round 3 minimize init 5499 final 4400 change 19.99%
 scorereport init 58927 final 4400 change 92.53%
 chromosomes_order = Scaffold_73,GROS_00421,GROS_00322,GROS_00228,GROS_00611,GROS_00287,GROS_00072,GROS_00410,GROS_00248,GROS_00016,GROS_00075,GROS_00279,GROS_00059,GROS_00261,GROS_00375,GROS_00002,GROS_00066,GROS_00451,GROS_00005,GROS_00858,GROS_00744,GROS_00505,GROS_00184,GROS_00137,GROS_00176,GROS_00321,GROS_00039,GROS_00277,GROS_00339,GROS_00168,GROS_00175,Scaffold_3,GROS_00238,GROS_00093,GROS_00179,GROS_00054,GROS_00210,GROS_00195,GROS_00035,GROS_00219,GROS_00267,GROS_00288,GROS_00290,GROS_00056,GROS_00098,GROS_00350,GROS_00538,GROS_00701,GROS_00132,GROS_00298,GROS_00198,GROS_00007,GROS_00124,GROS_00150,GROS_00138,GROS_00240,GROS_00094,GROS_00047,GROS_00614,GROS_00069,GROS_00068,GROS_00009,GROS_00110,Scaffold_60,GROS_00006,GROS_00028,GROS_00114,Scaffold_50,GROS_00025,GROS_00134,GROS_00246,GROS_00024,GROS_00431,GROS_00330,GROS_00064,GROS_00163,GROS_00499,GROS_00376,GROS_00190,GROS_00088,GROS_00011,GROS_00318,GROS_00060,GROS_00033,GROS_00020,GROS_00469,GROS_00063,GROS_00101,GROS_00386,GROS_00026,Scaffold_80,GROS_00500,GROS_00147,GROS_00247,GROS_00084,GROS_00360,GROS_00578,GROS_00346,GROS_00092,GROS_00249,GROS_00065,GROS_00099,GROS_00188,GROS_00251,GROS_00543,GROS_00524,GROS_00764,GROS_00120,GROS_00095,GROS_00332,GROS_00873,GROS_00472,GROS_00151,GROS_00270,GROS_00342,GROS_00651,GROS_00111,GROS_00012,Scaffold_25,GROS_00051,GROS_00042,GROS_00019,GROS_00479,GROS_00004,GROS_00340,GROS_00057,GROS_00391,GROS_00519,GROS_00468,GROS_00022,GROS_00036,GROS_00731,GROS_00214,GROS_00062,GROS_00960,GROS_00123,GROS_00070,GROS_01560,GROS_00378,GROS_00076,GROS_00336,GROS_00183,GROS_00032,GROS_00394,GROS_00121,Scaffold_44,GROS_00073,GROS_00215,GROS_00485,GROS_00044,GROS_00031,GROS_00126,GROS_00374,GROS_00136,GROS_00234,GROS_00324,GROS_00638,GROS_00097,GROS_00001,GROS_00055,GROS_00034,GROS_00199,GROS_00021,GROS_00105,GROS_01087,Scaffold_103,GROS_00038,GROS_00018,GROS_00052,GROS_00944,GROS_00128,GROS_00118,GROS_00408,GROS_00061,GROS_00135,GROS_00335,GROS_00087,GROS_00241,GROS_00041,GROS_00158,GROS_00815,GROS_00079,GROS_00107,GROS_00191,GROS_00233,GROS_00293,GROS_00122,GROS_00160,GROS_00522,Scaffold_19,GROS_00010,GROS_00525,GROS_00265,GROS_00366,GROS_00212,GROS_00555,GROS_00362,GROS_00222,GROS_00385,GROS_00882,GROS_00819,GROS_00177,GROS_00131,GROS_00907,GROS_01083,GROS_00268,GROS_00918,GROS_00169,GROS_00678,GROS_00145,GROS_00148,GROS_00648,GROS_00296,GROS_00014,GROS_00390,GROS_00345,GROS_00077,Scaffold_104,GROS_00480,GROS_00493,GROS_00027,GROS_00141,GROS_00125,GROS_00078,GROS_00517,GROS_00037,GROS_00231,GROS_01093,GROS_00294,GROS_00113,GROS_00628,GROS_00252,GROS_00225,Scaffold_13,GROS_00040,GROS_00008,GROS_00003,GROS_00043,GROS_00481,GROS_00312,Scaffold_72,Scaffold_61,GROS_00159,Scaffold_66,GROS_00029,Scaffold_90





 #circos.conf
 #############################################################################
 karyotype = ./karyotype.conf
 chromosomes_units = 100000
   <<include ideogram.conf>>
   <<include ticks.conf>>
   <<include bands.conf>>

   <links>
   <link>
     file=SyntenicRibbons.conf
     radius = 0.94r
     bezier_radius = 0.1r
     thickness = 1
     ribbon = yes
   </link>
   </links>



 <image>
   <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
 angle_offset* = -46
 </image>
 <<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
  <<include ./housekeeping.conf>>
 chromosomes_order = Scaffold_3,scaffold_85,scaffold_25,scaffold_110,scaffold_102,scaffold_222,scaffold_539,scaffold_72,scaffold_53,scaffold_87,scaffold_325,scaffold_444,scaffold_8,scaffold_43,scaffold_278,scaffold_623,scaffold_447,scaffold_412,scaffold_111,scaffold_163,scaffold_1578,scaffold_2159,scaffold_139,scaffold_41,scaffold_930,scaffold_146,scaffold_74,scaffold_619,scaffold_596,scaffold_91,scaffold_390,scaffold_1503,scaffold_233,scaffold_564,scaffold_483,scaffold_296,Scaffold_73,scaffold_826,scaffold_92,scaffold_34,scaffold_311,scaffold_224,scaffold_526,scaffold_196,scaffold_410,scaffold_30,scaffold_145,scaffold_499,scaffold_240,scaffold_47,scaffold_67,scaffold_15,scaffold_231,scaffold_293,scaffold_125,scaffold_257,scaffold_327,scaffold_27,scaffold_271,scaffold_446,scaffold_561,scaffold_226,scaffold_5,scaffold_42,scaffold_1020,scaffold_23,Scaffold_60,scaffold_147,scaffold_90,scaffold_1181,scaffold_155,scaffold_100,scaffold_6,scaffold_38,scaffold_14,scaffold_13,scaffold_20,scaffold_704,scaffold_11,scaffold_131,scaffold_171,scaffold_180,scaffold_50,scaffold_223,scaffold_103,scaffold_71,scaffold_194,scaffold_616,scaffold_573,Scaffold_25,scaffold_22,scaffold_381,scaffold_588,scaffold_3,scaffold_285,scaffold_785,scaffold_323,scaffold_540,scaffold_369,scaffold_193,scaffold_18,scaffold_228,scaffold_170,scaffold_148,scaffold_1364,scaffold_402,scaffold_1325,scaffold_373,scaffold_48,scaffold_31,scaffold_300,scaffold_114,scaffold_84,scaffold_557,scaffold_414,scaffold_107,scaffold_287,scaffold_136,Scaffold_44,scaffold_2,scaffold_443,scaffold_404,scaffold_109,scaffold_159,scaffold_527,scaffold_377,scaffold_260,scaffold_117,scaffold_112,scaffold_12,scaffold_81,scaffold_79,scaffold_1389,scaffold_40,scaffold_214,scaffold_309,scaffold_1446,scaffold_73,scaffold_77,Scaffold_19,scaffold_19,scaffold_328,scaffold_253,scaffold_683,scaffold_119,scaffold_207,scaffold_910,scaffold_335,scaffold_552,scaffold_168,scaffold_54,scaffold_56,scaffold_123,scaffold_441,scaffold_59,scaffold_153,scaffold_374,scaffold_246,scaffold_618,scaffold_143,scaffold_209,Scaffold_103,scaffold_129,scaffold_217,scaffold_44,scaffold_69,scaffold_382,scaffold_211,scaffold_329,scaffold_1166,scaffold_653,scaffold_399,scaffold_9,scaffold_284,scaffold_45,scaffold_36,scaffold_63,scaffold_545,Scaffold_80,scaffold_178,scaffold_505,scaffold_611,scaffold_160,scaffold_21,scaffold_68,scaffold_529,scaffold_530,scaffold_236,scaffold_1244,Scaffold_104,scaffold_94,scaffold_52,scaffold_225,scaffold_1394,scaffold_210,scaffold_29,scaffold_58,scaffold_386,Scaffold_13,scaffold_1,scaffold_35,Scaffold_70,Scaffold_51,scaffold_126,Scaffold_90

 #############################################################################

 ticks.conf
 ###############################################################################
 show_ticks          = yes
show_tick_labels    = yes
<ticks>
    radius           = 1r
    color            = black
    thickness        = 10p
    multiplier       = 1e-5
    format           = %d
 <tick>
    spacing        = 10u
    size           = 25p
    show_label     = yes
    label_size     = 25p
    label_offset   = 10p
    format         = %d
  </tick>

</ticks>

 ###############################################################################

 bands.conf
 ###############################################################################
 <bands>
    show_bands = yes
    fill_bands = yes
    band_transparency = 4
 </bands>
 ###############################################################################

ideogram.conf
###############################################################################
<ideogram>
  <spacing>
    default = 0.006r
    break   = 30u
    axis_break_at_edge = yes
    axis_break         = yes
    axis_break_style   = 2
    <break_style 1>
          stroke_color     = black
          thickness        = 0.45r
          stroke_thickness = 2p
    </break>
    <break_style 2>
          stroke_color     = black
          stroke_thickness = 5p
          thickness        = 4r
    </break>
  </spacing>
  radius           = 0.84r
  thickness        = 80p
  fill             = yes
  stroke_color     = white
  stroke_thickness = 4p
  fill_color       = black
  show_label       = yes
  label_font       = bold
  label_size       = 16
  label_parallel   = no

  label_radius = dims(ideogram,radius_outer) + 0.06r
</ideogram>
###############################################################################

 cp /work/GIF/software/programs/circos/0.69-4/etc/housekeeping.conf .
```
### Meloidogyne hapla
 ```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/04_hapla
  #grabbing only primary isoform protein fastas
  ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/DovetailSCNMaker4.all.maker.proteins.Isoform1Only.fasta
  ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/meloidogyne_hapla.PRJNA29083.WBPS10.protein.fa

  #/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/
  module load orthofinder
  orthofinder -f 04_hapla

#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/04_hapla/01_iadhore
 less ../Results_Jun21/Orthogroups.txt |tr " " "\n" |awk -v Family=0 '{if(substr($1,1,2)=="OG") {Family=$1} else {print $0,Family}}' |sed 's/://g' |cat - <(dos2unix ../Results_Jun21/Orthogroups_UnassignedGenes.csv|awk 'NR>1 {print $2,$1}') |tr " " "\t" >Orthologues.list


 ln -s ../../../12_MakerGenesOrthofinder/globodera_rostochiensis.PRJEB13504.WBPS10.annotations.gff3
 ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff
 mkdir query
 less ../DovetailSCNMaker4.all.NOFASTA.gff| awk '$3=="mRNA"' |sed 's/;/\t/g' |sed 's/ID=//g' |grep "RA"|awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'
 sed -i 's/ .*//g' *.lst
 ls *lst >input.txt
 paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini

  mkdir subject
  #make sure to get only the primary transcript
less ../meloidogyne_incognita.PRJEA28837.WBPS10.annotations.gff3 |awk '$3=="gene"' |sed 's/;/\t/g' |sed 's/ID=gene://g' |awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'
  sed -i 's/ .*//g' *.lst
  ls *lst >input.txt
  paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
  cd ..
  cat query/query.ini subject/subject.ini |tr "\t" " " >iadhore.ini

  blast_table=Orthologues.list
  table_type=family
  prob_cutoff=0.001
  anchor_points=3
  number_of_threads=16
  visualizeAlignment=true
  output_path= output
  alignment_method=gg2
  gap_size=15
  cluster_gap=20
  level_2_only=true
  q_value=.05
  ```
  Circos setup
  ```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/04_hapla/02_circos
 ln -s ../01_iadhore/output/segments.txt
 ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.genome.fasta
 ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff
 ln -s ../../../12_MakerGenesOrthofinder/meloidogyne_hapla.PRJNA29083.WBPS10.annotations.gff3
 ln -s ../../../12_MakerGenesOrthofinder/meloidogyne_hapla.PRJNA29083.WBPS10.genomic.fa

 sed 's/;/\t/g' DovetailSCNMaker4.all.NOFASTA.gff |sed 's/ID=//g' |awk '$3=="mRNA"' |grep "RA" >SCNGrepMod.gff

 less globodera_rostochiensis.PRJEB13504.WBPS10.annotations.gff3 |awk '$3=="gene"' |sed 's/;/\t/g' |sed 's/ID=gene://g' >MhapGrepMod.gff


  less segments.txt |awk 'NR>1'  |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="M.hapla") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line MhapGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
   less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="M.hapla") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line MhapGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
   less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="M.hapla") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
   less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="M.hapla") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8

  less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '{if($5=="M.hapla") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' |awk 'NR>1' |paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf



  bioawk -c fastx '{print $name,length($seq)}' meloidogyne_hapla.PRJNA29083.WBPS10.genomic.fa |awk '{print "chr","-",$1,$1,"0",$2,"blue"}'  >M.haplaKaryotype.conf

  bioawk -c fastx '{print $name,length($seq)}' DovetailSCNMaker4.genome.fasta |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >H.glycinesKaryotype.conf


  awk '{print $1}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' M.haplaKaryotype.conf >>tmpKaryotype.conf1";done >M.haplaKaryotype.sh
   sh M.haplaKaryotype.sh

  awk '{print $4}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' H.glycinesKaryotype.conf >>tmpKaryotype.conf2";done >H.glycinesKaryotype.sh
  sh H.glycinesKaryotype.sh
  cat <(sort tmpKaryotype.conf1 |uniq) <(sort tmpKaryotype.conf2 |uniq) >karyotype.conf
  ../../03_rostochiensis/02_circos/circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_15,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_4,Scaffold_58,Scaffold_60,Scaffold_61,Scaffold_62,Scaffold_63,Scaffold_66,Scaffold_69,Scaffold_70,Scaffold_71,Scaffold_72,Scaffold_73,Scaffold_80,Scaffold_90 -static_rx Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_15,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_4,Scaffold_58,Scaffold_60,Scaffold_61,Scaffold_62,Scaffold_63,Scaffold_66,Scaffold_69,Scaffold_70,Scaffold_71,Scaffold_72,Scaffold_73,Scaffold_80,Scaffold_90
calculating round 0
report round 0 minimize init 3105 final 397 change 87.21%
calculating round 1
report round 1 minimize init 397 final 268 change 32.49%
calculating round 2
report round 2 minimize init 268 final 132 change 50.75%
calculating round 3
report round 3 minimize init 132 final 74 change 43.94%
scorereport init 3105 final 74 change 97.62%
chromosomes_order = Scaffold_3,MhA1_Contig737,MhA1_Contig730,MhA1_Contig542,MhA1_Contig388,MhA1_Contig48,MhA1_Contig304,MhA1_Contig42,MhA1_Contig116,MhA1_Contig168,MhA1_Contig290,MhA1_Contig120,MhA1_Contig121,MhA1_Contig707,MhA1_Contig1,MhA1_Contig320,Scaffold_19,MhA1_Contig765,MhA1_Contig202,MhA1_Contig508,MhA1_Contig754,MhA1_Contig756,MhA1_Contig7,MhA1_Contig9,MhA1_Contig0,MhA1_Contig580,MhA1_Contig130,MhA1_Contig53,MhA1_Contig1668,MhA1_Contig855,Scaffold_25,MhA1_Contig884,MhA1_Contig684,MhA1_Contig656,MhA1_Contig253,MhA1_Contig1534,MhA1_Contig800,MhA1_Contig2086,MhA1_Contig1246,MhA1_Contig271,MhA1_Contig232,MhA1_Contig764,Scaffold_73,MhA1_Contig487,MhA1_Contig990,MhA1_Contig1305,MhA1_Contig1078,MhA1_Contig1555,MhA1_Contig131,MhA1_Contig2,MhA1_Contig877,MhA1_Contig325,MhA1_Contig30,MhA1_Contig2034,MhA1_Contig590,Scaffold_103,MhA1_Contig78,MhA1_Contig1462,MhA1_Contig1067,MhA1_Contig417,MhA1_Contig952,MhA1_Contig353,MhA1_Contig309,MhA1_Contig118,MhA1_Contig2021,MhA1_Contig933,Scaffold_60,MhA1_Contig671,MhA1_Contig492,MhA1_Contig912,MhA1_Contig1987,MhA1_Contig165,MhA1_Contig849,MhA1_Contig109,MhA1_Contig108,Scaffold_44,MhA1_Contig342,MhA1_Contig25,MhA1_Contig26,MhA1_Contig70,MhA1_Contig29,MhA1_Contig194,MhA1_Contig133,MhA1_Contig20,Scaffold_80,MhA1_Contig649,MhA1_Contig2087,MhA1_Contig484,MhA1_Contig1189,MhA1_Contig1000,MhA1_Contig1391,MhA1_Contig134,MhA1_Contig213,MhA1_Contig196,Scaffold_13,MhA1_Contig123,MhA1_Contig113,MhA1_Contig2334,MhA1_Contig874,Scaffold_104,MhA1_Contig354,Scaffold_61,MhA1_Contig1377,Scaffold_70






  #circos.conf
  #############################################################################
  karyotype = ./karyotype.conf
  chromosomes_units = 100000
    <<include ideogram.conf>>
    <<include ticks.conf>>
    <<include bands.conf>>

    <links>
    <link>
      file=SyntenicRibbons.conf
      radius = 0.94r
      bezier_radius = 0.1r
      thickness = 1
      ribbon = yes
    </link>
    </links>



  <image>
    <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
  angle_offset* = -46
  </image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
   <<include ./housekeeping.conf>>
   chromosomes_order = Scaffold_3,MhA1_Contig737,MhA1_Contig730,MhA1_Contig542,MhA1_Contig388,MhA1_Contig48,MhA1_Contig304,MhA1_Contig42,MhA1_Contig116,MhA1_Contig168,MhA1_Contig290,MhA1_Contig120,MhA1_Contig121,MhA1_Contig707,MhA1_Contig1,MhA1_Contig320,Scaffold_19,MhA1_Contig765,MhA1_Contig202,MhA1_Contig508,MhA1_Contig754,MhA1_Contig756,MhA1_Contig7,MhA1_Contig9,MhA1_Contig0,MhA1_Contig580,MhA1_Contig130,MhA1_Contig53,MhA1_Contig1668,MhA1_Contig855,Scaffold_25,MhA1_Contig884,MhA1_Contig684,MhA1_Contig656,MhA1_Contig253,MhA1_Contig1534,MhA1_Contig800,MhA1_Contig2086,MhA1_Contig1246,MhA1_Contig271,MhA1_Contig232,MhA1_Contig764,Scaffold_73,MhA1_Contig487,MhA1_Contig990,MhA1_Contig1305,MhA1_Contig1078,MhA1_Contig1555,MhA1_Contig131,MhA1_Contig2,MhA1_Contig877,MhA1_Contig325,MhA1_Contig30,MhA1_Contig2034,MhA1_Contig590,Scaffold_103,MhA1_Contig78,MhA1_Contig1462,MhA1_Contig1067,MhA1_Contig417,MhA1_Contig952,MhA1_Contig353,MhA1_Contig309,MhA1_Contig118,MhA1_Contig2021,MhA1_Contig933,Scaffold_60,MhA1_Contig671,MhA1_Contig492,MhA1_Contig912,MhA1_Contig1987,MhA1_Contig165,MhA1_Contig849,MhA1_Contig109,MhA1_Contig108,Scaffold_44,MhA1_Contig342,MhA1_Contig25,MhA1_Contig26,MhA1_Contig70,MhA1_Contig29,MhA1_Contig194,MhA1_Contig133,MhA1_Contig20,Scaffold_80,MhA1_Contig649,MhA1_Contig2087,MhA1_Contig484,MhA1_Contig1189,MhA1_Contig1000,MhA1_Contig1391,MhA1_Contig134,MhA1_Contig213,MhA1_Contig196,Scaffold_13,MhA1_Contig123,MhA1_Contig113,MhA1_Contig2334,MhA1_Contig874,Scaffold_104,MhA1_Contig354,Scaffold_61,MhA1_Contig1377,Scaffold_70


  #############################################################################

  ticks.conf
  ###############################################################################
  show_ticks          = yes
 show_tick_labels    = yes
 <ticks>
     radius           = 1r
     color            = black
     thickness        = 10p
     multiplier       = 1e-5
     format           = %d
  <tick>
     spacing        = 10u
     size           = 25p
     show_label     = yes
     label_size     = 25p
     label_offset   = 10p
     format         = %d
   </tick>

 </ticks>

  ###############################################################################

  bands.conf
  ###############################################################################
  <bands>
     show_bands = yes
     fill_bands = yes
     band_transparency = 4
  </bands>
  ###############################################################################

 ideogram.conf
 ###############################################################################
 <ideogram>
   <spacing>
     default = 0.006r
     break   = 30u
     axis_break_at_edge = yes
     axis_break         = yes
     axis_break_style   = 2
     <break_style 1>
           stroke_color     = black
           thickness        = 0.45r
           stroke_thickness = 2p
     </break>
     <break_style 2>
           stroke_color     = black
           stroke_thickness = 5p
           thickness        = 4r
     </break>
   </spacing>
   radius           = 0.84r
   thickness        = 80p
   fill             = yes
   stroke_color     = white
   stroke_thickness = 4p
   fill_color       = black
   show_label       = yes
   label_font       = bold
   label_size       = 16
   label_parallel   = no

   label_radius = dims(ideogram,radius_outer) + 0.06r
 </ideogram>
 ###############################################################################

  cp /work/GIF/software/programs/circos/0.69-4/etc/housekeeping.conf .
  ```

### Meloidogyne incognita
   ```
/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/05_incognita
#grabbing only primary isoform protein fastas
ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/DovetailSCNMaker4.all.maker.proteins.Isoform1Only.fasta
ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/meloidogyne_incognita.PRJEA28837.WBPS10.protein.fa

#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/
module load orthofinder
orthofinder -f 05_incognita

#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/05_incognita/01_iadhore
less ../Results_Jun21/Orthogroups.txt |tr " " "\n" |awk -v Family=0 '{if(substr($1,1,2)=="OG") {Family=$1} else {print $0,Family}}' |sed 's/://g' |cat - <(dos2unix ../Results_Jun21/Orthogroups_UnassignedGenes.csv|awk 'NR>1 {print $2,$1}') |tr " " "\t" >Orthologues.list
less Orthologues.list |grep -v "HetGly" |awk 'length($1)>9 {print $0}' |grep "a" |sed 's/a//g' >1isoform.list
less Orthologues.list |grep -v "HetGly" |awk 'length($1)==9 {print $0}' >2isoform.list
cat <(grep "HetGly" Orthologues.list ) 1isoform.list 2isoform.list |sort -k2,2V >FixedOrtholog.list

ln -s ../../../12_MakerGenesOrthofinder/globodera_rostochiensis.PRJEB13504.WBPS10.annotations.gff3
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff
mkdir query
less ../DovetailSCNMaker4.all.NOFASTA.gff| awk '$3=="mRNA"' |sed 's/;/\t/g' |sed 's/ID=//g' |grep "RA"|awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini

mkdir subject
#make sure to get only the primary transcript
less ../meloidogyne_incognita.PRJEA28837.WBPS10.annotations.gff3 |awk '$3=="gene"' |sed 's/;/\t/g' |sed 's/ID=gene://g' |awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
cd ..
cat query/query.ini subject/subject.ini |tr "\t" " " >iadhore.ini

blast_table=Orthologues.list
table_type=family
prob_cutoff=0.001
anchor_points=3
number_of_threads=16
visualizeAlignment=true
output_path= output
alignment_method=gg2
gap_size=15
cluster_gap=20
level_2_only=true
q_value=.05
```
Circos setup
```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/04_hapla/02_circos
ln -s ../../../12_MakerGenesOrthofinder/meloidogyne_incognita.PRJEA28837.WBPS10.annotations.gff3
ln -s ../../../12_MakerGenesOrthofinder/meloidogyne_incognita.PRJEA28837.WBPS10.genomic.fa
ln -s ../01_iadhore/output/segments.txt
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.genome.fasta
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff

sed 's/;/\t/g' DovetailSCNMaker4.all.NOFASTA.gff |sed 's/ID=//g' |awk '$3=="mRNA"' |grep "RA" >SCNGrepMod.gff

less meloidogyne_incognita.PRJEA28837.WBPS10.annotations.gff3 |awk '$3=="gene"' |sed 's/;/\t/g' |sed 's/ID=gene://g' >MincGrepMod.gff


less segments.txt |awk 'NR>1'  |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="M.incognita") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line MincGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
 less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="M.incognita") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line MincGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
 less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="M.incognita") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
 less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="M.incognita") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8

less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '{if($5=="M.incognita") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' |awk 'NR>1' |paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf



bioawk -c fastx '{print $name,length($seq)}' meloidogyne_incognita.PRJEA28837.WBPS10.genomic.fa |awk '{print "chr","-",$1,$1,"0",$2,"blue"}'  >M.incognitaKaryotype.conf

bioawk -c fastx '{print $name,length($seq)}' DovetailSCNMaker4.genome.fasta |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >H.glycinesKaryotype.conf


awk '{print $1}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' M.haplaKaryotype.conf >>tmpKaryotype.conf1";done >M.haplaKaryotype.sh
 sh M.haplaKaryotype.sh

awk '{print $4}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' H.glycinesKaryotype.conf >>tmpKaryotype.conf2";done >H.glycinesKaryotype.sh
sh H.glycinesKaryotype.sh
cat <(sort tmpKaryotype.conf1 |uniq) <(sort tmpKaryotype.conf2 |uniq) >karyotype.conf
../../03_rostochiensis/02_circos/circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_15,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_4,Scaffold_58,Scaffold_60,Scaffold_61,Scaffold_62,Scaffold_63,Scaffold_66,Scaffold_69,Scaffold_70,Scaffold_71,Scaffold_72,Scaffold_73,Scaffold_80,Scaffold_90 -static_rx Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_15,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_4,Scaffold_58,Scaffold_60,Scaffold_61,Scaffold_62,Scaffold_63,Scaffold_66,Scaffold_69,Scaffold_70,Scaffold_71,Scaffold_72,Scaffold_73,Scaffold_80,Scaffold_90
calculating round 0
report round 0 minimize init 2301 final 291 change 87.35%
calculating round 1
report round 1 minimize init 291 final 171 change 41.24%
calculating round 2
report round 2 minimize init 171 final 72 change 57.89%
calculating round 3
report round 3 minimize init 72 final 39 change 45.83%
scorereport init 2301 final 39 change 98.31%
chromosomes_order = Scaffold_3,MiV1ctg564,MiV1ctg276,MiV1ctg709,MiV1ctg64,MiV1ctg128,MiV1ctg101,MiV1ctg69,MiV1ctg97,MiV1ctg443,MiV1ctg159,MiV1ctg307,MiV1ctg66,MiV1ctg22,MiV1ctg288,MiV1ctg16,Scaffold_44,MiV1ctg23,MiV1ctg10,MiV1ctg1080,MiV1ctg117,MiV1ctg1025,MiV1ctg719,MiV1ctg769,MiV1ctg114,MiV1ctg327,MiV1ctg370,MiV1ctg338,MiV1ctg39,MiV1ctg60,MiV1ctg18,MiV1ctg118,MiV1ctg296,Scaffold_25,MiV1ctg1213,MiV1ctg344,MiV1ctg12,MiV1ctg57,MiV1ctg70,MiV1ctg107,MiV1ctg4,MiV1ctg464,MiV1ctg1068,MiV1ctg33,Scaffold_80,MiV1ctg244,MiV1ctg477,MiV1ctg934,MiV1ctg490,MiV1ctg21,MiV1ctg100,MiV1ctg218,MiV1ctg603,Scaffold_19,MiV1ctg120,MiV1ctg232,MiV1ctg371,MiV1ctg381,MiV1ctg17,MiV1ctg493,MiV1ctg105,Scaffold_73,MiV1ctg349,MiV1ctg77,MiV1ctg199,MiV1ctg431,MiV1ctg471,MiV1ctg174,MiV1ctg365,Scaffold_103,MiV1ctg13,MiV1ctg281,MiV1ctg28,MiV1ctg44,MiV1ctg615,Scaffold_60,MiV1ctg707,MiV1ctg388,MiV1ctg494,MiV1ctg419,MiV1ctg41,MiV1ctg144,Scaffold_13,MiV1ctg273,MiV1ctg400,MiV1ctg164,Scaffold_61







#circos.conf
#############################################################################
karyotype = ./karyotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
chromosomes_order = Scaffold_3,MiV1ctg564,MiV1ctg276,MiV1ctg709,MiV1ctg64,MiV1ctg128,MiV1ctg101,MiV1ctg69,MiV1ctg97,MiV1ctg443,MiV1ctg159,MiV1ctg307,MiV1ctg66,MiV1ctg22,MiV1ctg288,MiV1ctg16,Scaffold_44,MiV1ctg23,MiV1ctg10,MiV1ctg1080,MiV1ctg117,MiV1ctg1025,MiV1ctg719,MiV1ctg769,MiV1ctg114,MiV1ctg327,MiV1ctg370,MiV1ctg338,MiV1ctg39,MiV1ctg60,MiV1ctg18,MiV1ctg118,MiV1ctg296,Scaffold_25,MiV1ctg1213,MiV1ctg344,MiV1ctg12,MiV1ctg57,MiV1ctg70,MiV1ctg107,MiV1ctg4,MiV1ctg464,MiV1ctg1068,MiV1ctg33,Scaffold_80,MiV1ctg244,MiV1ctg477,MiV1ctg934,MiV1ctg490,MiV1ctg21,MiV1ctg100,MiV1ctg218,MiV1ctg603,Scaffold_19,MiV1ctg120,MiV1ctg232,MiV1ctg371,MiV1ctg381,MiV1ctg17,MiV1ctg493,MiV1ctg105,Scaffold_73,MiV1ctg349,MiV1ctg77,MiV1ctg199,MiV1ctg431,MiV1ctg471,MiV1ctg174,MiV1ctg365,Scaffold_103,MiV1ctg13,MiV1ctg281,MiV1ctg28,MiV1ctg44,MiV1ctg615,Scaffold_60,MiV1ctg707,MiV1ctg388,MiV1ctg494,MiV1ctg419,MiV1ctg41,MiV1ctg144,Scaffold_13,MiV1ctg273,MiV1ctg400,MiV1ctg164,Scaffold_61


#############################################################################

ticks.conf
###############################################################################
show_ticks          = yes
show_tick_labels    = yes
<ticks>
   radius           = 1r
   color            = black
   thickness        = 10p
   multiplier       = 1e-5
   format           = %d
<tick>
   spacing        = 10u
   size           = 25p
   show_label     = yes
   label_size     = 25p
   label_offset   = 10p
   format         = %d
 </tick>

</ticks>

###############################################################################

bands.conf
###############################################################################
<bands>
   show_bands = yes
   fill_bands = yes
   band_transparency = 4
</bands>
###############################################################################

ideogram.conf
###############################################################################

 <ideogram>
   <spacing>
     default = 0.006r
     break   = 30u
     axis_break_at_edge = yes
     axis_break         = yes
     axis_break_style   = 2
     <break_style 1>
           stroke_color     = black
           thickness        = 0.45r
           stroke_thickness = 2p
     </break>
     <break_style 2>
           stroke_color     = black
           stroke_thickness = 5p
           thickness        = 4r
     </break>
   </spacing>
   radius           = 0.84r
   thickness        = 80p
   fill             = yes
   stroke_color     = white
   stroke_thickness = 4p
   fill_color       = black
   show_label       = yes
   label_font       = bold
   label_size       = 16
   label_parallel   = no

   label_radius = dims(ideogram,radius_outer) + 0.06r
 </ideogram>
 ###############################################################################

  cp /work/GIF/software/programs/circos/0.69-4/etc/housekeeping.conf .
```

### B. xylophilus
```
/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/06_xylophilus
#grabbing only primary isoform protein fastas
ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/DovetailSCNMaker4.all.maker.proteins.Isoform1Only.fasta
ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/bursaphelenchus_xylophilus.PRJEA64437.WBPS10.protein.Isoform1Only.fasta

#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/
module load orthofinder
orthofinder -f 06_xylophilus

/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/06_xylophilus/01_iadhore
less ../Results_Jun21/Orthogroups.txt |tr " " "\n" |awk -v Family=0 '{if(substr($1,1,2)=="OG") {Family=$1} else {print $0,Family}}' |sed 's/://g' |cat - <(dos2unix ../Results_Jun21/Orthogroups_UnassignedGenes.csv|awk 'NR>1 {print $2,$1}') |tr " " "\t" >Orthologues.list


ln -s ../../../12_MakerGenesOrthofinder/globodera_rostochiensis.PRJEB13504.WBPS10.annotations.gff3
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff
mkdir query
cd query
less ../DovetailSCNMaker4.all.NOFASTA.gff| awk '$3=="mRNA"' |sed 's/;/\t/g' |sed 's/ID=//g' |grep "RA"|awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini
cd ../
mkdir subject
cd subject
#make sure to get only the primary transcript
less ../bursaphelenchus_xylophilus.PRJEA64437.WBPS10.annotations.gff3 |awk '$3=="mRNA"' |sed 's/ID=transcript://g' |sed 's/;/\t/g' |awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}''
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
cd ..
cat query/query.ini subject/subject.ini |tr "\t" " " >iadhore.ini

blast_table=Orthologues.list
table_type=family
prob_cutoff=0.001
anchor_points=3
number_of_threads=16
visualizeAlignment=true
output_path= output
alignment_method=gg2
gap_size=15
cluster_gap=20
level_2_only=true
q_value=.05
```
Circos setup
```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/06_xylophilus/02_circos
ln -s ../../../12_MakerGenesOrthofinder/bursaphelenchus_xylophilus.PRJEA64437.WBPS10.annotations.gff3
ln -s ../../../12_MakerGenesOrthofinder/bursaphelenchus_xylophilus.PRJEA64437.WBPS10.genomic.fa
ln -s ../01_iadhore/output/segments.txt
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.genome.fasta
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff

sed 's/;/\t/g' DovetailSCNMaker4.all.NOFASTA.gff |sed 's/ID=//g' |awk '$3=="mRNA"' |grep "RA" >SCNGrepMod.gff

less ../bursaphelenchus_xylophilus.PRJEA64437.WBPS10.annotations.gff3 |awk '$3=="mRNA"' |sed 's/ID=transcript://g' |sed 's/;/\t/g' >BxylGrepMod.gff

less segments.txt |awk 'NR>1'  |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="B.xylophilus") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line BxylGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="B.xylophilus") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line BxylGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="B.xylophilus") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="B.xylophilus") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8

less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '{if($5=="B.xylophilus") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' |awk 'NR>1' |paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf



bioawk -c fastx '{print $name,length($seq)}' bursaphelenchus_xylophilus.PRJEA64437.WBPS10.genomic.fa |awk '{print "chr","-",$1,$1,"0",$2,"blue"}'  >B.xylophilusKaryotype.conf

bioawk -c fastx '{print $name,length($seq)}' DovetailSCNMaker4.genome.fasta |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >H.glycinesKaryotype.conf


awk '{print $1}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' B.xylophilusKaryotype.conf >>tmpKaryotype.conf1";done >B.xylophilusKaryotype.sh
sh B.xylophilusKaryotype.sh

awk '{print $4}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' H.glycinesKaryotype.conf >>tmpKaryotype.conf2";done >H.glycinesKaryotype.sh
sh H.glycinesKaryotype.sh
cat <(sort tmpKaryotype.conf1 |uniq) <(sort tmpKaryotype.conf2 |uniq) >karyotype.conf
../../03_rostochiensis/02_circos/circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_15,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_4,Scaffold_58,Scaffold_60,Scaffold_61,Scaffold_62,Scaffold_63,Scaffold_66,Scaffold_69,Scaffold_70,Scaffold_71,Scaffold_72,Scaffold_73,Scaffold_80,Scaffold_90 -static_rx Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_15,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_4,Scaffold_58,Scaffold_60,Scaffold_61,Scaffold_62,Scaffold_63,Scaffold_66,Scaffold_69,Scaffold_70,Scaffold_71,Scaffold_72,Scaffold_73,Scaffold_80,Scaffold_90
calculating round 0
report round 0 minimize init 126 final 29 change 76.98%
calculating round 1
report round 1 minimize init 29 final 9 change 68.97%
calculating round 2
report round 2 minimize init 9 final 9 change 0.00%
calculating round 3
report round 3 minimize init 9 final 9 change 0.00%
scorereport init 126 final 9 change 92.86%
chromosomes_order = Scaffold_19,scaffold01144,Scaffold_1,scaffold00422,scaffold00460,Scaffold_104,Scaffold_44,Scaffold_3,Scaffold_25,Scaffold_73,scaffold00055,scaffold01109,scaffold01513,scaffold00116,scaffold00713,Scaffold_60,scaffold01143,scaffold00298,scaffold01653,scaffold01078,Scaffold_103







#circos.conf
#############################################################################
karyotype = ./karyotype.conf
chromosomes_units = 100000
<<include ideogram.conf>>
<<include ticks.conf>>
<<include bands.conf>>

<links>
<link>
 file=SyntenicRibbons.conf
 radius = 0.94r
 bezier_radius = 0.1r
 thickness = 1
 ribbon = yes
</link>
</links>



<image>
<<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
<<include ./housekeeping.conf>>
chromosomes_order = Scaffold_19,scaffold01144,Scaffold_1,scaffold00422,scaffold00460,Scaffold_104,Scaffold_44,Scaffold_3,Scaffold_25,Scaffold_73,scaffold00055,scaffold01109,scaffold01513,scaffold00116,scaffold00713,Scaffold_60,scaffold01143,scaffold00298,scaffold01653,scaffold01078,Scaffold_103


#############################################################################

ticks.conf
###############################################################################
show_ticks          = yes
show_tick_labels    = yes
<ticks>
radius           = 1r
color            = black
thickness        = 10p
multiplier       = 1e-5
format           = %d
<tick>
spacing        = 10u
size           = 25p
show_label     = yes
label_size     = 25p
label_offset   = 10p
format         = %d
</tick>

</ticks>

###############################################################################

bands.conf
###############################################################################
<bands>
show_bands = yes
fill_bands = yes
band_transparency = 4
</bands>
###############################################################################

ideogram.conf
###############################################################################

<ideogram>
<spacing>
  default = 0.006r
  break   = 30u
  axis_break_at_edge = yes
  axis_break         = yes
  axis_break_style   = 2
  <break_style 1>
        stroke_color     = black
        thickness        = 0.45r
        stroke_thickness = 2p
  </break>
  <break_style 2>
        stroke_color     = black
        stroke_thickness = 5p
        thickness        = 4r
  </break>
</spacing>
radius           = 0.84r
thickness        = 80p
fill             = yes
stroke_color     = white
stroke_thickness = 4p
fill_color       = black
show_label       = yes
label_font       = bold
label_size       = 16
label_parallel   = no

label_radius = dims(ideogram,radius_outer) + 0.06r
</ideogram>
###############################################################################

cp /work/GIF/software/programs/circos/0.69-4/etc/housekeeping.conf .
```

### D. destructor
```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/07_destructor
#grabbing only primary isoform protein fastas
ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/DovetailSCNMaker4.all.maker.proteins.Isoform1Only.fasta
ln -s ../../12_MakerGenesOrthofinder/02_Orthofinder1isoform/ditylenchus_destructor.PRJNA312427.WBPS10.protein.fa

#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/
module load orthofinder
orthofinder -f 07_destructor

/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/07_destructor/01_iadhore
less ../Results_Jun21/Orthogroups.txt |tr " " "\n" |awk -v Family=0 '{if(substr($1,1,2)=="OG") {Family=$1} else {print $0,Family}}' |sed 's/://g' |cat - <(dos2unix ../Results_Jun21/Orthogroups_UnassignedGenes.csv|awk 'NR>1 {print $2,$1}') |tr " " "\t" >Orthologues.list


ln -s ../../../12_MakerGenesOrthofinder/ditylenchus_destructor.PRJNA312427.WBPS10.annotations.gff3
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff
mkdir query
cd query
less ../DovetailSCNMaker4.all.NOFASTA.gff| awk '$3=="mRNA"' |sed 's/;/\t/g' |sed 's/ID=//g' |grep "RA"|awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini
cd ../
mkdir subject
cd subject
#make sure to get only the primary transcript
less ../ditylenchus_destructor.PRJNA312427.WBPS10.annotations.gff3 |awk '$3=="gene"' |sed 's/ID=gene://g' |sed 's/;/\t/g' |awk '{print $9$7,$1}' |awk '{print >> $2 ".lst"; close($2)}''
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
cd ..
cat query/query.ini subject/subject.ini |tr "\t" " " >iadhore.ini

blast_table=Orthologues.list
table_type=family
prob_cutoff=0.001
anchor_points=3
number_of_threads=16
visualizeAlignment=true
output_path= output
alignment_method=gg2
gap_size=15
cluster_gap=20
level_2_only=true
q_value=.05
```
Circos setup
```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/15_iadhore/06_xylophilus/02_circos
ln -s ../../../12_MakerGenesOrthofinder/ditylenchus_destructor.PRJNA312427.WBPS10.genomic.fa
ln -s ../../../12_MakerGenesOrthofinder/ditylenchus_destructor.PRJNA312427.WBPS10.annotations.gff3
ln -s ../01_iadhore/output/segments.txt
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.genome.fasta
ln -s ../../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff

sed 's/;/\t/g' DovetailSCNMaker4.all.NOFASTA.gff |sed 's/ID=//g' |awk '$3=="mRNA"' |grep "RA" >SCNGrepMod.gff

less ditylenchus_destructor.PRJNA312427.WBPS10.annotations.gff3 |awk '$3=="gene"' |sed 's/ID=gene://g' |sed 's/;/\t/g' >DdesGrepMod.gff

less segments.txt |awk 'NR>1'  |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="D.destructor") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line DdesGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col3
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="D.destructor") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line DdesGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col4
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="D.destructor") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $5} else {print $4}}' >Col7
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="D.destructor") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -w $line SCNGrepMod.gff; done |awk '{if($7=="+") {print $4} else {print $5}}' >Col8

less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '{if($5=="D.destructor") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' |awk 'NR>1' |paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf



bioawk -c fastx '{print $name,length($seq)}' ditylenchus_destructor.PRJNA312427.WBPS10.genomic.fa |awk '{print "chr","-",$1,$1,"0",$2,"blue"}'  >D.destructorKaryotype.conf

bioawk -c fastx '{print $name,length($seq)}' DovetailSCNMaker4.genome.fasta |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >H.glycinesKaryotype.conf


awk '{print $1}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' D.destructorKaryotype.conf >>tmpKaryotype.conf1";done >D.destructorKaryotype.sh
sh D.destructorKaryotype.sh

awk '{print $4}' SyntenicRibbons.conf|while read line; do echo "awk '\$3==\""$line"\"' H.glycinesKaryotype.conf >>tmpKaryotype.conf2";done >H.glycinesKaryotype.sh
sh H.glycinesKaryotype.sh
cat <(sort tmpKaryotype.conf1 |uniq) <(sort tmpKaryotype.conf2 |uniq) >karyotype.conf
../../03_rostochiensis/02_circos/circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf - Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_15,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_4,Scaffold_58,Scaffold_60,Scaffold_61,Scaffold_62,Scaffold_63,Scaffold_66,Scaffold_69,Scaffold_70,Scaffold_71,Scaffold_72,Scaffold_73,Scaffold_80,Scaffold_90 -static_rx Scaffold_103,Scaffold_104,Scaffold_13,Scaffold_15,Scaffold_19,Scaffold_25,Scaffold_3,Scaffold_44,Scaffold_4,Scaffold_58,Scaffold_60,Scaffold_61,Scaffold_62,Scaffold_63,Scaffold_66,Scaffold_69,Scaffold_70,Scaffold_71,Scaffold_72,Scaffold_73,Scaffold_80,Scaffold_90
calculating round 0
report round 0 minimize init 1668 final 227 change 86.39%
calculating round 1
report round 1 minimize init 227 final 126 change 44.49%
calculating round 2
report round 2 minimize init 126 final 78 change 38.10%
calculating round 3
report round 3 minimize init 78 final 72 change 7.69%
scorereport init 1668 final 72 change 95.68%
chromosomes_order = Scaffold_3,scaffold15,scaffold61,scaffold176,scaffold10,scaffold3,scaffold144,scaffold111,scaffold42,scaffold23,Scaffold_19,scaffold14,scaffold74,scaffold292,scaffold25,Scaffold_28,scaffold39,scaffold36,Scaffold_80,scaffold66,scaffold6,scaffold22,scaffold20,Scaffold_73,scaffold51,scaffold26,scaffold24,scaffold73,scaffold104,scaffold29,Scaffold_60,scaffold1,Scaffold_50,scaffold28,scaffold43,Scaffold_44,scaffold38,scaffold12,scaffold2,scaffold31,Scaffold_25,scaffold89,scaffold72,scaffold18,Scaffold_103,scaffold326,scaffold30,scaffold13,scaffold16,scaffold8

#circos.conf
#############################################################################
karyotype = ./karyotype.conf
chromosomes_units = 100000
<<include ideogram.conf>>
<<include ticks.conf>>
<<include bands.conf>>

<links>
<link>
 file=SyntenicRibbons.conf
 radius = 0.94r
 bezier_radius = 0.1r
 thickness = 1
 ribbon = yes
</link>
</links>



<image>
<<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
<<include ./housekeeping.conf>>
chromosomes_order = Scaffold_3,scaffold15,scaffold61,scaffold176,scaffold10,scaffold3,scaffold144,scaffold111,scaffold42,scaffold23,Scaffold_19,scaffold14,scaffold74,scaffold292,scaffold25,Scaffold_28,scaffold39,scaffold36,Scaffold_80,scaffold66,scaffold6,scaffold22,scaffold20,Scaffold_73,scaffold51,scaffold26,scaffold24,scaffold73,scaffold104,scaffold29,Scaffold_60,scaffold1,Scaffold_50,scaffold28,scaffold43,Scaffold_44,scaffold38,scaffold12,scaffold2,scaffold31,Scaffold_25,scaffold89,scaffold72,scaffold18,Scaffold_103,scaffold326,scaffold30,scaffold13,scaffold16,scaffold8


#############################################################################

ticks.conf
###############################################################################
show_ticks          = yes
show_tick_labels    = yes
<ticks>
radius           = 1r
color            = black
thickness        = 10p
multiplier       = 1e-5
format           = %d
<tick>
spacing        = 10u
size           = 25p
show_label     = yes
label_size     = 25p
label_offset   = 10p
format         = %d
</tick>

</ticks>

###############################################################################

bands.conf
###############################################################################
<bands>
show_bands = yes
fill_bands = yes
band_transparency = 4
</bands>
###############################################################################

ideogram.conf
###############################################################################

<ideogram>
<spacing>
  default = 0.006r
  break   = 30u
  axis_break_at_edge = yes
  axis_break         = yes
  axis_break_style   = 2
  <break_style 1>
        stroke_color     = black
        thickness        = 0.45r
        stroke_thickness = 2p
  </break>
  <break_style 2>
        stroke_color     = black
        stroke_thickness = 5p
        thickness        = 4r
  </break>
</spacing>
radius           = 0.84r
thickness        = 80p
fill             = yes
stroke_color     = white
stroke_thickness = 4p
fill_color       = black
show_label       = yes
label_font       = bold
label_size       = 16
label_parallel   = no

label_radius = dims(ideogram,radius_outer) + 0.06r
</ideogram>
###############################################################################

cp /work/GIF/software/programs/circos/0.69-4/etc/housekeeping.conf .
```
