# Need to run synteny to determine relative contiguity


### X12 genome assembly
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/01_X12

ln -s ../../28_X12GenomeComparison/pasa2.longest.filter.pep
ln -s ../../28_X12GenomeComparison/pasa2.longest.filter.gff3
ln -s ../../10_RepeatModeler/SCNgenome.fasta
ln -s ../../25_AnnotateGenes/mikado.loci.gff3
ln -s ../../25_AnnotateGenes/mikado_proteinsFixed.fasta

ml blast-plus
makeblastdb -in mikado_proteinsFixed.fasta -dbtype prot -out mikado_proteinsFixed

echo "ml blast-plus; blastp -db mikado_proteinsFixed -query pasa2.longest.filter.pep -outfmt 6 -num_threads 16 -out pasa2mikado.blastout" >blastp.sh

#filtering the blast
awk '$4>70' ../../Revised.pasa2mikado.blastout >xyz.blast

awk '$3=="mRNA" {print $1,$4,$5,$9}' mikado.loci.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print "10"$1,$4,$2,$3}' |tr " " "\t" >mikadoMCFormat.gff


#get scaffolds larger than 100kb
bioawk -c fastx '{print $name,length($seq)}' ../../../../28_X12GenomeComparison/SCN_genome.fa |awk '$2>100000 {print $1}' |sort|uniq|tr "\n" "," |less


java -classpath /opt/rit/app/mcscanx/20170403/bin/  dot_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o dotplotX12.png

java -classpath /opt/rit/app/mcscanx/20170403/bin/  bar_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o BarX12.png




control_file.ctl
##################################################################
3000
3000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9
12chr1,12chr2,12chr3,12chr4,12chr5,12chr6,12chr7,12chr8,12chr9,12chr_4_pilon,12chr_18_pilon,12chr_21_pilon,12chr_26_pilon,12chr_27_pilon,12chr_30_pilon,12chr_154_pilon,12chr_236_pilon,12chr_244_pilon,12chr_274_pilon,12chr_287_pilon,12chr_288_pilon,12chr_292_pilon,12chr_303_pilon,12chr_323_pilon,12chr_329_pilon,12chr_338_pilon,12chr_402_pilon,12chr_407_pilon,12chr_410_pilon,12chr_415_pilon,12chr_443_pilon,12chr_481_pilon,12chr_496_pilon,12chr_508_pilon,12chr_528_pilon,12chr_542_pilon,12chr_567_pilon,12chr_576_pilon,12chr_638_pilon,12chr_644_pilon,12chr_651_pilon,12chr_675_pilon,12chr_725_pilon,12chr_726_pilon
##################################################################


 java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c Circlecontrol_file.ctl -o CircleX12.png
Circlecontrol_file.ctl
##############################################################
3000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9,12chr1,12chr2,12chr3,12chr4,12chr5,12chr6,12chr7,12chr8,12chr9,12chr_4_pilon,12chr_18_pilon,12chr_21_pilon,12chr_26_pilon,12chr_27_pilon,12chr_30_pilon,12chr_154_pilon,12chr_236_pilon,12chr_244_pilon,12chr_274_pilon,12chr_287_pilon,12chr_288_pilon,12chr_292_pilon,12chr_303_pilon,12chr_323_pilon,12chr_329_pilon,12chr_338_pilon,12chr_402_pilon,12chr_407_pilon,12chr_410_pilon,12chr_415_pilon,12chr_443_pilon,12chr_481_pilon,12chr_496_pilon,12chr_508_pilon,12chr_528_pilon,12chr_542_pilon,12chr_567_pilon,12chr_576_pilon,12chr_638_pilon,12chr_644_pilon,12chr_651_pilon,12chr_675_pilon,12chr_725_pilon,12chr_726_pilon
##############################################################  

#number of syntenic blocks
less xyz.collinearity |awk 'substr($1,1,3)=="##" && NR>3' |wc
    367    2936   27970


```

### SCN genome assembly v1
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/02_738Assembly/01_mcscanx

ln -s /work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/1_genomeNgff/genome738sl.polished.mitoFixed.fa
ln -s /work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/1_genomeNgff/fixed.augustus.gff3
ln -s /work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/1_genomeNgff/augustus.aa
ln -s ../01_X12/mikado_proteinsFixed.fasta
ln -s ../01_X12/mikado_proteinsFixed.phr
ln -s ../01_X12/mikado_proteinsFixed.pin
ln -s ../01_X12/mikado_proteinsFixed.psq
ln -s ../01_X12/mikadoMCFormat.gff
cp ../01_X12/blastp.sh .

ml blast-plus; blastp -db mikado_proteinsFixed -query augustus.aa -outfmt 6 -num_threads 16 -out 7382mikado.blastout
less ../7382mikado.blastout |awk '$3>70' >xyz.blast
sed -i 's/t//g' xyz.blast

#Create gff
awk '$3=="mRNA" {print $1,$4,$5,$9}' fixed.augustus.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print "V1"$1,$4,$2,$3}' |tr " " "\t" >738MCFormat.gff
cat ../mikadoMCFormat.gff ../738MCFormat.gff >xyz.gff
sed -i 's/Hetgly.T/Hetgly.G/g' xyz.gff
MCScanX xyz

java -classpath /opt/rit/app/mcscanx/20170403/bin/  dot_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o dotplotV1.png
java -classpath /opt/rit/app/mcscanx/20170403/bin/  bar_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o BarV1.png


control_file.ctl
#################################################################################
3000
3000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9
V1scaffold_1,V1scaffold_2,V1scaffold_3,V1scaffold_4,V1scaffold_5,V1scaffold_6,V1scaffold_7,V1scaffold_8,V1scaffold_9,V1scaffold_10,V1scaffold_11,V1scaffold_12,V1scaffold_13,V1scaffold_14,V1scaffold_15,V1scaffold_16,V1scaffold_17,V1scaffold_18,V1scaffold_19,V1scaffold_20,V1scaffold_21,V1scaffold_22,V1scaffold_23,V1scaffold_24,V1scaffold_25,V1scaffold_26,V1scaffold_27,V1scaffold_28,V1scaffold_29,V1scaffold_30,V1scaffold_31,V1scaffold_32,V1scaffold_33,V1scaffold_34,V1scaffold_35,V1scaffold_36,V1scaffold_37,V1scaffold_38,V1scaffold_39,V1scaffold_40,V1scaffold_41,V1scaffold_42,V1scaffold_43,V1scaffold_44,V1scaffold_45,V1scaffold_46,V1scaffold_47,V1scaffold_48,V1scaffold_49,V1scaffold_50,V1scaffold_51
#################################################################################

java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c Circlecontrol_file.ctl -o CircleV1.png
Circlecontrol_file.ctl
#################################################################################
3000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9,V1scaffold_1,V1scaffold_2,V1scaffold_3,V1scaffold_4,V1scaffold_5,V1scaffold_6,V1scaffold_7,V1scaffold_8,V1scaffold_9,V1scaffold_10,V1scaffold_11,V1scaffold_12,V1scaffold_13,V1scaffold_14,V1scaffold_15,V1scaffold_16,V1scaffold_17,V1scaffold_18,V1scaffold_19,V1scaffold_20,V1scaffold_21,V1scaffold_22,V1scaffold_23,V1scaffold_24,V1scaffold_25,V1scaffold_26,V1scaffold_27,V1scaffold_28,V1scaffold_29,V1scaffold_30,V1scaffold_31,V1scaffold_32,V1scaffold_33,V1scaffold_34,V1scaffold_35,V1scaffold_36,V1scaffold_37,V1scaffold_38,V1scaffold_39,V1scaffold_40,V1scaffold_41,V1scaffold_42,V1scaffold_43,V1scaffold_44,V1scaffold_45,V1scaffold_46,V1scaffold_47,V1scaffold_48,V1scaffold_49,V1scaffold_50,V1scaffold_51
#################################################################################

#number of syntenic blocks
less xyz.collinearity |awk 'substr($1,1,3)=="##" && NR>3' |wc
   1727   13816  144318



```


### G. pallida
```

ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/globodera_pallida.PRJEB123.WBPS10.protein.fa
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/globodera_pallida.PRJEB123.WBPS10.genomic.fa
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/globodera_pallida.PRJEB123.WBPS10.annotations.gff3
ln -s ../01_X12/mikado_proteinsFixed.fasta
ln -s ../01_X12/mikado_proteinsFixed.phr
ln -s ../01_X12/mikado_proteinsFixed.pin
ln -s ../01_X12/mikado_proteinsFixed.psq
ln -s ../01_X12/mikadoMCFormat.gff
cp ../01_X12/blastp.sh .
cp ../02_738Assembly/01_mcscanx/control_file.ctl .

#filtering the blast
less Gpal2mikado.blastout |awk '$3>70'  >xyz.blast


awk '$3=="gene" {print $1,$4,$5,$9}' globodera_pallida.PRJEB123.WBPS10.annotations.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print "10"$1,$4,$2,$3}' |tr " " "\t" |sed 's/gene://g' >GpalMCFormat.gff

bioawk -c fastx '{print $name,length($seq)}' globodera_pallida.PRJEB123.WBPS10.genomic.fa |sort -k2,2nr |awk 'NR<51 {print "Gp"$1}' |tr "\n" "," |less

cat mikadoMCFormat.gff GpalMCFormat.gff >xyz.gff

MCScanX xyz

java -classpath /opt/rit/app/mcscanx/20170403/bin/  dot_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o dotplotGpal.png

java -classpath /opt/rit/app/mcscanx/20170403/bin/  bar_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o BarGpal.png

control_file.ctl
##################################################################
3000
3000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9
Gppathogens_Gpal_scaffold_1,Gppathogens_Gpal_scaffold_2,Gppathogens_Gpal_scaffold_3,Gppathogens_Gpal_scaffold_5,Gppathogens_Gpal_scaffold_4,Gppathogens_Gpal_scaffold_6,Gppathogens_Gpal_scaffold_7,Gppathogens_Gpal_scaffold_8,Gppathogens_Gpal_scaffold_9,Gppathogens_Gpal_scaffold_10,Gppathogens_Gpal_scaffold_11,Gppathogens_Gpal_scaffold_12,Gppathogens_Gpal_scaffold_14,Gppathogens_Gpal_scaffold_13,Gppathogens_Gpal_scaffold_15,Gppathogens_Gpal_scaffold_16,Gppathogens_Gpal_scaffold_17,Gppathogens_Gpal_scaffold_18,Gppathogens_Gpal_scaffold_19,Gppathogens_Gpal_scaffold_20,Gppathogens_Gpal_scaffold_21,Gppathogens_Gpal_scaffold_22,Gppathogens_Gpal_scaffold_23,Gppathogens_Gpal_scaffold_25,Gppathogens_Gpal_scaffold_26,Gppathogens_Gpal_scaffold_27,Gppathogens_Gpal_scaffold_28,Gppathogens_Gpal_scaffold_30,Gppathogens_Gpal_scaffold_24,Gppathogens_Gpal_scaffold_29,Gppathogens_Gpal_scaffold_32,Gppathogens_Gpal_scaffold_31,Gppathogens_Gpal_scaffold_33,Gppathogens_Gpal_scaffold_34,Gppathogens_Gpal_scaffold_35,Gppathogens_Gpal_scaffold_36,Gppathogens_Gpal_scaffold_38,Gppathogens_Gpal_scaffold_37,Gppathogens_Gpal_scaffold_39,Gppathogens_Gpal_scaffold_41,Gppathogens_Gpal_scaffold_42,Gppathogens_Gpal_scaffold_40,Gppathogens_Gpal_scaffold_43,Gppathogens_Gpal_scaffold_44,Gppathogens_Gpal_scaffold_45,Gppathogens_Gpal_scaffold_46,Gppathogens_Gpal_scaffold_48,Gppathogens_Gpal_scaffold_47,Gppathogens_Gpal_scaffold_51,Gppathogens_Gpal_scaffold_49
##################################################################


 java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c Circlecontrol_file.ctl -o CircleHsc.png
Circlecontrol_file.ctl
##############################################################
3000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9,Gppathogens_Gpal_scaffold_1,Gppathogens_Gpal_scaffold_2,Gppathogens_Gpal_scaffold_3,Gppathogens_Gpal_scaffold_5,Gppathogens_Gpal_scaffold_4,Gppathogens_Gpal_scaffold_6,Gppathogens_Gpal_scaffold_7,Gppathogens_Gpal_scaffold_8,Gppathogens_Gpal_scaffold_9,Gppathogens_Gpal_scaffold_10,Gppathogens_Gpal_scaffold_11,Gppathogens_Gpal_scaffold_12,Gppathogens_Gpal_scaffold_14,Gppathogens_Gpal_scaffold_13,Gppathogens_Gpal_scaffold_15,Gppathogens_Gpal_scaffold_16,Gppathogens_Gpal_scaffold_17,Gppathogens_Gpal_scaffold_18,Gppathogens_Gpal_scaffold_19,Gppathogens_Gpal_scaffold_20,Gppathogens_Gpal_scaffold_21,Gppathogens_Gpal_scaffold_22,Gppathogens_Gpal_scaffold_23,Gppathogens_Gpal_scaffold_25,Gppathogens_Gpal_scaffold_26,Gppathogens_Gpal_scaffold_27,Gppathogens_Gpal_scaffold_28,Gppathogens_Gpal_scaffold_30,Gppathogens_Gpal_scaffold_24,Gppathogens_Gpal_scaffold_29,Gppathogens_Gpal_scaffold_32,Gppathogens_Gpal_scaffold_31,Gppathogens_Gpal_scaffold_33,Gppathogens_Gpal_scaffold_34,Gppathogens_Gpal_scaffold_35,Gppathogens_Gpal_scaffold_36,Gppathogens_Gpal_scaffold_38,Gppathogens_Gpal_scaffold_37,Gppathogens_Gpal_scaffold_39,Gppathogens_Gpal_scaffold_41,Gppathogens_Gpal_scaffold_42,Gppathogens_Gpal_scaffold_40,Gppathogens_Gpal_scaffold_43,Gppathogens_Gpal_scaffold_44,Gppathogens_Gpal_scaffold_45,Gppathogens_Gpal_scaffold_46,Gppathogens_Gpal_scaffold_48,Gppathogens_Gpal_scaffold_47,Gppathogens_Gpal_scaffold_51,Gppathogens_Gpal_scaffold_49
##############################################################  

#number of syntenic blocks
less xyz.collinearity |awk 'substr($1,1,3)=="##" && NR>3' |wc
   607    4856   59461


```

### G. ellingtonae
```
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/G.ellingtonae.genomic.fa
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/G.ellingtonae.gff3
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/G.ellingtonae_protein.fasta
ln -s ../01_X12/mikado_proteinsFixed.fasta
ln -s ../01_X12/mikado_proteinsFixed.phr
ln -s ../01_X12/mikado_proteinsFixed.pin
ln -s ../01_X12/mikado_proteinsFixed.psq
ln -s ../01_X12/mikadoMCFormat.gff
cp ../01_X12/blastp.sh .

less Gell2mikado.blastout |awk '$3>70'  >xyz.blast


awk '$3=="mRNA" {print $1,$4,$5,$9}' G.ellingtonae.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print "GE"$1,$4,$2,$3}' |tr " " "\t" >GellMCFormat.gff

bioawk -c fastx '{print $name,length($seq)}' G.ellingtonae.genomic.fa |sort -k2,2nr |awk 'NR<51 {print "GE"$1}' |tr "\n" "," |less

cat mikadoMCFormat.gff GellMCFormat.gff >xyz.gff

MCScanX xyz

java -classpath /opt/rit/app/mcscanx/20170403/bin/  dot_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o dotplotGell.png

java -classpath /opt/rit/app/mcscanx/20170403/bin/  bar_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o BarGell.png

control_file.ctl
##################################################################
4000
4000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9
GEMEIZ01000001.1,GEMEIZ01000002.1,GEMEIZ01000003.1,GEMEIZ01000004.1,GEMEIZ01000005.1,GEMEIZ01000006.1,GEMEIZ01000007.1,GEMEIZ01000009.1,GEMEIZ01000008.1,GEMEIZ01000010.1,GEMEIZ01000011.1,GEMEIZ01000013.1,GEMEIZ01000014.1,GEMEIZ01000012.1,GEMEIZ01000016.1,GEMEIZ01000015.1,GEMEIZ01000017.1,GEMEIZ01000020.1,GEMEIZ01000019.1,GEMEIZ01000021.1,GEMEIZ01000018.1,GEMEIZ01000023.1,GEMEIZ01000024.1,GEMEIZ01000022.1,GEMEIZ01000029.1,GEMEIZ01000028.1,GEMEIZ01000027.1,GEMEIZ01000036.1,GEMEIZ01000035.1,GEMEIZ01000030.1,GEMEIZ01000034.1,GEMEIZ01000025.1,GEMEIZ01000031.1,GEMEIZ01000032.1,GEMEIZ01000026.1,GEMEIZ01000033.1,GEMEIZ01000041.1,GEMEIZ01000040.1,GEMEIZ01000037.1,GEMEIZ01000048.1,GEMEIZ01000046.1,GEMEIZ01000043.1,GEMEIZ01000039.1,GEMEIZ01000044.1,GEMEIZ01000042.1,GEMEIZ01000053.1,GEMEIZ01000051.1,GEMEIZ01000052.1,GEMEIZ01000038.1,GEMEIZ01000049.1
##################################################################


 java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c Circlecontrol_file.ctl -o CircleGell.png
Circlecontrol_file.ctl
##############################################################
3000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9,GEMEIZ01000001.1,GEMEIZ01000002.1,GEMEIZ01000003.1,GEMEIZ01000004.1,GEMEIZ01000005.1,GEMEIZ01000006.1,GEMEIZ01000007.1,GEMEIZ01000009.1,GEMEIZ01000008.1,GEMEIZ01000010.1,GEMEIZ01000011.1,GEMEIZ01000013.1,GEMEIZ01000014.1,GEMEIZ01000012.1,GEMEIZ01000016.1,GEMEIZ01000015.1,GEMEIZ01000017.1,GEMEIZ01000020.1,GEMEIZ01000019.1,GEMEIZ01000021.1,GEMEIZ01000018.1,GEMEIZ01000023.1,GEMEIZ01000024.1,GEMEIZ01000022.1,GEMEIZ01000029.1,GEMEIZ01000028.1,GEMEIZ01000027.1,GEMEIZ01000036.1,GEMEIZ01000035.1,GEMEIZ01000030.1,GEMEIZ01000034.1,GEMEIZ01000025.1,GEMEIZ01000031.1,GEMEIZ01000032.1,GEMEIZ01000026.1,GEMEIZ01000033.1,GEMEIZ01000041.1,GEMEIZ01000040.1,GEMEIZ01000037.1,GEMEIZ01000048.1,GEMEIZ01000046.1,GEMEIZ01000043.1,GEMEIZ01000039.1,GEMEIZ01000044.1,GEMEIZ01000042.1,GEMEIZ01000053.1,GEMEIZ01000051.1,GEMEIZ01000052.1,GEMEIZ01000038.1,GEMEIZ01000049.1
##############################################################  

#number of syntenic blocks
less xyz.collinearity |awk 'substr($1,1,3)=="##" && NR>3' |wc
    367    2936   27970


```


### G. rostochiensis
```
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/globodera_rostochiensis.PRJEB13504.WBPS10.annotations.gff3
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/globodera_rostochiensis.PRJEB13504.WBPS10.genomic.fa
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/globodera_rostochiensis.PRJEB13504.WBPS10.protein.fa
ln -s ../01_X12/mikado_proteinsFixed.fasta
ln -s ../01_X12/mikado_proteinsFixed.phr
ln -s ../01_X12/mikado_proteinsFixed.pin
ln -s ../01_X12/mikado_proteinsFixed.psq
ln -s ../01_X12/mikadoMCFormat.gff
cp ../01_X12/blastp.sh .

less Gros2mikado.blastout |awk '$3>70'  >xyz.blast


awk '$3=="mRNA" {print $1,$4,$5,$9}' globodera_rostochiensis.PRJEB13504.WBPS10.annotations.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print "GR"$1,$4,$2,$3}' |tr " " "\t" |sed 's/transcript://g'  >GrosMCFormat.gff

bioawk -c fastx '{print $name,length($seq)}' globodera_rostochiensis.PRJEB13504.WBPS10.genomic.fa |sort -k2,2nr |awk 'NR<51 {print "GR"$1}' |tr "\n" "," |less

cat mikadoMCFormat.gff GrosMCFormat.gff >xyz.gff

MCScanX xyz

java -classpath /opt/rit/app/mcscanx/20170403/bin/  dot_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o dotplotGros.png

java -classpath /opt/rit/app/mcscanx/20170403/bin/  bar_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o BarGros.png

control_file.ctl
##################################################################
4000
4000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9
GRGROS_00001,GRGROS_00002,GRGROS_00003,GRGROS_00005,GRGROS_00004,GRGROS_00006,GRGROS_00007,GRGROS_00008,GRGROS_00010,GRGROS_00009,GRGROS_00011,GRGROS_00012,GRGROS_00013,GRGROS_00014,GRGROS_00016,GRGROS_00015,GRGROS_00018,GRGROS_00017,GRGROS_00019,GRGROS_00020,GRGROS_00023,GRGROS_00021,GRGROS_00022,GRGROS_00024,GRGROS_00025,GRGROS_00026,GRGROS_00027,GRGROS_00028,GRGROS_00031,GRGROS_00032,GRGROS_00029,GRGROS_00030,GRGROS_00033,GRGROS_00034,GRGROS_00037,GRGROS_00036,GRGROS_00035,GRGROS_00038,GRGROS_00040,GRGROS_00039,GRGROS_00042,GRGROS_00041,GRGROS_00043,GRGROS_00048,GRGROS_00049,GRGROS_00046,GRGROS_00047,GRGROS_00045,GRGROS_00044,GRGROS_00051
##################################################################


 java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c Circlecontrol_file.ctl -o CircleGros.png
Circlecontrol_file.ctl
##############################################################
3000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9,GRGROS_00001,GRGROS_00002,GRGROS_00003,GRGROS_00005,GRGROS_00004,GRGROS_00006,GRGROS_00007,GRGROS_00008,GRGROS_00010,GRGROS_00009,GRGROS_00011,GRGROS_00012,GRGROS_00013,GRGROS_00014,GRGROS_00016,GRGROS_00015,GRGROS_00018,GRGROS_00017,GRGROS_00019,GRGROS_00020,GRGROS_00023,GRGROS_00021,GRGROS_00022,GRGROS_00024,GRGROS_00025,GRGROS_00026,GRGROS_00027,GRGROS_00028,GRGROS_00031,GRGROS_00032,GRGROS_00029,GRGROS_00030,GRGROS_00033,GRGROS_00034,GRGROS_00037,GRGROS_00036,GRGROS_00035,GRGROS_00038,GRGROS_00040,GRGROS_00039,GRGROS_00042,GRGROS_00041,GRGROS_00043,GRGROS_00048,GRGROS_00049,GRGROS_00046,GRGROS_00047,GRGROS_00045,GRGROS_00044,GRGROS_00051
##############################################################  

#number of syntenic blocks
less xyz.collinearity |awk 'substr($1,1,3)=="##" && NR>3' |wc
    266    2128   21631



```

### M incognita
```
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/meloidogyne_incognita.PRJEA28837.WBPS10.annotations.gff3
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/meloidogyne_incognita.PRJEA28837.WBPS10.genomic.fa
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/meloidogyne_incognita.PRJEA28837.WBPS10.protein.fa
ln -s ../01_X12/mikado_proteinsFixed.fasta
ln -s ../01_X12/mikado_proteinsFixed.phr
ln -s ../01_X12/mikado_proteinsFixed.pin
ln -s ../01_X12/mikado_proteinsFixed.psq
ln -s ../01_X12/mikadoMCFormat.gff
cp ../01_X12/blastp.sh .

less Minc2mikado.blastout |awk '$3>70'  >xyz.blast


awk '$3=="mRNA" {print $1,$4,$5,$9}' meloidogyne_incognita.PRJEA28837.WBPS10.annotations.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print "MI"$1,$4,$2,$3}' |tr " " "\t" |sed 's/transcript://g'  >MincMCFormat.gff

cat mikadoMCFormat.gff MincMCFormat.gff >xyz.gff

MCScanX xyz

bioawk -c fastx '{print $name,length($seq)}' meloidogyne_incognita.PRJEA28837.WBPS10.genomic.fa |sort -k2,2nr |awk 'NR<51 {print "MI"$1}' |tr "\n" "," |less



vi control_file.ctl
##################################################################
4000
4000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9
MIMiV1ctg0,MIMiV1ctg2,MIMiV1ctg1,MIMiV1ctg3,MIMiV1ctg7,MIMiV1ctg5,MIMiV1ctg8,MIMiV1ctg4,MIMiV1ctg6,MIMiV1ctg10,MIMiV1ctg17,MIMiV1ctg9,MIMiV1ctg13,MIMiV1ctg12,MIMiV1ctg11,MIMiV1ctg16,MIMiV1ctg15,MIMiV1ctg14,MIMiV1ctg18,MIMiV1ctg19,MIMiV1ctg21,MIMiV1ctg20,MIMiV1ctg26,MIMiV1ctg25,MIMiV1ctg27,MIMiV1ctg34,MIMiV1ctg24,MIMiV1ctg36,MIMiV1ctg22,MIMiV1ctg23,MIMiV1ctg30,MIMiV1ctg28,MIMiV1ctg32,MIMiV1ctg29,MIMiV1ctg33,MIMiV1ctg31,MIMiV1ctg35,MIMiV1ctg37,MIMiV1ctg43,MIMiV1ctg39,MIMiV1ctg38,MIMiV1ctg40,MIMiV1ctg42,MIMiV1ctg45,MIMiV1ctg51,MIMiV1ctg41,MIMiV1ctg46,MIMiV1ctg44,MIMiV1ctg48,MIMiV1ctg54
##################################################################

java -classpath /opt/rit/app/mcscanx/20170403/bin/  dot_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o dotplotMinc.png

java -classpath /opt/rit/app/mcscanx/20170403/bin/  bar_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o BarMinc.png

vi Circlecontrol_file.ctl
##############################################################
3000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9,MIMiV1ctg0,MIMiV1ctg2,MIMiV1ctg1,MIMiV1ctg3,MIMiV1ctg7,MIMiV1ctg5,MIMiV1ctg8,MIMiV1ctg4,MIMiV1ctg6,MIMiV1ctg10,MIMiV1ctg17,MIMiV1ctg9,MIMiV1ctg13,MIMiV1ctg12,MIMiV1ctg11,MIMiV1ctg16,MIMiV1ctg15,MIMiV1ctg14,MIMiV1ctg18,MIMiV1ctg19,MIMiV1ctg21,MIMiV1ctg20,MIMiV1ctg26,MIMiV1ctg25,MIMiV1ctg27,MIMiV1ctg34,MIMiV1ctg24,MIMiV1ctg36,MIMiV1ctg22,MIMiV1ctg23,MIMiV1ctg30,MIMiV1ctg28,MIMiV1ctg32,MIMiV1ctg29,MIMiV1ctg33,MIMiV1ctg31,MIMiV1ctg35,MIMiV1ctg37,MIMiV1ctg43,MIMiV1ctg39,MIMiV1ctg38,MIMiV1ctg40,MIMiV1ctg42,MIMiV1ctg45,MIMiV1ctg51,MIMiV1ctg41,MIMiV1ctg46,MIMiV1ctg44,MIMiV1ctg48,MIMiV1ctg54
##############################################################  

 java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c Circlecontrol_file.ctl -o CircleGros.png

 #numbers of syntenic blocks
 less xyz.collinearity |awk 'substr($1,1,3)=="##" && NR>3' |wc
     6      48     470

```

### M hapla
```
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/meloidogyne_hapla.PRJNA29083.WBPS10.annotations.gff3
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/meloidogyne_hapla.PRJNA29083.WBPS10.genomic.fa
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/meloidogyne_hapla.PRJNA29083.WBPS10.protein.fa
ln -s ../01_X12/mikado_proteinsFixed.fasta
ln -s ../01_X12/mikado_proteinsFixed.phr
ln -s ../01_X12/mikado_proteinsFixed.pin
ln -s ../01_X12/mikado_proteinsFixed.psq
ln -s ../01_X12/mikadoMCFormat.gff
cp ../01_X12/blastp.sh .

less Mhap2mikado.blastout |awk '$3>70'  >xyz.blast


awk '$3=="mRNA" {print $1,$4,$5,$9}' meloidogyne_hapla.PRJNA29083.WBPS10.annotations.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print "MH"$1,$4,$2,$3}' |tr " " "\t" |sed 's/transcript://g'  >MhapMCFormat.gff

cat mikadoMCFormat.gff MhapMCFormat.gff >xyz.gff

MCScanX xyz

bioawk -c fastx '{print $name,length($seq)}' meloidogyne_hapla.PRJNA29083.WBPS10.genomic.fa |sort -k2,2nr |awk 'NR<51 {print "MH"$1}' |tr "\n" "," |less



vi control_file.ctl
##################################################################
4000
4000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9
MHMhA1_Contig88,MHMhA1_Contig123,MHMhA1_Contig118,MHMhA1_Contig147,MHMhA1_Contig0,MHMhA1_Contig127,MHMhA1_Contig121,MHMhA1_Contig138,MHMhA1_Contig168,MHMhA1_Contig194,MHMhA1_Contig309,MHMhA1_Contig42,MHMhA1_Contig134,MHMhA1_Contig20,MHMhA1_Contig210,MHMhA1_Contig163,MHMhA1_Contig342,MHMhA1_Contig353,MHMhA1_Contig53,MHMhA1_Contig108,MHMhA1_Contig7,MHMhA1_Contig354,MHMhA1_Contig130,MHMhA1_Contig373,MHMhA1_Contig70,MHMhA1_Contig261,MHMhA1_Contig280,MHMhA1_Contig271,MHMhA1_Contig91,MHMhA1_Contig120,MHMhA1_Contig253,MHMhA1_Contig320,MHMhA1_Contig213,MHMhA1_Contig34,MHMhA1_Contig484,MHMhA1_Contig310,MHMhA1_Contig338,MHMhA1_Contig334,MHMhA1_Contig185,MHMhA1_Contig239,MHMhA1_Contig107,MHMhA1_Contig30,MHMhA1_Contig468,MHMhA1_Contig357,MHMhA1_Contig325,MHMhA1_Contig515,MHMhA1_Contig96,MHMhA1_Contig104,MHMhA1_Contig65,MHMhA1_Contig50
##################################################################

java -classpath /opt/rit/app/mcscanx/20170403/bin/  dot_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o dotplotMhap.png

java -classpath /opt/rit/app/mcscanx/20170403/bin/  bar_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o BarMhap.png

vi Circlecontrol_file.ctl
##############################################################
3000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9,MHMhA1_Contig88,MHMhA1_Contig123,MHMhA1_Contig118,MHMhA1_Contig147,MHMhA1_Contig0,MHMhA1_Contig127,MHMhA1_Contig121,MHMhA1_Contig138,MHMhA1_Contig168,MHMhA1_Contig194,MHMhA1_Contig309,MHMhA1_Contig42,MHMhA1_Contig134,MHMhA1_Contig20,MHMhA1_Contig210,MHMhA1_Contig163,MHMhA1_Contig342,MHMhA1_Contig353,MHMhA1_Contig53,MHMhA1_Contig108,MHMhA1_Contig7,MHMhA1_Contig354,MHMhA1_Contig130,MHMhA1_Contig373,MHMhA1_Contig70,MHMhA1_Contig261,MHMhA1_Contig280,MHMhA1_Contig271,MHMhA1_Contig91,MHMhA1_Contig120,MHMhA1_Contig253,MHMhA1_Contig320,MHMhA1_Contig213,MHMhA1_Contig34,MHMhA1_Contig484,MHMhA1_Contig310,MHMhA1_Contig338,MHMhA1_Contig334,MHMhA1_Contig185,MHMhA1_Contig239,MHMhA1_Contig107,MHMhA1_Contig30,MHMhA1_Contig468,MHMhA1_Contig357,MHMhA1_Contig325,MHMhA1_Contig515,MHMhA1_Contig96,MHMhA1_Contig104,MHMhA1_Contig65,MHMhA1_Contig50

 java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c Circlecontrol_file.ctl -o CircleMhap.png

#number of syntenic blocks
less xyz.collinearity |awk 'substr($1,1,3)=="##" && NR>3' |wc
      8      64     665


```

### D dipsaci
```
ln -s /work/GIF/remkv6/Baum/06_D_dipsaci/10_SubtractContamination/augustus.pepRevised.fasta
ln -s /work/GIF/remkv6/Baum/06_D_dipsaci/10_SubtractContamination/augustusRevised.gff3
ln -s /work/GIF/remkv6/Baum/06_D_dipsaci/10_SubtractContamination/genome.fa
ln -s ../01_X12/mikado_proteinsFixed.fasta
ln -s ../01_X12/mikado_proteinsFixed.phr
ln -s ../01_X12/mikado_proteinsFixed.pin
ln -s ../01_X12/mikado_proteinsFixed.psq
ln -s ../01_X12/mikadoMCFormat.gff
cp ../01_X12/blastp.sh .

less Ddip2mikado.blastout |awk '$3>70'  >xyz.blast


awk '$3=="mRNA" {print $1,$4,$5,$9}' augustusRevised.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print "D1"$1,$4,$2,$3}' |tr " " "\t"   >DdipMCFormat.gff

cat mikadoMCFormat.gff DdipMCFormat.gff >xyz.gff

MCScanX xyz -k 300 -e 1e-03 -m 50

#only two syntenic blocks, not worth pursuing plots
less xyz.collinearity |awk 'substr($1,1,3)=="##" && NR>3' |wc
2
```

### D destructor
```
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/ditylenchus_destructor.PRJNA312427.WBPS10.annotations.gff3
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/ditylenchus_destructor.PRJNA312427.WBPS10.protein.fa
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/ditylenchus_destructor.PRJNA312427.WBPS10.genomic.fa
ln -s ../01_X12/mikado_proteinsFixed.fasta
ln -s ../01_X12/mikado_proteinsFixed.phr
ln -s ../01_X12/mikado_proteinsFixed.pin
ln -s ../01_X12/mikado_proteinsFixed.psq
ln -s ../01_X12/mikadoMCFormat.gff
cp ../01_X12/blastp.sh .

less Ddes2mikado.blastout |awk '$3>70'  >xyz.blast


awk '$3=="mRNA" {print $1,$4,$5,$9}' ditylenchus_destructor.PRJNA312427.WBPS10.annotations.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print "D2"$1,$4,$2,$3}' |tr " " "\t" |sed 's/transcript://g'   >DdesMCFormat.gff

cat mikadoMCFormat.gff DdesMCFormat.gff >xyz.gff

MCScanX xyz -k 300 -e 1e-03 -m 50

less xyz.collinearity |awk 'substr($1,1,2)=="##"' |awk '{print $7}' |sed 's/&/\t/g' |awk '{print $2}' |tr "\n" "," |less

vi control_file.ctl
##################################################################
4000
4000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9
D2scaffold1,D2scaffold14,D2scaffold26,D2scaffold3,D2scaffold42,D2scaffold6
##################################################################

java -classpath /opt/rit/app/mcscanx/20170403/bin/  dot_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o dotplotDdes.png

java -classpath /opt/rit/app/mcscanx/20170403/bin/  bar_plotter -g xyz.gff -s xyz.collinearity -c control_file.ctl -o BarDdes.png

vi Circlecontrol_file.ctl
##############################################################
3000
10Scaffold_1,10Scaffold_2,10Scaffold_3,10Scaffold_4,10Scaffold_5,10Scaffold_6,10Scaffold_7,10Scaffold_8,10Scaffold_9,D2scaffold1,D2scaffold14,D2scaffold26,D2scaffold3,D2scaffold42,D2scaffold6

 java -classpath /opt/rit/app/mcscanx/20170403/bin/  circle_plotter -g xyz.gff -s xyz.collinearity -c Circlecontrol_file.ctl -o CircleDdes.png


#number of syntenic blocks
less xyz.collinearity |awk 'substr($1,1,3)=="##" && NR>3' |wc
      9      72     717



```

### B xylophilus
```
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/bursaphelenchus_xylophilus.PRJEA64437.WBPS10.annotations.gff3
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/bursaphelenchus_xylophilus.PRJEA64437.WBPS10.genomic.fa
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/bursaphelenchus_xylophilus.PRJEA64437.WBPS10.protein.fa
ls
ln -s ../01_X12/mikado_proteinsFixed.fasta
ln -s ../01_X12/mikado_proteinsFixed.phr
ln -s ../01_X12/mikado_proteinsFixed.pin
ln -s ../01_X12/mikado_proteinsFixed.psq
ln -s ../01_X12/mikadoMCFormat.gff
cp ../01_X12/blastp.sh .

less Bxyl2mikado.blastout |awk '$3>70'  >xyz.blast


awk '$3=="mRNA" {print $1,$4,$5,$9}' bursaphelenchus_xylophilus.PRJEA64437.WBPS10.annotations.gff3 |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print "BX"$1,$4,$2,$3}' |tr " " "\t" |sed 's/transcript://g'   >BxylMCFormat.gff

cat mikadoMCFormat.gff BxylMCFormat.gff >xyz.gff

MCScanX xyz -k 300 -e 1e-03 -m 50

# zero syntenic blocks
````
