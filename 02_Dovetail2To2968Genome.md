# How did Dovetail scaffold the second scn dovetail assembly?

```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/18_2692MummerDovetail2

module load GIF/mummer/4.0.0beta2
mummer -mum -l 1000 -b -n -qthreads 16 nematode_sp_19Jul2018_IbtP1.fasta 2968Genome.fasta >Dove2to2968Mummer.out

#How many scaffolds represented from 2968 assembly
less Dove2to2968Mummer.out |awk -v source=0 '{if(substr($1,1,1)==">") {source=$0;} else {print source,$0}}' |awk '{if($3=="Reverse") {print $2,$6,$6+$7"\t"$4,$5,$5+$7} else {print $2,$5,$5+$6"\t"$3,$4,$4+$6}}' |sort -k1,1V -k2,3n |awk '{print $1}' |sort |uniq |wc
   2450    2450   19495

#How many scaffolds represented in the second dovetail assembly
less Dove2to2968Mummer.out |awk -v source=0 '{if(substr($1,1,1)==">") {source=$0;} else {print source,$0}}' |awk '{if($3=="Reverse") {print $2,$6,$6+$7"\t"$4,$5,$5+$7} else {print $2,$5,$5+$6"\t"$3,$4,$4+$6}}' |sort -k1,1V -k2,3n |awk '{print $4}' |sort |uniq |wc
       587     587   13982

less Dove2to2968Mummer.out |awk -v source=0 '{if(substr($1,1,1)==">") {source=$0;} else {print source,$0}}' |awk '{if($3=="Reverse") {print $2,$6,$6+$7"\t"$4,$5,$5+$7} else {print $2,$5,$5+$6"\t"$3,$4,$4+$6}}' |sort -k1,1V -k2,3n >SyntenicRibbons.txt


 bioawk -c fastx '{print $name,length($seq)}' 2968Genome.fasta  |awk '{print "chr","-",$1,$1,"0",$2,"blue"}'  >scn1Karyotype.conf
bioawk -c fastx '{print $name,length($seq)}' nematode_sp_19Jul2018_IbtP1.fasta |sed 's/>/>SCN2/g' |sed 's/\.3//g' |sed 's/;/\t/g' |awk '{print "chr","-",$1,$1,"0",$3,"green"}' >scn2Karyotype.conf

cat <(sort scn1Karyotype.conf |uniq) <(sort scn2Karyotype.conf |uniq) >karyotype.conf
cat <(awk '{print $1}' SyntenicRibbons.txt |sort|uniq) <(awk '{print $4}' SyntenicRibbons.txt |sort|uniq) |grep -f - karyotype.conf >karyotypeRevised.conf
```
