#  This new genome turned out from Dovetail looks very different from our previous assembly.   This scaffolded assembly originated from the 2692 contigs that were generated directly out of falcon.


### Map the maker genes to the new genome
```
#rename the proteins so that iadhore doesnt get confused as to which protein belongs to which genome.
sed 's/>/>NEW/g'  ../09_Maker/01_maker/DovetailSCNMaker4.all.maker.transcripts.fasta|awk '{print $1}'  >NewTranscripts.fasta

#map these proteins to the genome
sh runGmap.sh SCN3 /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/17_Synteny2DovetailGenomes nematode_sp_19Jul2018_IbtP1.fasta NewTranscripts.fasta

```

### run Orthofinder
```
#create folder with proteins for all species, just primary isoform
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/17_Synteny2DovetailGenomes/Orthofinder
less ../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.maker.proteins.fasta |awk 'substr ($1,15,17) =="RA" {print $1}' |grep ">" |sed 's/>//g'|cdbyank ../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.maker.proteins.fasta.cidx |awk '{print $1}' |sed 's/>/>New/g' >NewPrimaryIsoform.fasta
less ../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.maker.proteins.fasta |awk 'substr ($1,15,17) =="RA" {print $1}' |grep ">" |sed 's/>//g'|cdbyank ../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.maker.proteins.fasta.cidx |awk '{print $1}' >PrimaryIsoform.fasta


module load orthofinder
orthofinder -t 12 -f Orthofinder


#Hmm. some of the same exact proteins are not found in the same family, need to fix that.
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/17_Synteny2DovetailGenomes/Orthofinder/Results_Jul20_1
less Orthogroups.txt |grep -v "New" |awk '{print $2}' |grep -f - Orthogroups.txt |sed 's/HetGly\./HetGly\.\t/1' |sed 's/-RA/\t-RA/1' |sort -k3,3nr |sed 's/HetGly\.\t/HetGly\./1' |sed 's/\t-RA/-RA/1' |awk '{print $1}' |sed '$!N;s/\n/ /' |awk '{print "sed -i s/"$2"/"$1"/g"}' |sed "s/s/'s/2" |sed "s/g/g' Orthogroups.txt /g" >FixOrthogroups.sh
sh FixOrthogroups.sh

#create ortholog family list for iadhore
less Orthogroups.txt |tr " " "\n" |awk -v Family=0 '{if(substr($1,1,2)=="OG") {Family=$1} else {print $0,Family}}' |sed 's/://g' |cat - <(dos2unix Or
thogroups_UnassignedGenes.csv|awk 'NR>1 {print $2,$1}') |tr " " "\t" >Orthologues.list


```

### Create Files for iadhore
```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/17_Synteny2DovetailGenomes/01_Synteny
less  ../../09_Maker/01_maker/DovetailSCNMaker4.all.NOFASTA.gff |awk '$3=="mRNA"' |sed 's/;/\t/g' |sed 's/ID=//g' |awk 'substr($9,14,15)=="RA"' |wc
  22856  352588 4918458
grep -c ">" ../Orthofinder/PrimaryIsoform.fasta
  22856

mkdir query
cd query/
less  ../../../09_Maker/01_maker/DovetailSCNMaker4.all.NOFASTA.gff |awk '$3=="mRNA"' |sed 's/;/\t/g' |sed 's/ID=//g' |awk 'substr($9,14,15)=="RA" {print $
  9$7,$1}' |sort|uniq|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini

mkdir subject
cd subject/
less ../SCN3.NewTranscripts.gff3 |awk '$3=="mRNA" ' |sed 's/;/\t/g' |sed 's/ID=//g' |sed 's/^/SCN2/g' |awk  '{print $9,$7,$1}' |grep "RA\.mrna1" |sed 's/\.mrna1//g' |awk '{print $1$2,$3}' |sort|uniq|sed 's/NEW/New/g' |awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/ .*//g' *.lst
ls *lst >input.txt
paste <(cut -f 1,2 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
cat query/query.ini subject/subject.ini |tr "\t" " " >iadhore.ini
vi iadhore.ini
i-adhore iadhore.ini  
```
### Create circos plot
```
mkdir 02_circos
cd 02_circos/
less DovetailSCNMaker4.all.NOFASTA.gff |awk '$3=="mRNA"' |sed 's/;/\t/g' |sed 's/ID=//g' |awk 'substr($9,14,16)=="RA"' >SCN1Grepmod.gff

#having some issues with col3 and col4 being in the wrong order.
#less SCN3.NewTranscripts.gff3 |awk '$3=="mRNA"' |sed 's/;/\t/g' |sed 's/ID=//g' |sed 's/\.mrna1//g' |awk 'substr($10,17,19)=="RA"' |sed 's/NEW/New/g' |sed 's/^/SCN2/g' |cut -f 1,3- >SCN2Grepmod.gff
less SCN3.NewTranscripts.gff3 |awk '$3=="mRNA" ' |sed 's/;/\t/g' |sed 's/ID=//g' |sed 's/^/SCN2/g' |grep "RA\.mrna1" |sed 's/\.mrna1//g' |cut -f 1,3- >SCN2Grepmod.gff


ln -s ../nematode_sp_19Jul2018_IbtP1.fasta
ln -s ../../12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff
ln -s ../../12_MakerGenesOrthofinder/DovetailSCNMaker4.genome.fasta
ln -s ../Orthofinder/PrimaryIsoform.fasta
ln -s ../Orthofinder/NewPrimaryIsoform.fasta
ln -s ../SCN3.NewTranscripts.gff3
ln -s ../01_Synteny/output/segments.txt

### The + and - did not matter here, so had to remove and just print.
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn2") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line SCN2Grepmod.gff; done |awk '{print $4}' >Col3
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn2") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line SCN2Grepmod.gff; done |awk '{print $5}' >Col4
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn2") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -w $line SCN1Grepmod.gff; done |awk '{print $4}' >Col7
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn2") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -w $line SCN1Grepmod.gff; done |awk '{print $5}' >Col8

less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '$3!=$10 {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '{if($5=="scn2") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '{print $2,$6}' |awk 'NR>1' |paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf

#cant figure out why col3 and col2 are not always set in the correct direction, so fixed it later
awk '{if($3>$2) {print $0} else {print $1,$3,$2,$4,$5,$6}}' SyntenicRibbons.conf >RevisedSyntenicRibbons.conf
###  Interestingly, there are only 443 multiplicons here, when there were ~781 multiplicons total reported in segments.txt

bioawk -c fastx '{print $name,length($seq)}' DovetailSCNMaker4.genome.fasta  |awk '{print "chr","-",$1,$1,"0",$2,"blue"}'  >scn1Karyotype.conf
bioawk -c fastx '{print $name,length($seq)}' nematode_sp_19Jul2018_IbtP1.fasta |sed 's/>/>SCN2/g' |sed 's/\.3//g' |sed 's/;/\t/g' |cut -f 1,3 |awk '{print "chr","-","SCN2"$1,"SCN2"$1,"0",$2,"green"}' >scn2Karyotype.conf

#How large are the multiplicons between the two genomes?
#dovetail genome 1
less SyntenicRibbons.conf |awk '{print $6-$5}' |summary.sh
Total:  106,843,754
Count:  443
Mean:   241,182
Median: 57,581
Min:    3,313
Max:    5,633,211


#dovetail genome 2
less RevisedSyntenicRibbons.conf |awk '{print $3-$2}' |summary.sh
Total:  186,426,551
Count:  443
Mean:   420,827
Median: 48,333
Min:    191
Max:    18,461,083

#this is pretty weird, since the genome is only 160MB.

```
### Investigate missing multiplicons
```
# multiplicons that are internal to dovetail genome 1
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '{print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn2") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '$1=="scn1" && $5=="scn1"' |wc
    225    1800   22147
    less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '{print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn1") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '$1=="scn1" && $5=="scn1"' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line SCN1Grepmod.gff; done |awk '{print $4}' >Col3
    less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '{print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn1") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '$1=="scn1" && $5=="scn1"' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line SCN1Grepmod.gff; done |awk '{print $5}' >Col4
    less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '{print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn1") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '$1=="scn1" && $5=="scn1"' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -w $line SCN1Grepmod.gff; done |awk '{print $4}' >Col7
    less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '{print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn1") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '$1=="scn1" && $5=="scn1"' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -w $line SCN1Grepmod.gff; done |awk '{print $5}' >Col8

    less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk ' {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '{if($5=="scn1") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '$1=="scn1" && $5=="scn1"' |awk '{print $2,$6}' |awk 'NR>1' |paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SyntenicRibbons.conf


#multiplicons that are internal to dovetail genome 2
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '{print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn2") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '$1=="scn2" && $5=="scn2"' |wc
   113     904   13537

less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '{print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn2") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '$1=="scn2" && $5=="scn2"' |awk '{print $3}' |sed '/^$/d' |while read line; do grep -w $line SCN2Grepmod.gff; done |awk '{print $4}' >Col3
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '{print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn1") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '$1=="scn2" && $5=="scn2"' |awk '{print $4}' |sed '/^$/d' |while read line; do grep -w $line SCN2Grepmod.gff; done |awk '{print $5}' >Col4
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '{print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn1") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '$1=="scn2" && $5=="scn2"' |awk '{print $7}' |sed '/^$/d' |while read line; do grep -w $line SCN2Grepmod.gff; done |awk '{print $4}' >Col7
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk '{print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}' |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" " " |tr "#" "\n" |awk '{if($5=="scn1") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '$1=="scn2" && $5=="scn2"' |awk '{print $8}' |sed '/^$/d' |while read line; do grep -w $line SCN2Grepmod.gff; done |awk '{print $5}' >Col8
less segments.txt |awk 'NR>1' |sed 'N;s/\n/ /' |awk ' {print $1,$2,$3,$4,$5,$6,$7"\n"$8,$9,$10,$11,$12,$13,$14}'  |awk '{if(NR%2) {print "#"$3,$4,$5,$6}else {print $3,$4,$5,$6}}' |tr "\n" "\t" |tr "#" "\n" |awk '{if($5=="scn2") {print $5,$6,$7,$8,$1,$2,$3,$4} else {print $1,$2,$3,$4,$5,$6,$7,$8}}' |awk '$1=="scn2" && $5=="scn2"' |awk '{print $2,$6}' |paste - Col3 Col4 Col7 Col8 |awk '{print $1,$3,$4,$2,$5,$6}' >SCN2SyntenicRibbons.conf
```
