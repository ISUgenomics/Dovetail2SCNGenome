# This is a run of Mummer to determine the extent of internal duplication found in the SCN genome.

## Run mummer
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/50_tandemMummer
module load mummer
run-mummer3 genome738sl.polished.mitoFixed.fa genome738sl.polished.mitoFixed.fa self
show-coords -r seq_seq.delta > seq_seq.coords

less seq_seq.coords |awk 'NR>4' |awk ' $7>1000' |awk '$12!=$13' >diff.scaff.tand.coords
```

## Create a bed file that contains all overlapping coordinates.
```
less seq_seq.coords |awk 'NR>5' |awk '{if ($1==$4 && $2==$5 && $12==$13) {next;} else {print $0}}'  |awk '{print $12,$1,$2"\n"$13,$4,$5}' |awk '{if($2<$3) {print $0}else {p
rint $1,$3,$2}}' |tr " " "\t" |sort --parallel=16 -k1,1V -k2,3n |bedtools merge >mummer.bed
```
### What sizes and how much duplication is present?
```
[remkv6@condo135 50_tandemMummer]$ less mummer.bed |awk '{print $3-$2}' |summary.sh
Total:  85,199,835
Count:  32,624
Mean:   2,611
Median: 490
Min:    34
Max:    347,817
[remkv6@condo135 50_tandemMummer]$ less mummer.bed |awk '{print $3-$2}' |awk '$1>1000'|summary.sh
Total:  76,238,094
Count:  9,561
Mean:   7,973
Median: 2,925
Min:    1,001
Max:    347,817
[remkv6@condo135 50_tandemMummer]$ less mummer.bed |awk '{print $3-$2}' |awk '$1>5000'|summary.sh
Total:  63,198,527
Count:  3,451
Mean:   18,313
Median: 11,800
Min:    5,004
Max:    347,817
[remkv6@condo135 50_tandemMummer]$ less mummer.bed |awk '{print $3-$2}' |awk '$1>10000'|summary.sh
Total:  53,939,814
Count:  2,150
Mean:   25,088
Median: 17,908
Min:    10,009
Max:    347,817
```

#Create gff for overlap analyses
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/32_MummerTandem
ln -s ../../50_tandemMummer/mummer.bed
less scaffold.renamer.sh |sed 's|/|\t|g' |awk '{print $1,$2,$3,$5"/"$6"/"$7"/"$8,"mummer.bed"}' >scaffoldRenamerMummer.sh
sh scaffoldRenamerMummer.sh
awk '{print $1,"Mummer3","duplication", $2,$3,".","+",".",$1"_"$2"_"$3}' mummer.bed |tr " " "\t" >mummer.gff
```
