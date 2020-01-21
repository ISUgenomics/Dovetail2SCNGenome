# Need a better estimate of synteny that does not rely only on gene order, i.e. Mummer

### Draft 738 to pseudomolecule tn10
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/34_Mummer

ln -s ../10_RepeatModeler/SCNgenome.fasta
ln -s ../../CamTechGenomeComparison/58_Renamatorium/1_genomeNgff/genome738sl.polished.mitoFixed.fa

ml bedtools2/2.27.1-s2mtpsu
ml mummer/3.23

 run-mummer3  genome738sl.polished.mitoFixed.fa SCNgenome.fasta 7382pseudo

 #500bp segments minimum, redundancy removed with bedtools
  less 7382pseudo.out |awk '$4>500 {print $1,$3,$3+$4}' |tr " " "\t" |sort -k1,1V -k2,3n |bedtools merge -i - |awk '{print $3-$2}'
  |summary.sh
 Total:  92,397,944
 Count:  57,232
 Mean:   1,614
 Median: 904
 Min:    501
 Max:    50,951



```
### X12 assembly
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/34_Mummer

ln -s ../28_X12GenomeComparison/SCN_genome.fa X12genome.fasta

 run-mummer3  X12genome.fasta SCNgenome.fasta X122pseudo

#500bp segments minimum, redundancy removed with bedtools
 less X122pseudo.out |awk '$4>500 {print $1,$3,$3+$4}' |tr " " "\t" |sort -k1,1V -k2,3n |bedtools merge -i - |awk '{print $3-$2}'
  |summary.sh
 Total:  43,435,617
 Count:  42,141
 Mean:   1,030
 Median: 757
 Min:    501
 Max:    28,247


```
