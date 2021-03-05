# Checking tandem duplications in the dovetail genome

```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/13_Mummer4Self


#this only looks for those duplications over 1000bp
module load GIF/mummer/4.0.0beta2
mummer -mum -l 1000 -b -n -qthreads 12 DovetailSCNMaker4.genome.fasta DovetailSCNMaker4.genome.fasta >SelfAlignmentMummer.out

#what is the total duplication length?
awk -v scaff=1 '{if($1==">") {print $0,scaff=$2} else if($1==scaff  && $2!=$3) {print $0}}' SelfAlignmentMummer.out |awk '$1!=">" {print $4}' |summary.sh
Total:  5,089,970
Count:  2,391
Mean:   2,128
Median: 1,481
Min:    1,000
Max:    54,425

#this doesnt give a real bed file.  it is in mummer format still: start, start, length
awk -v scaff=1 '{if($1==">") {print $0,scaff=$2} else if($1==scaff  && $2!=$3) {print $0}}' SelfAlignmentMummer.out |awk '$1!=">" '>SelfMummer.bed



```
