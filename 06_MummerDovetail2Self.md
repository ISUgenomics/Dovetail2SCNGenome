# How many duplications are present in the dovetail 2 genomes


### run mummer
```
sed  's/>/>2/g' nematode_sp_19Jul2018_IbtP1.fasta >Genome2.fasta
ln -s ../16_NewDovetailGenome/nematode_sp_19Jul2018_IbtP1.fasta

module load GIF/mummer/4.0.0beta2
mummer -mum -l 1000 -b -n -qthreads 12  nematode_sp_19Jul2018_IbtP1.fasta Genome2.fasta >2vs2.mummer.out
```

### How larger are the duplications
```
less 2vs2.mummer.out |sed 's/2Scaffold/Scaffold/g' |awk -v scaff=1 '{if($1==">") {print $0,scaff=$2} else if($1==scaff  && $2!=$3) {print $0}}' |awk '$1!=">" {print $4}' |summary.sh
Total:  3,430,752
Count:  1,909
Mean:   1,797
Median: 1,373
Min:    1,000
Max:    19,730


#here is an actual representation of the duplications in the genome, in a bed format.  However, it will need to be modified to show which duplicated fragment belongs to which.  
less 2vs2.mummer.out |sed 's/2Scaffold/Scaffold/g' |awk -v scaff=1 '{if($1==">") {print $0,scaff=$2} else if($1==scaff  && $2!=$3) {print $0}}' |awk '$1!=">" ' |awk -v scaff=1 '{print $1,$2,$2+$4"\n"$1,$3,$3+$4}' |sort -k1,1V -k2,3n | tr " " "\t" |bedtools merge -d 500 >SelfMummer.bed

less SelfMummer.bed |awk '{print $3-$2}' |summary.sh
Total:  6,560,828
Count:  2,761
Mean:   2,376
Median: 1,565
Min:    1,000
Max:    26,508

```
