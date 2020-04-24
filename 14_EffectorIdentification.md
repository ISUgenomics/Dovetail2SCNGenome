#How many wet-bench verified effectors are present in the pseudomolecule genome.  


### Gmap the known effectors
```

#runGmap.sh at 50% min query coverage
##########################################################
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
#gmap_build -d $dbname  -D $dbloc $dbfasta
gmap -D $dbloc -d $dbname -B 5 -t 16  --min-trimmed-coverage=0.5 --input-buffer-size=1000000 --output-buffer-size=1000000 -f gff3_gene  $query >${dbname%.*}.${query%.*}.gff3
#################################################################

#how many alignment of 50% query coverage did I get?
awk '$3=="gene"' SCNgenome.effector.gff3 |wc
   126    1134   18667

#prep tab file of gene"\t"effector
less SCNgenome.effector.gff3 |awk '$3=="CDS"' |bedtools intersect -wo -b - -a OrderedSCNGenePredictions.gff3  |sed 's/ID=//g' |sed 's/;/\t/g' |awk '$3=="mRNA"' |cut -f 9,23 |sed 's/Target=//g' |sed 's/Name=//g' |sed 's/Parent=//g' |awk '{print $1"\t"$2}' |sort -k1,1 -u |sed 's/lcl|//g' >GmapEffectorsGenes.tab
wc GmapEffectorsGenes.tab # these are mRNAs
 156  312 7331 GmapEffectorsGenes.tab

#How many in the old genome?
awk '$3=="gene"' Sortedeffectors.gmapped.gff3 |awk '{print $1,$4,$5}' |tr " " "\t" |sort -k1,1V -k2,3n |bedtools merge -i - |wc
    100     300    2546
```


### use diamond to find conserved effector proteins
```
ml diamond/0.9.23-xqnzcyt
diamond makedb --in OrderedSCNGenePredictionsVHEJ_proteins.fasta -d OrderedSCNGenePredictionsVHEJ_proteins

#blastx with 50% min query cover and 50% min identity
diamond blastx --query effector.fa  -d OrderedSCNGenePredictionsVHEJ_proteins --in OrderedSCNGenePredictionsVHEJ_proteins.fasta --strand both --query-cover .5 --id 0.5 >diamond.out

# how many mRNA's meet the above criteria?
awk '{print $2}' diamond.out |sort|uniq|wc
    431     431    4178

awk '{print $2"\t"$1}' diamond.out |sort -u  -k1,1 >NamedDiamondEffectors.tab
```
