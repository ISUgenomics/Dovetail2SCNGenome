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

```




### use diamond to find conserved effector proteins
```


ml diamond/0.9.23-xqnzcyt

#blastx with 50% min query cover and 50% min identity
diamond blastx --query effector.fa  -d mikado_proteinsFixed --in mikado_proteinsFixed.fasta --strand both --query-cover .5
  --id 0.5 >diamond.out

# how many primary isoforms meet the above criteria?
less diamond.out |awk 'substr($2,length($2),2)=="1"' |wc

```
### Identify secreted proteins
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/32_SignalP
ml dafoam/1.0
ml signalp/4.1

signalp -f summary mikado_proteinsFixed.fasta >mikado_proteinsFixed.fasta.out

#how many secreted without transmembrane domains
less mikado_proteinsFixed.fasta.out |grep "Name=" |grep "YES" |grep "noTM" |wc
   3267   42471  420764
```