# Annotating Ribosomal arrays in the dovetail scaffolded scn genome

Downloaded all EST's that were annotated as H. glycines rDNA

```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/08_ribosomalArrays

#Used the runGMAP.sh script
#!/bin/bash

#Makes a database and searches your sequences.
#sh runGmap.sh <database name> <folder of database file ending with a "/"> <Fasta file> <query file>

#examples
#sh run_gmap.sh red_abalone_02Jun2017_5fUJu /work/GIF/remkv6/Serb/03_DavideGMAP/ red_abalone_02Jun2017_5fUJu.fasta DavideQuerydna.fasta
#sh run_gmap.sh  m.yessoensisGenome /work/GIF/remkv6/Serb/03_DavideGMAP DavideQuerydna.fasta
#sh run_gmap.sh Crassostreagigasgenome /work/GIF/remkv6/Serb/03_DavideGMAP Crassostreagigasgenome.fa DavideQuerydna.fasta


module load gsnap
dbname=$1
dbloc=$2
dbfasta=$3
query=$4
gmap_build -d $dbname  -D $dbloc $dbfasta
gmap -D $dbloc -d $dbname -B 5 -t 12  --input-buffer-size=1000000 --output-buffer-size=1000000 -f gff3_gene  $query >${dbname%.*}.${query%.*}.gff


sh run_gmap.sh nematodeDB /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/08_ribosomalArrays/ nematode_sp._22Aug2017_DZkUC.fasta RDNANucleotide.fasta


Manually merged the output of the above to get three loci
#less RibosomalCoord.bed
Scaffold_25;HRSCAF=85   2038969 2039583
Scaffold_25;HRSCAF=85   3870261 3911828
Scaffold_70;HRSCAF=275  8174    33130
```
