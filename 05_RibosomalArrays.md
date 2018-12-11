#  Need to identify all of the ribosomal arrays in the new Dovetail2 genome

```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/21_Dovetail2ribosomalArrays

Downloaded all EST's that were annotated as H. glycines rDNA

sh runGmap.sh nematode.db /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/21_Dovetail2ribosomalArrays nematode_sp_19Jul2018_IbtP1.fasta RDNANucleotide.fasta


#########################################################################################################################
#!/bin/bash

#Makes a database and searches your sequences.
#sh runGmap.sh <database name> <folder of database file ending with a "/"> <Fasta file> <query file>

#examples
#sh run_gmap.sh red_abalone_02Jun2017_5fUJu /work/GIF/remkv6/Serb/03_DavideGMAP/ red_abalone_02Jun2017_5fUJu.fasta DavideQuerydna.fasta
#sh run_gmap.sh  m.yessoensisGenome /work/GIF/remkv6/Serb/03_DavideGMAP DavideQuerydna.fasta
#sh run_gmap.sh Crassostreagigasgenome /work/GIF/remkv6/Serb/03_DavideGMAP Crassostreagigasgenome.fa DavideQuerydna.fasta

module use /opt/rit/modules/
module load gsnap/20151120
dbname=$1
dbloc=$2
dbfasta=$3
query=$4
gmap_build -d $dbname  -D $dbloc $dbfasta
gmap -D $dbloc -d $dbname -B 5 -t 12  --input-buffer-size=1000000 --output-buffer-size=1000000 -f gff3_gene  $query >${dbname%.*}.${query%.*}.gff

#########################################################################################################################
```

### Dissect where the ribosomal arrays actually are.
```
less nematode.RDNANucleotide.gff |awk '$3=="gene"' |awk '{print $1,$4,$5}' |sort -u -k1,1 |less

########################################
Scaffold_127;HRSCAF=254 3350004 3354134
Scaffold_254;HRSCAF=414 6123 7145
Scaffold_339;HRSCAF=514 3640 4664
Scaffold_370;HRSCAF=552 12294 13317
Scaffold_372;HRSCAF=554 10308 11327
Scaffold_45;HRSCAF=141 45082 46098
Scaffold_474;HRSCAF=809 4221823 4222846
########################################

```
