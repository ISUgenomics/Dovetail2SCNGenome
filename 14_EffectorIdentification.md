#How many wet-bench verified effectors are present in the pseudomolecule genome.  


### Gmap the known effectors
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/29_Effectors

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
less SCNgenome.effector.gff3 |awk '$3=="CDS"' |bedtools intersect -wo -b - -a mikado.loci.ancestral.gff3  |sed 's/ID=//g' |sed 's/;/\t/g' |awk '$3=="mRNA"' |cut -f 9,23 |sed 's/Target=//g' |sed 's/Name=//g' |sed 's/Parent=//g' |awk '{print $1"\t"$2}' |sort -k1,1 -u |sed 's/lcl|//g' >GmapEffectorsGenes.tab
wc GmapEffectorsGenes.tab # these are mRNAs
125  250 8249 GmapEffectorsGenes.tab


#How many in the old genome?
awk '$3=="gene"' Sortedeffectors.gmapped.gff3 |awk '{print $1,$4,$5}' |tr " " "\t" |sort -k1,1V -k2,3n |bedtools merge -i - |wc
    100     300    2546
```


### use diamond to find conserved effector proteins
```
ml diamond/0.9.23-xqnzcyt
diamond makedb --in mikado.loci.ancestralVHEJ_proteins.fasta -d mikado.loci.ancestralVHEJ_proteins

#blastx with 50% min query cover and 50% min identity
diamond blastx --query effector.fa  -d mikado.loci.ancestralVHEJ_proteins --in mikado.loci.ancestralVHEJ_proteins.fasta --strand both --query-cover .5 --id 0.5 >diamond.out

# how many mRNA's meet the above criteria?
awk '{print $2}' diamond.out |sort|uniq|wc
362     362    8989


awk '{print $2"\t"$1}' diamond.out |sort -u  -k1,1 >NamedDiamondEffectors.tab


less ../mikado.loci.ancestral.gff3 |sed 's/ID=/ID=\t/g' |sed 's/transcripts=/transcripts=\t/g' |sed 's/;/\t/1' >mikado.loci.ancestral.GFFGREPMOD

#get all mrna, cds, exons for diamond effectors
awk '{print $1}' NamedDiamondEffectors.tab|while read line; do  grep -w  $line mikado.loci.ancestral.GFFGREPMOD >>DiamondEffectors.gffmod;done &
#grab genes
less DiamondEffectors.gffmod |awk '$3=="locus"' |paste - NamedDiamondEffectors.tab |cut -f 1-12,14 |awk -F"\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9$10":"$11$12";Note="$13}' |tr " " "\t" >DiamondEffectorsGene.gff

#concat genes & all mrna, cds, exons for diamond effectors
cat DiamondEffectorsGene.gff <(awk '$3!="locus"' DiamondEffectors.gffmod) |sed 's/ID=\t/ID=/g' |sed 's/Parent=\t/Parent=/g' >UnsortedDiamondEffectors.gff


perl gff3sort/gff3sort.pl --precise --chr_order natural TidyUnsortedDiamondEffectors.gff >SortedTidyUnsortedDiamondEffectors.gff
sh ~/common_scripts/runTabix.sh SortedTidyUnsortedDiamondEffectors.gff


less SortedTidyUnsortedDiamondEffectors.gff |awk '$3=="gene"' |wc
    402    3618   42698

```
