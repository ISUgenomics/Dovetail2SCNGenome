#need some comparisons of busco for PAG 2019 between the four assemblies.


These will be just genome runs, no augustus on the nematoda dataset
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison
mkdir 01_738
mkdir 02_104D1
mkdir 03_1544FFU
mkdir 04_590D2
```

738 genome
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/01_738
ln -s ../../../CamTechGenomeComparison/58_Renamatorium/1_genomeNgff/genome738sl.polished.mitoFixed.fa

module use /work/GIF/software/modules
module load GIF2/ncbi-blast/2.2.30+
module load hmmer
module unload augustus/3.3-py2-cuda9-openmpi3-fimdyeu
module load GIF/augustus/3.3.2

/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/04_590D2/busco-3.0.1-ze7lkiedvzma2wiiehfdwa7usmcgk5wi/bin/run_BUSCO.py -i genome738sl.polished.mitoFixed.fa  -l /work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison//04_590D2/busco-3.0.1-ze7lkiedvzma2wiiehfdwa7usmcgk5wi/nematoda_odb9 -o test -m geno -c 16 -sp Hglycines


INFO    C:66.3%[S:54.7%,D:11.6%],F:7.6%,M:26.1%,n:982
INFO    651 Complete BUSCOs (C)
INFO    537 Complete and single-copy BUSCOs (S)
INFO    114 Complete and duplicated BUSCOs (D)
INFO    75 Fragmented BUSCOs (F)
INFO    256 Missing BUSCOs (M)
INFO    982 Total BUSCO groups searched


```

scaffolded 738 genome -- dovetail 1
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/02_104D1
 ln -s ../../../01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/DovetailSCNMaker4.genome.fasta
 module use /work/GIF/software/modules
 module load GIF2/ncbi-blast/2.2.30+
 module load hmmer
 module unload augustus/3.3-py2-cuda9-openmpi3-fimdyeu
 module load GIF/augustus/3.3.2

 /work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/04_590D2/busco-3.0.1-ze7lkiedvzma2wiiehfdwa7usmcgk5wi/bin/run_BUSCO.py -i DovetailSCNMaker4.genome.fasta  -l /work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison//04_590D2/busco-3.0.1-ze7lkiedvzma2wiiehfdwa7usmcgk5wi/nematoda_odb9 -o test -m geno -c 15 -sp Hglycines

 C:65.8%[S:57.3%,D:8.5%],F:7.8%,M:26.4%,n:982
 646 Complete BUSCOs (C)
 563 Complete and single-copy BUSCOs (S)
 83 Complete and duplicated BUSCOs (D)
 77 Fragmented BUSCOs (F)
 259 Missing BUSCOs (M)
 982 Total BUSCO groups searched
```

Falcon -- Falcon Unzip assembly
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/03_1544FFU

module use /work/GIF/software/modules
module load GIF2/ncbi-blast/2.2.30+
module load hmmer
module unload augustus/3.3-py2-cuda9-openmpi3-fimdyeu
module load GIF/augustus/3.3.2

/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/04_590D2/busco-3.0.1-ze7lkiedvzma2wiiehfdwa7usmcgk5wi/bin/run_BUSCO.py -i all_p_ctg_preprocess.fa  -l /work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison//04_590D2/busco-3.0.1-ze7lkiedvzma2wiiehfdwa7usmcgk5wi/nematoda_odb9 -o test1 -m geno -c 15 -sp Hglycines -f

C:63.3%[S:54.8%,D:8.5%],F:6.8%,M:29.9%,n:982
621 Complete BUSCOs (C)
538 Complete and single-copy BUSCOs (S)
83 Complete and duplicated BUSCOs (D)
67 Fragmented BUSCOs (F)
294 Missing BUSCOs (M)
982 Total BUSCO groups searched
```

#Dovetail2 -- Falcon/Falcon_Unzip assembly
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/04_590D2
ln -s ../../01_GenomeDownload/nematode_sp_19Jul2018_IbtP1.fasta

module load GIF2/ncbi-blast/2.2.30+
module load hmmer
module unload augustus/3.3-py2-cuda9-openmpi3-fimdyeu
module load GIF/augustus/3.3.2

#copied entire busco3 directory to my work directory, as I did not have write permission.  Had to make a dozen changes to config.ini.default


busco-3.0.1-ze7lkiedvzma2wiiehfdwa7usmcgk5wi/bin/run_BUSCO.py -i nematode_sp_19Jul2018_IbtP1.fasta  -l busco-3.0.1-ze7lkiedvzma2wiiehfdwa7usmcgk5wi/nematoda_odb9 -o test -m geno -c 15 -sp Hglycines -f

C:62.6%[S:56.4%,D:6.2%],F:7.2%,M:30.2%,n:982
615 Complete BUSCOs (C)
554 Complete and single-copy BUSCOs (S)
61 Complete and duplicated BUSCOs (D)
71 Fragmented BUSCOs (F)
296 Missing BUSCOs (M)
982 Total BUSCO groups searched

```

### Pseudomolecule assembly after juicebox and pilon iterations
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/05_pseudomolecule

ln -s ../../23_Mikado/SCNgenome.fasta

cp -rf /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/augustus-3.3-fimdyeurm35h63s5mq7lqyrxjryhn3ks/config/ .
cp -rf Hglycines2/ ../../../../09_BuscoComparison/05_pseudomolecule/config/species/.

ml miniconda2
conda create -n busco
source activate busco
conda install -c bioconda/label/cf201901 busco


ml miniconda2; source activate busco; export AUGUSTUS_CONFIG_PATH=/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/05_pseudomolecule/config; export ; run_BUSCO.py -i SCNgenome.fasta -l /work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison//04_590D2/busco-3.0.1-ze7lkiedvzma2wiiehfdwa7usmcgk5wi/nematoda_odb9 -o PseudoBUSCO -m geno -c 15 -s Hglycines2 -f

C:64.6%[S:59.9%,D:4.7%],F:8.8%,M:26.6%,n:982
634 Complete BUSCOs (C)
588 Complete and single-copy BUSCOs (S)
46 Complete and duplicated BUSCOs (D)
86 Fragmented BUSCOs (F)
262 Missing BUSCOs (M)
982 Total BUSCO groups searched
```

 ### Run busco on annotated proteins
 ```
 #/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/06_BuscoPseudomoleculeProteins

 ml cdbfasta
 cdbfasta mikado_proteinsFixed.fasta

less mikado_proteinsFixed.fasta |grep ".1" |awk '{print $1}' |sed 's/>//g'|cdbyank mikado_proteinsFixed.fasta.cidx >PrimaryIsoformsMikado_proteinsFixed.fasta

 module use /work/GIF/software/modules
 ml miniconda2; source activate busco; export AUGUSTUS_CONFIG_PATH=/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/05_pseudomolecule/config;  run_BUSCO.py -i PrimaryIsoformsMikado_proteinsFixed.fasta -l /work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison//04_590D2/busco-3.0.1-ze7lkiedvzma2wiiehfdwa7usmcgk5wi/nematoda_odb9 -o PseudoBUSCO -m prot -c 15 -s Hglycines2 -f

 C:68.6%[S:56.2%,D:12.4%],F:9.0%,M:22.4%,n:982
 INFO    674 Complete BUSCOs (C)
 INFO    552 Complete and single-copy BUSCOs (S)
 INFO    122 Complete and duplicated BUSCOs (D)
 INFO    88 Fragmented BUSCOs (F)
 INFO    220 Missing BUSCOs (M)
 INFO    982 Total BUSCO groups searched

 ```
### Check buscos via blast to see how many are missing
```
 ln -s  /work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison//04_590D2/busco-3.0.1-ze7lkiedvzma2wiiehfdwa7usmcgk5wi/nematoda_odb9/ancestral
ln -s ../06_BuscoPseudomoleculeProteins/mikado_proteinsFixed.fasta
ln -s ../06_BuscoPseudomoleculeProteins/mikado_proteinsFixed.fasta.cidx
makeblastdb -in mikado_proteinsFixed.fasta -dbtype prot

cdbfasta ancestral
blastp -db mikado_proteinsFixed.fasta -query ancestral -num_threads 16 -outfmt 6 -out ancestral2SCN.blastout

#create table adding protein lengths
sort -u -k1,1 ancestral2SCN.blastout |awk '{print $1}' |cdbyank ancestral.cidx |bioawk -c fastx '{print $name,length($seq)}' |paste <(sort -u -k1,1 ancestral2SCN.blastout) - >AddBUSCOlength.blastout

less AddBUSCOlength.blastout |awk '{print $2}' |while read line; do grep -m 1 -w $line mikado.lengths ;done|paste AddBUSCOlength.blastout - >AddBUSCOandMikadolength.blastout

 #how many busco genes were found with only an evalue constraint
 less AddBUSCOandMikadolength.blastout |awk '$11<.01' |wc
    948   15168  118326

#How many busco genes were found with an evalue constraint and a length constraint of 1.6x to .4x maximum?
less AddBUSCOandMikadolength.blastout |awk '$11<.01' |awk '{print $14/$16}' |awk '$1<1.6 && $1>.4' |wc
    788     788    6915



By just evalue constraints, missing only from protein busco..
less missing_busco_list_PseudoBUSCO.tsv |grep -v "#" |grep -f - ../../09_ManualBlastsOfBuscoGenes/AddBUSCOandMikadolength.blastout |awk '$11<.01' |wc
    187    2992   23656


What is the median evalue of those genes that were called as missing?
less missing_busco_list_PseudoBUSCO.tsv |grep -v "#" |grep -f - ../../09_ManualBlastsOfBuscoGenes/AddBUSCOandMikadolength.blastout |awk '$11<.01' |awk '{print $11*1000000000000000000000000000000}' |summary.sh
Total:  24,926,427,861,843,795,432,106,360,832
Count:  90
Mean:   276,960,309,576,042,158,487,502,848
Median: 5,800,000,000,000
Min:    0
Max:    8,999,999,999,999,999,844,710,088,704

awk '{print 5800000000000/1000000000000000000000000000000}' missing_busco_list_PseudoBUSCO.tsv |less
5.8e-18

#Busco total by addition of .01 and 1.6x and 0.4x length
less missing_busco_list_PseudoBUSCO.tsv |grep -v "#" |grep -f - ../../09_ManualBlastsOfBuscoGenes/AddBUSCOandMikadolength.blastout |awk '$11<.01' |awk '{print $14/$16}' |awk '$1<1.6 && $1>.4' |wc
   141     141    1225
 141 + 674 = 815 == 815/982 == 82.99%
```
