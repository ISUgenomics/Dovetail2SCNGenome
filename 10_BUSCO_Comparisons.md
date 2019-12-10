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


#run busco
echo " ml busco/3.0.1-py2-cuda9-openmpi3-ze7lkie; ml augustus/3.3-py2-cuda9-openmpi3-fimdyeu;ml dafoam/1.0; ml ncbi-blast/2.4.0+;ml hmmer/3.1b2-cuda9-openmpi3-4ab6zzt; export AUGUSTUS_CONFIG_PATH=/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/05_pseudomolecule/config; export BUSCO_CONFIG_FILE=/work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison/05_pseudomolecule/my-busco.conf; run_BUSCO.py -i SCNgenome.fasta -l /work/GIF/remkv6/Baum/04_Dovetail2Restart/09_BuscoComparison//04_590D2/busco-3.0.1-ze7lkiedvzma2wiiehfdwa7usmcgk5wi/nematoda_odb9 -o PseudoBUSCO -m geno -c 15 -s Hglycines2 -f ">busco.sh

#my-busco.conf
################################################################################
# BUSCO specific configuration
# It overrides default values in code and dataset cfg, and is overridden by arguments in command line
# Uncomment lines when appropriate
[busco]
# Input file
;in = ./sample_data/target.fa
# Run name, used in output files and folder
;out = SAMPLE
# Where to store the output directory
;out_path = ./sample_data
# Path to the BUSCO dataset
;lineage_path = ./sample_data/example
# Which mode to run (genome / protein / transcriptome)
;mode = genome
# How many threads to use for multithreaded steps
;cpu = 1
# Domain for augustus retraining, eukaryota or prokaryota
# Do not change this unless you know exactly why !!!
;domain = eukaryota
# Force rewrite if files already exist (True/False)
;force = False
# Restart mode (True/False)
;restart = False
# Blast e-value
;evalue = 1e-3
# Species to use with augustus, for old datasets only
;species = fly
# Augustus extra parameters
# Use single quotes, like this: '--param1=1 --param2=2'
;augustus_parameters = ''
# Tmp folder
;tmp_path = ./tmp/
# How many candidate regions (contigs, scaffolds) to consider for each BUSCO
;limit = 3
# Augustus long mode for retraining (True/False)
;long = False
# Quiet mode (True/False)
;quiet = False
# Debug logs (True/False), it needs Quiet to be False
;debug = True
# tar gzip output files (True/False)
;gzip = False
# Force single core for the tblastn step
;blast_single_core = True

[tblastn]
# path to tblastn
path =  /opt/rit/app/ncbi-blast/2.4.0+/bin/

[makeblastdb]
# path to makeblastdb
path =  /opt/rit/app/ncbi-blast/2.4.0+/bin/

[augustus]
# path to augustus
path = /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/augustus-3.3-fimdyeurm35h63s5mq7lqyrxjryhn3ks/bin/

[etraining]
# path to augustus etraining
path = /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/augustus-3.3-fimdyeurm35h63s5mq7lqyrxjryhn3ks/bin/

# path to augustus perl scripts, redeclare it for each new script
[gff2gbSmallDNA.pl]
path = /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/augustus-3.3-fimdyeurm35h63s5mq7lqyrxjryhn3ks/scripts/
[new_species.pl]
path = /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/augustus-3.3-fimdyeurm35h63s5mq7lqyrxjryhn3ks/scripts/
[optimize_augustus.pl]
path = /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/augustus-3.3-fimdyeurm35h63s5mq7lqyrxjryhn3ks/scripts/

[hmmsearch]
# path to HMMsearch executable
path = /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/hmmer-3.1b2-4ab6zzt2z5x2xp2osta2cobfb4dbh3b7/bin/

[Rscript]
# path to Rscript, if you wish to use the plot tool
path = /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/r-3.4.3-kaltwmmc7x5pobe6lzwecoicwod5ntpm/bin/
################################################################################

```
