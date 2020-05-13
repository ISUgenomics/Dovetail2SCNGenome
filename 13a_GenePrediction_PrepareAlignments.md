## Prepare all alignments and masking for gene prediction

### Repeatmodeler and Repeatmasker
```
# /work/GIF/remkv6/Baum/04_Dovetail2Restart/10_RepeatModeler

module use /work/GIF/software/modules
sh runRepeatModeler.sh SCNgenome.fasta

################################################################################
#!/bin/bash
# runs repeat masking for the genome after constructing custom repeat library
# uses repeat modeler for building custom db and RepeatMasking for masking
# run it as:
# runRepeatModeler.sh Genome.fasta
# based on Rick's guide https://intranet.gif.biotech.iastate.edu/doku.php/people:remkv6:genome738polished_repeatmodeler_--de_novo_repeat_identification

if [ $# -lt 1 ] ; then
        echo "usage: runRepeatModeler <genome.fasta>"
        echo ""
        echo "To build custom repeat library and mask the repeats of the genome"
        echo ""
exit 0
fi


GENOME="$1"
module use /shared/software/GIF/modules/
module purge
module load parallel

module load GIF2/repeatmasker/4.0.6
module load GIF2/repeatmodeler/1.0.8
module load GIF2/perl/5.22.1
DATABASE="$(basename ${GENOME%.*}).DB"
BuildDatabase -name ${DATABASE} -engine ncbi ${GENOME}
RepeatModeler -database ${DATABASE}  -engine ncbi -pa 16
ln -s $(find $(pwd) -name "consensi.fa.classified")
RepeatMasker -pa 16 -gff -lib consensi.fa.classified ${GENOME}
################################################################################

#Results of repeatmodeler/masker
==================================================
file name: SCNgenome.fasta
sequences:             9
total length:  157982452 bp  (156301001 bp excl N/X-runs)
GC level:         36.66 %
bases masked:   61438776 bp ( 38.89 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
SINEs:               67        13915 bp    0.01 %
      ALUs            0            0 bp    0.00 %
      MIRs            0            0 bp    0.00 %

LINEs:             5375      1795726 bp    1.14 %
      LINE1         422        69365 bp    0.04 %
      LINE2           0            0 bp    0.00 %
      L3/CR1       2654      1265055 bp    0.80 %

LTR elements:     10759      6004426 bp    3.80 %
      ERVL            0            0 bp    0.00 %
      ERVL-MaLRs      0            0 bp    0.00 %
      ERV_classI     35        49865 bp    0.03 %
      ERV_classII   757       149391 bp    0.09 %

DNA elements:     70403     14691523 bp    9.30 %
     hAT-Charlie      0            0 bp    0.00 %
     TcMar-Tigger    23        23125 bp    0.01 %

Unclassified:    154224     35100826 bp   22.22 %

Total interspersed repeats: 57606416 bp   36.46 %


Small RNA:          295       107683 bp    0.07 %

Satellites:         775       200791 bp    0.13 %
Simple repeats:   53062      2878802 bp    1.82 %
Low complexity:   17729      1257503 bp    0.80 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element


The query species was assumed to be homo
RepeatMasker version open-4.0.6 , default mode

run with rmblastn version 2.2.27+
The query was compared to classified sequences in "consensi.fa.classified"
RepBase Update 20160829, RM database version 20160829
```

### Align RNAseq reads
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA

for f in /work/GIF/remkv6/Baum/03_GlandRepeat738/01_Alignment2Genome/*val_*.fq.gz;do ln -s $f;done
for f in /work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm/*val_*.fq.gz;do ln -s $f;done
ln -s ../10_RepeatModeler/SCNgenome.fasta

paste <(ls -1 *val_1*)  <(ls -1 *val_2*) |while read a b ; do echo "sh runFeatureCounts.sh "$a" "$b" /work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA SCNgenome.fasta";done >align.sh

ls -l *val_1* |awk '{print $11}' |sed  's|/|\t|g' |cut -f 8 |sed 's/\.gz/_sorted.bam/g' |tr "\n" " " |awk '{print "samtools merge AllRNASEQ_sorted.bam "$0}' |less

samtools merge AllRNASEQ_sorted.bam 1703FL-02-01_S1_L001_R1_001_val_1.fq_sorted.bam 1703FL-02-02_S2_L001_R1_001_val_1.fq_sorted.bam 1703FL-02-03_S3_L001_R1_001_val_1.fq_sorted.bam 1703FL-02-04_S4_L001_R1_001_val_1.fq_sorted.bam 1703-TM101_S0_L003_R1_001_val_1.fq_sorted.bam 1703-TM102_S0_L003_R1_001_val_1.fq_sorted.bam SRR6230579_1_val_1.fq_sorted.bam SRR6230580_1_val_1.fq_sorted.bam SRR6230581_1_val_1.fq_sorted.bam SRR6230582_1_val_1.fq_sorted.bam SRR6230583_1_val_1.fq_sorted.bam SRR6230584_1_val_1.fq_sorted.bam SRR6230585_1_val_1.fq_sorted.bam SRR6230586_1_val_1.fq_sorted.bam SRR6230587_1_val_1.fq_sorted.bam
```

### Align isoseq and ESTs
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA

ln -s /work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/6_isoseq/consensus_isoforms.fasta
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/09_Maker/H.glycinesEST.fasta
sh runGmap.sh SCNgenome /work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA SCNgenome.fasta consensus_isoforms.fasta
sh runGmap.sh SCNgenome /work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA SCNgenome.fasta H.glycinesEST.fasta
```


### Filter the rRNA out of the bam file, as it is excessive in gland reads
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/20_ribosomalArrayID

#ribosomal sequences from ncbi nucleotide (tylenchida) and fastas from the old genome
cat ../../01_SCNDovetailScaffolding/08_ribosomalArrays/RDNANucleotide.fasta ../../01_SCNDovetailScaffolding/08_ribosomalArrays/SCNSpecificRibosomalSeq.fasta >RibosomalNCBINOldGenome.fasta

#get prospective coordinates for rDNA arrays
module load blast-plus
makeblastdb -in  SCNgenome.fasta -dbtype nucl -out SCNgenome.blastdb
blastn -db SCNgenome.blastdb -query RibosomalNCBINOldGenome.fasta -num_threads 16 -outfmt 6 -out ribo2scn.blastout

module load bedtools2
less ribo2scn.blastout |awk '$12>200 {print $2,$9,$10}'|awk '{if($2>$3) {print $1,$3,$2} else {print $0}}' |sort|uniq|tr " " "\t" |bedtools merge -d 2000 -i - >RiboCoords.bed

module load samtools
samtools view ../19_brakerMasked/AllRNASEQ_sorted.bam -b -h -o rdnaReads.bam -U NonRiboReads.bam -L RiboCoords.bed
```


### Create Tylenchida EST database
```
sh runGmap.sh SCNgenome /work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn/ SCNgenome.fasta TylenchidaESTNotH.glycines.fasta

#runGmap.sh
###################################################################################
#!/bin/bash

#Makes a database and searches your sequences.
#sh runGmap.sh <database name> <folder of database file ending with a "/"> <Fasta file> <query file>

#examples
#sh run_gmap.sh red_abalone_02Jun2017_5fUJu /work/GIF/remkv6/Serb/03_DavideGMAP/ red_abalone_02Jun2017_5fUJu.fasta DavideQuerydna.fasta
#sh run_gmap.sh  m.yessoensisGenome /work/GIF/remkv6/Serb/03_DavideGMAP DavideQuerydna.fasta
#sh run_gmap.sh Crassostreagigasgenome /work/GIF/remkv6/Serb/03_DavideGMAP Crassostreagigasgenome.fa DavideQuerydna.fasta


#module load gsnap
dbname=$1
dbloc=$2
dbfasta=$3
query=$4
gmap_build -d $dbname  -D $dbloc $dbfasta
gmap -D $dbloc -d $dbname -B 5 -t 16  --input-buffer-size=1000000 --output-buffer-size=1000000 -f gff3_gene $query >${dbname%.*}.${query%.*}.gff3

###################################################################################

#How many EST's were there?
grep -c ">"   TylenchidaESTNotH.glycines.fasta
591134
How many mapped?
awk '$3=="gene"' SCNgenome.TylenchidaESTNotH.glycines.gff3 |wc
  61094  549846 5384817

#These are the counts from the related species gene preidctions I used in the tylenchida.  
awk '$3=="gene"' SCNgenome.TylenchidaESTNotH.glycines.gff3 |cut -f 9 |sed 's/ID=//g' |sed 's/_/\t/1' |cut -f 1 |sort|uniq -c |sort -k1,1nr |less
  8305 GPLIN
  8248 GROS
   372 MhA1
   248 BXY
   181 augustus
   150 snap
   134 Dd
...
#Another 43000 that were actual ESTs from NCBI's library
awk '$3=="gene"' SCNgenome.TylenchidaESTNotH.glycines.gff3 |cut -f 9 |sed 's/ID=//g' |sed 's/_/\t/1' |cut -f 1 |sort|uniq -c |sort -k1,1nr |wc
  43187
```
### Gmap alignments of previous 368 scaffold assembly's, Maker annotation's transcripts

```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/24_MapDove1MakerGenes

ln -s ../../01_SCNDovetailScaffolding/09_Maker/01_maker/DovetailSCNMaker4.all.maker.transcripts.fasta
ln -s ../10_RepeatModeler/SCNgenome.fasta

sh runGmap.sh SCNgenome /work/GIF/remkv6/Baum/04_Dovetail2Restart/24_MapDove1MakerGenes/ SCNgenome.fasta DovetailSCNMaker4.all.maker.transcripts.fasta


awk '$3=="mRNA"' SCNgenome.DovetailSCNMaker4.all.maker.transcripts.gff3 |grep -v "#" |wc
 78024  702216 15488896



```
