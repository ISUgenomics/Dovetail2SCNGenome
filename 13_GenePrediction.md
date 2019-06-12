# Need to predict genes for the SCN TN10 pseudomolecule assembly

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

# Align reads
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA

for f in /work/GIF/remkv6/Baum/03_GlandRepeat738/01_Alignment2Genome/*val_*.fq.gz;do ln -s $f;done
for f in /work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm/*val_*.fq.gz;do ln -s $f;done
ln -s ../10_RepeatModeler/SCNgenome.fasta


 paste <(ls -1 *val_1*)  <(ls -1 *val_2*) |while read a b ; do echo "sh runFeatureCounts.sh 16 "$a" "$b" /work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA SCNgenome.fasta";done >align.sh

ls -l *val_1* |awk '{print $11}' |sed  's|/|\t|g' |cut -f 8 |sed 's/\.gz/_sorted.bam/g' |tr "\n" " " |awk '{print "samtools merge AllRNASEQ_sorted.bam "$0}' |less
 samtools merge AllRNASEQ_sorted.bam 1703FL-02-01_S1_L001_R1_001_val_1.fq_sorted.bam 1703FL-02-02_S2_L001_R1_001_val_1.fq_sorted.bam 1703FL-02-03_S3_L001_R1_001_val_1.fq_sorted.bam 1703FL-02-04_S4_L001_R1_001_val_1.fq_sorted.bam 1703-TM101_S0_L003_R1_001_val_1.fq_sorted.bam 1703-TM102_S0_L003_R1_001_val_1.fq_sorted.bam SRR6230579_1_val_1.fq_sorted.bam SRR6230580_1_val_1.fq_sorted.bam SRR6230581_1_val_1.fq_sorted.bam SRR6230582_1_val_1.fq_sorted.bam SRR6230583_1_val_1.fq_sorted.bam SRR6230584_1_val_1.fq_sorted.bam SRR6230585_1_val_1.fq_sorted.bam SRR6230586_1_val_1.fq_sorted.bam SRR6230587_1_val_1.fq_sorted.bam

```
