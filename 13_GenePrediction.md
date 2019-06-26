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

### Align reads
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA

for f in /work/GIF/remkv6/Baum/03_GlandRepeat738/01_Alignment2Genome/*val_*.fq.gz;do ln -s $f;done
for f in /work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm/*val_*.fq.gz;do ln -s $f;done
ln -s ../10_RepeatModeler/SCNgenome.fasta


 paste <(ls -1 *val_1*)  <(ls -1 *val_2*) |while read a b ; do echo "sh runFeatureCounts.sh "$a" "$b" /work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA SCNgenome.fasta";done >align.sh


ls -l *val_1* |awk '{print $11}' |sed  's|/|\t|g' |cut -f 8 |sed 's/\.gz/_sorted.bam/g' |tr "\n" " " |awk '{print "samtools merge AllRNASEQ_sorted.bam "$0}' |less
 samtools merge AllRNASEQ_sorted.bam 1703FL-02-01_S1_L001_R1_001_val_1.fq_sorted.bam 1703FL-02-02_S2_L001_R1_001_val_1.fq_sorted.bam 1703FL-02-03_S3_L001_R1_001_val_1.fq_sorted.bam 1703FL-02-04_S4_L001_R1_001_val_1.fq_sorted.bam 1703-TM101_S0_L003_R1_001_val_1.fq_sorted.bam 1703-TM102_S0_L003_R1_001_val_1.fq_sorted.bam SRR6230579_1_val_1.fq_sorted.bam SRR6230580_1_val_1.fq_sorted.bam SRR6230581_1_val_1.fq_sorted.bam SRR6230582_1_val_1.fq_sorted.bam SRR6230583_1_val_1.fq_sorted.bam SRR6230584_1_val_1.fq_sorted.bam SRR6230585_1_val_1.fq_sorted.bam SRR6230586_1_val_1.fq_sorted.bam SRR6230587_1_val_1.fq_sorted.bam
```

### align isoseq and ESTs
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA

ln -s /work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/6_isoseq/consensus_isoforms.fasta
ln -s /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/09_Maker/H.glycinesEST.fasta
sh runGmap.sh SCNgenome /work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA SCNgenome.fasta consensus_isoforms.fasta
sh runGmap.sh SCNgenome /work/GIF/remkv6/Baum/04_Dovetail2Restart/11_AlignRNA SCNgenome.fasta H.glycinesEST.fasta
```


### filter the rRNA out of the bam file, as it is excessive in gland reads
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

### set up trinity
```
#Trinity will not run in a satisfactory time with the rDNA reads, so running with them removed.
ln -s ../11_AlignRNA/SCNgenome.fasta
ln -s ../20_ribosomalArrayID/NonRiboReads.bam

sh runTrinity.sh NonRiboReads.bam


#runTrinity.sh
################################################################################
#!/bin/bash

module load GIF2/trinity

bam="$1"
out=$(basename ${bam%.*} |cut -f 1 -d "_")
Trinity \
   --genome_guided_bam ${bam} \
   --max_memory 110G \
   --genome_guided_max_intron 30000 \
   --full_cleanup \
--CPU 16
################################################################################

#finished in 32hrs
#how many transcripts did we get?
grep -c ">" Trinity-GG.fasta
110683
```

### set up braker on unmasked genome
```
#note this was ran with all rRNA included
ln -s ../11_AlignRNA/AllRNASEQ_sorted.bam
ln -s ../11_AlignRNA/SCNgenome.fasta
ln -s ../11_AlignRNA/SCNgenome.consensus_isoforms_sorted.bam
ln -s ../11_AlignRNA/SCNgenome.H.glycinesEST_sorted.bam

module use /work/GIF/software/modules
module load GIF/braker/2.1.0
braker.pl --species=Hglycines2 --genome=SCNgenome.fasta --bam=SCNgenome.consensus_isoforms_sorted.bam,SCNgenome.H.glycinesEST_sorted.bam,AllRNASEQ_sorted.bam

#how many genes?
grep -v "#" augustus.hints.gff |awk '$3=="gene"' |wc
  35514  319626 2025377

```

### run class2 to generate transcripts
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/15_class2
ln -s ../11_AlignRNA/AllRNASEQ_sorted.bam
ln -s ../10_RepeatModeler/SCNgenome.fasta

module load GIF/class2
run_class.pl -a AllRNASEQ_sorted.bam -o AllRNASEQClass2.gtf -p 16 --verbose --clean

#how many transcripts did we generate?
less transcripts.gff |awk '$3=="transcript"' |wc
  60967  853538 8375727
```

### set up stringtie
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/16_Stringtie
#Tried to run class2 with all of the rnaseq, but times out without finishing, so used rrna removed bam.


ln -s ../10_RepeatModeler/SCNgenome.fasta
ln -s ../20_ribosomalArrayID/NonRiboReads.bam

module load stringtie
stringtie NonRiboReads.bam -j 5 -p 16 -v -o NonRrnaRNASEQ_stringtie.gtf

#how many transcripts did we generate?
less NonRrnaRNASEQ_stringtie.gtf |awk '$3=="transcript"' |wc
 41504  747072 6605478
```

### Set up portcullis to analyze splice junctions
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/17_Portcullis
#This initially on the bam file with the rdna reads. Even after removing the rRNA reads, I had to run this on a high mem node, used ~180GB ram.

ln -s ../10_RepeatModeler/SCNgenome.fasta
ln -s ../20_ribosomalArrayID/NonRiboReads.bam

module load portcullis
portcullis full --threads 9 --verbose --use_csi --output portcullis_out --orientation FR SCNgenome.fasta NonRiboReads.bam
```

### Set up braker on a masked genome, excluding the simple repeats
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/19_brakerMasked
ln -s ../11_AlignRNA/AllRNASEQ_sorted.bam
ln -s ../13_braker/SCNgenome.consensus_isoforms_sorted.bam
ln -s ../13_braker/SCNgenome.H.glycinesEST_sorted.bam
ln -s ../10_RepeatModeler/SCNgenome.fasta.masked

module load GIF/braker/2.1.0
braker.pl --species=Hglycines3 --genome=SCNgenome.fasta.masked --bam=SCNgenome.consensus_isoforms_sorted.bam,SCNgenome.H.glycinesEST_sorted.bam,AllRNASEQ_sorted.bam

#how many genes?
awk '$3=="gene"' augustus.hints.gff |grep -v "#" |wc
  22408  201672 1272920

#translate the gff to fasta for use with spades below
~/common_scripts/gff2fasta.pl ../../../10_RepeatModeler/SCNgenome.fasta augustus.hints.gff 4SPADES
```

### Since cufflinks is deprecated to stringtie, I decided to use another assembler for transcripts
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/22_spadesTranscripts
ln -s ../10_RepeatModeler/SCNgenome.fasta
ln -s ../20_ribosomalArrayID/NonRiboReads.bam
ln -s ../19_brakerMasked/braker/Hglycines3/4SPADES.cdna.fasta

#extract only those reads that aligned to the genome
module load picard/2.17.0-ft5qztz; java -Xmx120G -Xms50G -jar  /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/picard-2.17.0-ft5qztzntoymuxiqt3b6yi6uqcmgzmds/bin/picard.jar SamToFastq I=NonRiboReads.bam FASTQ=NonRiboReads_R1.fq F2=NonRiboReads_R2.fq FU=NonRiboReads_unpaired.fq

#Picard failed on some of the reads, due to a mate being removed in ribosomal areas.  ignored this error.

#run spades with extracted reads and using trusted contigs from the masked braker gene prediction
module load spades; spades.py --trusted-contigs 4SPADES.cdna.fasta --rna -m 120 -t 16 -1 NonRiboReads_R1.fq -2 NonRiboReads_R2.fq -o SpadesAssembly
```
