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

#map trinity transcripts back to genome
sh runGmap.sh SCNgenome /work/GIF/remkv6/Baum/04_Dovetail2Restart/12_Trinity/trinity_out_dir/ ../SCNgenome.fasta Trinity-GG.fasta

How many gene entries did we get in the gff?
less SCNgenome.Trinity-GG.gff3 |awk '$3=="gene"' |wc
 188834 1699506 21760750


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

### Spades Assembly
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/22_spadesTranscripts

SCNgenome.fasta -> ../10_RepeatModeler/SCNgenome.fasta
NonRiboReads.bam -> ../20_ribosomalArrayID/NonRiboReads.bam
4SPADES.cdna.fasta -> ../19_brakerMasked/braker/Hglycines3/4SPADES.cdna.fasta

#used the cdna's predicted from braker using isoseq, est's, and rnaseq with a masked genome.
SPAdes-3.13.1-Linux/bin/spades.py --trusted-contigs 4SPADES.cdna.fasta --rna -m 300 -t 32 -1 NonRiboReads_R1.fq -2 NonRiboReads_R2.fq -o SpadesAssemblyFat

#/work/GIF/remkv6/Baum/04_Dovetail2Restart/22_spadesTranscripts/SpadesAssemblyFat
grep -c ">" transcripts.fasta
100402

#need to align these so they are available as a gff for mikado


```

## Mikado Setup
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/23_Mikado
To get blast annotations I downloaded all Tylenchina sequences from Trembl (179,064 sequences), as there were only ~5000 that have been manually reviewed for Nematoda in uniprot.
cat uniprot-reviewed\:* >uniprotCombined.fasta

ln -s /work/GIF/remkv6/Baum/04_Dovetail2Restart/12_Trinity/trinity_out_dir/SCNgenome.Trinity-GG.gff3
ln -s /work/GIF/remkv6/Baum/04_Dovetail2Restart/15_class2/AllRNASEQClass2.gtf
ln -s /work/GIF/remkv6/Baum/04_Dovetail2Restart/16_Stringtie/NonRrnaRNASEQ_stringtie.gtf
ln -s /work/GIF/remkv6/Baum/04_Dovetail2Restart/17_Portcullis/portcullis_out/3-filt/portcullis_filtered.pass.junctions.bed
ln -s /work/GIF/remkv6/Baum/04_Dovetail2Restart/22_spadesTranscripts/SpadesAssemblyFat/01_GmapAlign/SCNgenome.transcripts.gff3
ln -s ../10_RepeatModeler/SCNgenome.fasta
ln -s ../15_class2/AllRNASEQ_sorted.bam

vi list.txt
###############################################################################
AllRNASEQClass2.gtf     cl      True
NonRrnaRNASEQ_stringtie.gtf     st      True
SCNgenome.transcripts.gff3      sp      False
SCNgenome.Trinity-GG.gff3       tr      False
###############################################################################


###############################################################################
#!/bin/bash
#SBATCH -A its-hpc-condo-las-free
#SBATCH -N 1
#SBATCH -p freecompute
#SBATCH --ntasks-per-node=16
#SBATCH -t 96:00:00
#SBATCH -J MikadoScript_0
#SBATCH -o MikadoScript_0.o%j
#SBATCH -e MikadoScript_0.e%j
#SBATCH --mail-user=remkv6@istate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
ulimit -s unlimited
cd /work/GIF/remkv6/Baum/04_Dovetail2Restart/23_Mikado
conda activate mikado
#!/bin/bash
#setup variables
genome=SCNgenome.fasta
bam="AllRNASEQ_sorted.bam"
list="list.txt"
#run splice junction prediction
junctions="portcullis_filtered.pass.junctions.bed"
#configure
mikado configure \
   --list $list \
   --reference $genome \
   --mode permissive \
   --scoring worm.yaml \
   --copy-scoring worm.yaml \
   --junctions $junctions configuration.yaml
#prepare
mikado prepare \
   --json-conf configuration.yaml
#blast db
makeblastdb \
   -in uniprotCombined.fasta \
   -dbtype prot \
   -parse_seqids
#blast
blastx \
   -max_target_seqs 5 \
   -num_threads 16 \
   -query mikado_prepared.fasta \
   -outfmt 5 \
   -db uniprotCombined.fasta \
   -evalue 0.000001 2> blast.log | sed '/^$/d' > mikado.blast.xml
blastxml=mikado.blast.xml
#transdecoder
TransDecoder.LongOrfs \
   -t mikado_prepared.fasta
TransDecoder.Predict \
   -t mikado_prepared.fasta \
   --cpu 16
orfs=$(find $(pwd) -name "mikado_prepared.fasta.transdecoder.bed")
#serialise
mikado serialise \
   --start-method spawn \
   --procs 16 \
   --blast_targets ${genome} \
   --json-conf configuration.yaml \
   --xml ${blastxml} \
   --orfs ${orfs}
#pick
mikado pick \
   --start-method spawn \
   --procs 16 \
   --json-conf configuration.yaml \
   --subloci_out mikado.subloci.gff3
scontrol show job $SLURM_JOB_ID
###############################################################################
```

### Second mikado run
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/23_Mikado/02_Round02

ln -s ../SCNgenome.fasta
ln -s ../SCNgenome.DovetailSCNMaker4.all.maker.transcripts.gff3
ln -s ../brakermasked.gff
ln -s ../uniprotCombined.fasta
ln -s ../portcullis_filtered.pass.junctions.bed

vi list2.txt
###############################################################################
SCNgenome.DovetailSCNMaker4.all.maker.transcripts.gff3  MA      True
brakermasked.gff        BR      True
../01_Round1/mikado.loci.gff3   MI      True
###############################################################################

###############################################################################
cd /work/GIF/remkv6/Baum/04_Dovetail2Restart/23_Mikado/02_Round02
conda init bash
conda activate mikado
#!/bin/bash
#setup variables
genome=SCNgenome.fasta
bam="AllRNASEQ_sorted.bam"
list="list2.txt"
#run splice junction prediction
junctions="portcullis_filtered.pass.junctions.bed"
#configure
mikado configure \
   --list $list \
   --reference $genome \
   --mode stringent \
   --scoring worm.yaml \
   --copy-scoring worm.yaml \
   --junctions $junctions configuration.yaml
#prepare
mikado prepare \
   --json-conf configuration.yaml
#blast db
#makeblastdb \
   -in uniprotCombined.fasta \
   -dbtype prot \
   -parse_seqids
#blast
blastx \
   -max_target_seqs 5 \
   -num_threads 16 \
   -query mikado_prepared.fasta \
   -outfmt 5 \
   -db uniprotCombined.fasta \
   -evalue 0.000001 2> blast.log | sed '/^$/d' > mikado.blast.xml
blastxml=mikado.blast.xml
#transdecoder
TransDecoder.LongOrfs \
   -t mikado_prepared.fasta
TransDecoder.Predict \
   -t mikado_prepared.fasta \
   --cpu 16
orfs=$(find $(pwd) -name "mikado_prepared.fasta.transdecoder.bed")
#serialise
mikado serialise \
   --start-method spawn \
   --procs 16 \
   --blast_targets ${genome} \
   --json-conf configuration.yaml \
   --xml ${blastxml} \
   --orfs ${orfs}
#pick
mikado pick \
   --start-method spawn \
   --procs 16 \
   --json-conf configuration.yaml \
   --subloci_out mikado.subloci.gff3
###############################################################################
```

## Compare output with existing annotations
### Schachtii gene call stats
```
#gene length
less ../../../12_SchachtiiSynteny/H_sch_gene_calls_v1_CP.gff |awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  78,025,824
Count:  26,739
Mean:   2,918
Median: 2,129
Min:    252
Max:    56,948

#CDS length
(mikado) [remkv6@condofree032 02_Round02]$ less ../../../12_SchachtiiSynteny/H_sch_gene_calls_v1_CP.gff |awk '$3=="CDS"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  36,075,582
Count:  222,897
Mean:   161
Median: 125
Min:    1
Max:    13,080

#transcript length
less ../../../12_SchachtiiSynteny/H_sch_gene_calls_v1_CP.gff |awk '$3=="transcript"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  101,553,992
Count:  32,624
Mean:   3,112
Median: 2,275
Min:    252
Max:    56,948

#exons per transcript
less ../../../12_SchachtiiSynteny/H_sch_gene_calls_v1_CP.gff |awk '$3=="CDS" ' |cut -f 9 |awk '{print $1}' |sort |uniq -c |awk '{print $1}' |summary.sh
Total:  222,897
Count:  32,624
Mean:   6
Median: 5
Min:    1
Max:    113


```

### Mikado.loci.gff3 round2, including repeat genes
```
#gene length
less mikado.loci.gff3 |awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  120,922,982
Count:  38,646
Mean:   3,128
Median: 821
Min:    199
Max:    465,216

#transcript length
less mikado.loci.gff3 |awk '$3=="mRNA"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  142,902,933
Count:  45,286
Mean:   3,155
Median: 1,034
Min:    199
Max:    465,216

#exons per transcript
less mikado.loci.gff3 |awk '$3=="CDS" ' |cut -f 9 |sed 's/\./\t/3' |awk '{print $1}' |sort|uniq -c |awk '{print $1}' |summary.sh
Total:  258,802
Count:  45,286
Mean:   5
Median: 3
Min:    1
Max:    194

```

### mikado.loci.gff3
```
#exons per transcript
less mikado.loci.gff3 |awk '$3=="CDS"' |bedtools intersect -v -wo -a - -b ../../10_RepeatModeler/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/\./\t/3' |awk '{print $1}' |sort|uniq -c |awk '{print $1}' |summary.sh
Total:  228,567
Count:  31,566
Mean:   7
Median: 5
Min:    1
Max:    194

#how many genes are not repetitive
less mikado.loci.gff3 |awk '$3=="CDS"' |bedtools intersect -v -wo -a - -b ../../10_RepeatModeler/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/\./\t/2' |awk '{print $1}' |sort|uniq|wc
25180   25180  650224

#how many are repetitive then?
By subtraction that is: 13,466 genes

#gene sizes that are not repetitive
less NonRepetitiveGenes.gff3 |awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  104,488,696
Count:  25,180
Mean:   4,149
Median: 1,726
Min:    199
Max:    465,216

```
### gene stats of braker masked annotation
```
less augustus.hints.gff |awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  67,011,405
Count:  22,408
Mean:   2,990
Median: 2,015
Min:    200
Max:    52,212

```

### gene stats of braker unmasked annotation
```
less augustus.hints.gff |awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  83,144,932
Count:  35,514
Mean:   2,341
Median: 1,546
Min:    200
Max:    43,263
```

### gene stats of maker annotation on old dovetail genome
```
/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding

less 12_MakerGenesOrthofinder/DovetailSCNMaker4.all.NOFASTA.gff|awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  67,168,685
Count:  22,856
Mean:   2,938
Median: 1,811
Min:    5
Max:    65,204

```


## Round2 of mikado using all ESTs from tylenchida
```
less mikado.loci.gff3| awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  131,514,741
Count:  39,516
Mean:   3,328
Median: 920
Min:    199
Max:    465,216
```

### gene calls from 738 genome
```

less ../CamTechGenomeComparison/58_Renamatorium/1_genomeNgff/fixed.augustus.gff3|awk '$3=="gene"' |grep -v "#"|awk '{if($4>$5){print $4-$5} else {print $5-$4}}' |summary.sh
Total:  71,832,760
Count:  29,769
Mean:   2,413
Median: 1,603
Min:    79
Max:    65,717

```

### Functional annotation stats
```
#interproscan
less 01_Interpro/interproAnnot.tsv |awk '{print $1}' |sort|uniq|wc
  25779   25779  639842

#proteins to uniprot
less 04_ProtsUniprot/mikado_proteins.vs.uniprot_sprot.cul5.1e5.blastp.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1 }' |sort|uniq|wc
   1604    1604   39727

#transcripts to uniref
less 05_TransUniprot/mikado_transcripts.vs.uniprot_sprot.cul5.1e5.blastx.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1 }' |sort|uniq|wc
  12796   12796  317987

#prots to nr
less 02_Prots2Nr/mikado_proteinsFixed.vs.
nr.cul5.1e5.blastp.out |grep -v "hypothetical" |grep -v "uncharacterized" |aw
k '{print $1}' |sort|uniq|wc
   3056    3056   75814

#transcripts to nt
ess 03_Transcrips2Nt/mikado_transcripts.vs.nt.cul5.1e5.blastn.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1}' |sort|uniq|wc
   2266    2266   56298



#All databases together
cat <(less 01_Interpro/interproAnnot.tsv |awk '{print $1}' |sort|uniq) <(less 02_Prots2Nr/mikado_proteinsFixed.vs.nr.cul5.1e5.blastp.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1}' |sort|uniq) <( less 03_Transcrips2Nt/mikado_transcripts.vs.nt.cul5.1e5.blastn.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1}' |sort|uniq) <(less 04_ProtsUniprot/mikado_proteins.vs.uniprot_sprot.cul5.1e5.blastp.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1 }' |sort|uniq) <(less 05_TransUniprot/mikado_transcripts.vs.uniprot_sprot.cul5.1e5.blastx.out |grep -v "hypothetical" |grep -v "uncharacterized" |awk '{print $1 }' |sort|uniq) |sort|uniq|wc
  26951   26951  668844

````

##  Filter Mikado Genes by expression and Repeats to obtain a final count
```
# Every genes name
awk '$3=="gene"' mikado.loci.gff3 |cut -f 9 |sort|uniq|sed 's/ID=//g' |sed 's/;/\t/1' |cut -f 1 >AllGenes.list
wc AllGenes.list
 39516  39516 900462 AllGenes.list

# Genes that have a 50% overlap with a repeat from EDTA
bedtools intersect -f .5 -wo -a <(awk '$3=="gene"' mikado.loci.gff3) -b ../33_EDTA/EDTA/SCNgenome.fasta.EDTA.TEanno.gff |cut -f 9 |sort|uniq|sed 's/ID=//g' |sed 's/;/\t/1' |cut -f 1 >RepeatGenes.list
wc RepeatGenes.list
8316   8316 189033 RepeatGenes.list

less ../38_Expression/GeneCounts |awk '$7<2' |cut -f 1 |sed 's/\./\t/2' |cut -f 1 |sort|uniq >GenesNotExpressed.list

 wc GenesNotExpressed.list
 19278  19280 439092 GenesNotExpressed.list

#How many genes are expressed and have less than 50% repeat overlap with the gene space?
 cat AllGenes.list RepeatGenes.list GenesNotExpressed.list |sort|uniq -c |awk '$1==1 {print $2}'  >NonRepeatExpressedGenes.list
wc NonRepeatExpressedGenes.list
 17431 17431 34862 NonRepeatExpressedGenes.list

#How many of these 17,431 genes have a functional annotation?
less 06_Combine/SCNgenomeFunctionalGeneAnnotations.gff3 |awk '$3=="mRNA"' |sed 's/;/\t/1' |sed 's/ID=//g' |cut -f 9- |sed 's/\./\t/2' |grep -w -f NonRepeatExpressedGenes.list  - |grep -c "Note"
 15968

# Separate high confidence genes from others (repetitive and not expressed)
less 06_Combine/SCNgenomeFunctionalGeneAnnotations.gff3 |awk '$3=="gene"' |sed 's/ID=//g' |sed 's/;/\t;/g' |grep -w -f NonRepeatExpressedGenes.list - |awk '{print $1,$2,$3,$4,$5,$6,$7,$8,"ID="$9$10$11}' |tr " " "\t" >SCNgenomeFunctionalGeneAnnotationsHighConfidenceGenesOnly.gff3 &

awk '$3=="mRNA"' mikado.loci.gff3 |sed 's/ID=//g' |sed 's/;/\t;/1' |cut -f 9 |sed 's/\./\t\./2' >AllmRNA.mikado.loci.list





less 06_Combine/SCNgenomeFunctionalGeneAnnotations.gff3  |sed 's/\.CDS/\t\.CDS/1' |sed 's/\.exon/\t\.exon/1' |sed 's/\.three/\t\.three/1' |sed 's/;/\t/1' |sed 's/\.five/\t\.five/1' |sed 's/ID=/ID=\t/1' >mikado.loci.Grepmod.gff3


less NonRepeatExpressedmRNAs.list |while read line; do echo "awk '\$10==\""$line"\"' mikado.loci.Grepmod.gff3 >>SCNgenomeFunctionalGeneAnnotationsHighConfidenceRemainingFeatures.gff3" ;done >GetRemainingFeatures.sh
 sh GetRemainingFeatures.sh &

 less SCNgenomeFunctionalGeneAnnotationsHighConfidenceRemainingFeatures.gff3 |awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9$10$11";"$12$13}' |tr " " "\t" |cat - SCNgenomeFunctionalGeneAnnotationsHighConfidenceGenesOnly.gff3 |sort -k1,1V -k4,5n >SCNgenomeFuctionalGeneAnnotationsHighConfidenceAllGenesUnsorted.gff3

perl gff3sort/gff3sort.pl --precise --chr_order natural SCNgenomeFuctionalGeneAnnotationsHighConfidenceAllGenesUnsorted.gff3 >SCNgenomeFunctionalGeneAnnotationsHighConfidenceAllFeatures.gff3


```



# Gene models were found to not have start codons about 50% of the time, Anju says the splicing is not right in some of the effectors.  

Need to troublshoot.
###  Copy data, create folders, install

```
#Move data from condo to nova,  
from condo
 /work/GIF/remkv6/Baum/04_Dovetail2Restart/23_Mikado/04_Round2Redo
to nova
 /work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn

#had to delete all old mikado files for the rerun.








#kept getting this error
File "/home/remkv6/.conda/envs/mikado/lib/python3.6/site-packages/Mikado/serializers/blast_serializer/xml_serialiser.py", line 161, in run
  for pickled in self._pickler(filename):
File "/home/remkv6/.conda/envs/mikado/lib/python3.6/site-packages/Mikado/serializers/blast_serializer/xml_serialiser.py", line 105, in _pickler
  max_target_seqs=self.__max_target_seqs)
File "/home/remkv6/.conda/envs/mikado/lib/python3.6/site-packages/Mikado/serializers/blast_serializer/xml_serialiser.py", line 817, in objectify_record
  current_target = _get_target_for_blast(self, alignment)
File "/home/remkv6/.conda/envs/mikado/lib/python3.6/site-packages/Mikado/serializers/blast_serializer/xml_serialiser.py", line 780, in _get_target_for_blast
  raise KeyError("{} not found in the targets!".format(alignment.accession))
KeyError: '47529 not found in the targets!'

#tried lots of attempts
  checked correctness of list2.txt
  47529 appears to be a blast target fom the Tylenchida EST sequences.  
  serialize.log says some braker genes failed to be indexed
  tried adding --trancripts mikado_prepared.fasta to the mikado file, but did not have an effect
  redownloaded original augustus.hints.gtf from braker masked and converted to GFF3 via genometools gtf_2_gff3
ml genometools
gt gtf_to_gff3 -tidy -o augustus.hints.gff3 augustus.hints.gtf
version 1.5.9
  Checked the overlap between portcullis and the other transcript files.  It only seems to match up well with the braker predictions and the Round1mikado.loci.gff3.
bedtools intersect -wo -a <(awk '$3=="CDS"' SCNgenome.DovetailSCNMaker4.all.maker.transcripts.gff3)  -b portcullis_filtered.pass.junctions.bed |cut -f 1,4,5,10,16,17 |less
bedtools intersect -wo -a <(awk '$3=="CDS"' augustus.hints.gff3)  -b portcullis_filtered.pass.junctions.bed |cut -f 1,4,5,10,16,17 |less
bedtools intersect -wo -a <(awk '$3=="CDS"' Round1mikado.loci.gff3)  -b portcullis_filtered.pass.junctions.bed |cut -f 1,4,5,10,16,17|less  
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


### Mikado script evolution

Running this one step at a time to ensure I catch errors
```

#list3.txt
#################################################################################
SCNgenome.DovetailSCNMaker4.all.maker.transcripts.gff3  MA      True    -1
Round1mikado.loci.gff3  MI      True
SCNgenome.TylenchidaESTNotH.glycines.gff3       RE      False   -1
augustus.hints.gff3     BR      True    1
#################################################################################

#MikadoScript_0.sub
################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -A its-hpc-condo-las-free
#SBATCH -p freecompute
#SBATCH -t 96:00:00
#SBATCH -J MikadoScript_0
#SBATCH -o MikadoScript_0.o%j
#SBATCH -e MikadoScript_0.e%j
#SBATCH --mail-user=remkv6@istate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
ulimit -s unlimited
cd /work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn
source activate mikado
#!/bin/bash
#setup variables
genome=SCNgenome.fasta
bam="AllRNASEQ_sorted.bam"
list="list3.txt"
#run splice junction prediction
junctions="portcullis_filtered.pass.junctions.bed"
#configure
#mikado configure \
#   --list $list \
#   --reference $genome \
#   --mode stringent \
#   --scoring worm.yaml \
#   --copy-scoring worm.yaml \
#   --junctions $junctions configuration.yaml
#prepare
#mikado prepare \
#   --json-conf configuration.yaml
#blast db
makeblastdb \
   -in uniprotCombined.fasta \
   -dbtype prot \
   -parse_seqids
#blast
blastx \
   -max_target_seqs 5 \
   -num_threads 35 \
   -query mikado_prepared.fasta \
   -outfmt 5 \
   -db uniprotCombined.fasta \
   -evalue 0.000001 2> blast.log | sed '/^$/d' > mikado.blast.xml
blastxml=mikado.blast.xml
#transdecoder
#TransDecoder.LongOrfs \
#   -t mikado_prepared.fasta
#TransDecoder.Predict \
#   -t mikado_prepared.fasta \
#   --cpu 16
orfs=$(find $(pwd) -name "mikado_prepared.fasta.transdecoder.bed")
#serialise
#mikado serialise \
#   --start-method spawn \
#   --procs 16 \
#   --blast_targets ${genome} \
#   --json-conf configuration.yaml \
#   --xml ${blastxml} \
 #  --orfs ${orfs} \
#   -mr .5

#pick
#mikado pick \
#   --start-method spawn \
#   --procs 16 \
#   --json-conf configuration.yaml \
#   --subloci_out mikado.subloci.gff3 \
#   --pad
scontrol show job $SLURM_JOB_ID
```
### Post analysis of latest mikado run with all but ab inito
```
less mikado.loci.gff3 |sort -k1,1V -k4,5nr |uniq |grep -v "#" |sed 's/transcript/mRNA/1'>UniquedRound2mikado.loci.gff3
gt gff3 -sortlines -checkids -fixregionboundaries -tidy UniquedRound2mikado.loci.gff3 >tidiedUniquedRound2mikado.loci.gff3


```


### Since mikado fails to produce proteins with reliable splicing, start codons, and stop codons, filter all proteins that do not meet that criteria
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn/03_MethionineProteinsOnly


#braker genes, masked?
less ../augustus.hints.proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |wc
 24479   48958 10191962
less ../augustus.hints.proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |tr "\t" "\n" >Metaugustus.hints.proteins.fasta


#maker genes on 368 scaffold genome
less ../SCNgenome.DovetailSCNMaker4.all.maker.transcripts.proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |wc
  55950  111900 18316522
  less ../SCNgenome.DovetailSCNMaker4.all.maker.transcripts.proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |tr "\t" "\n" >MetSCNgenome.DovetailSCNMaker4.all.maker.transcripts.proteins.fasta

#spades transcripts
less ../SCNgenome.transcripts.proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |wc
42561   85122 10350491
less ../SCNgenome.transcripts.proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |tr "\t" "\n" >MetSpadesSCNgenome.transcripts.proteins.fasta

#Trinity transcripts
less ../SCNgenome.Trinity-GG.proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |wc
    72038  144076 14447267
less ../SCNgenome.Trinity-GG.proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |tr "\t" "\n" >MetSCNgenome.Trinity-GG.proteins.fasta

# NCBI EST proteins
less ../SCNgenome.TylenchidaESTNotH.glycines.proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |wc
  11857   23714 1639082
#totally odd that so many did not have start codons
grep -c ">" ../SCNgenome.TylenchidaESTNotH.glycines.proteins.fasta
  60218
less ../SCNgenome.TylenchidaESTNotH.glycines.proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |tr "\t" "\n" >MetSCNgenome.TylenchidaESTNotH.glycines.proteins.fasta

# Class2 transcripts
gt gtf_to_gff3 -tidy <(sort -k1,1V -k4,5nr AllRNASEQClass2.gtf|grep -v "#" |uniq) >AllRNASEQClass2.gff3
gt gff3 -sortlines -tidy -fixregionboundaries  AllRNASEQClass2.gff3 >tidiedAllRNASEQClass2.gff3
gt cds -startcodon -finalstopcodon -seqfile SCNgenome.fasta  -matchdesc tidiedAllRNASEQClass2.gff3 >FixedClass2.gff3
gffread FixedClass2.gff3 -g SCNgenome.fasta -t mRNA -x AllRNASEQClass2._transcripts.fasta -y AllRNASEQClass2._proteins.fasta

less ../AllRNASEQClass2._proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |wc
  36212   72424 12980842
less ../AllRNASEQClass2._proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |tr "\t" "\n" >MetAllRNASEQClass2._proteins.fast


#stringtie conversion

gt gtf_to_gff3 -tidy <(sort -k1,1V -k4,5nr NonRrnaRNASEQ_stringtie.gtf|grep -v "#" |uniq) >NonRrnaRNASEQ_stringtie.gff3
gt gff3 -sortlines -tidy -fixregionboundaries  NonRrnaRNASEQ_stringtie.gff3 >tidiedNonRrnaRNASEQ_stringtie.gff3
gt cds -startcodon -finalstopcodon -seqfile SCNgenome.fasta  -matchdesc tidiedNonRrnaRNASEQ_stringtie.gff3 >FixedNonRrnaRNASEQ_stringtie.gff3
gffread FixedNonRrnaRNASEQ_stringtie.gff3 -g SCNgenome.fasta -t mRNA -x NonRrnaRNASEQ_stringtie_transcripts.fasta -y NonRrnaRNASEQ_stringtie_proteins.fasta


less ../NonRrnaRNASEQ_stringtie_proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |wc
  29063   58126 10216583

less ../NonRrnaRNASEQ_stringtie_proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |tr "\t" "\n" >MetNonRrnaRNASEQ_stringtie_proteins.fasta


#braker unmasked
BrakerUnmasked.gff


gt gtf_to_gff3 -tidy <(sort -k1,1V -k4,5nr BrakerUnmasked.gff|grep -v "#" |uniq) >BrakerUnmasked.gff3
gt gff3 -sortlines -tidy -fixregionboundaries  BrakerUnmasked.gff3 >tidiedBrakerUnmasked.gff3
gt cds -startcodon -finalstopcodon -seqfile SCNgenome.fasta  -matchdesc tidiedBrakerUnmasked.gff3 >FixedBrakerUnmasked.gff3
gffread FixedBrakerUnmasked.gff3 -g SCNgenome.fasta -t mRNA -x BrakerUnmasked_transcripts.fasta -y BrakerUnmasked_proteins.fasta


less ../BrakerUnmasked_proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |wc
37894   75788 15010594

less ../BrakerUnmasked_proteins.fasta |awk '{print $1}' |tr "\n" "\t" |sed 's/>/\n>/g'  |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |awk 'substr($2,1,1)=="M"'  |tr "\t" "\n" >MetBrakerUnmasked_proteins.fasta




for f in *fasta; do echo "mkdir "$f"dir; cd "$f"dir; ln -s ../SCNgenome.fasta; ln -s ../"$f" ; ml miniconda3; source activate ge
nomethreader;gth -genomic SCNgenome.fasta -protein "$f" -gff3out -species nematode -skipalignmentout -o "${f%.*}"aln -force";done  >gth.sh

```
