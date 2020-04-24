
# Prepare all input assemblies and predictions for gene prediction

### Trinity
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
