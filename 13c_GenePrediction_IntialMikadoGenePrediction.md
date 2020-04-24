
# Need to predict genes for the SCN TN10 pseudomolecule assembly


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
