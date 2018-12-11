# Need to determine if contamination is present in the Dovetail2 assembly, run blobtools


### Setting up megablast
```
/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/19_Dovetail2Blobtools/01_MegaBlast

sh runMegablast.sh nematode_sp_19Jul2018_IbtP1.fasta

runMegablast.sh
#################################################################################
#!/bin/bash
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
#tar -zxvf taxdb.tar.gz

module load GIF/ncbi-blast

FASTA="$1"
blastn \
-task megablast \
-query ${FASTA} \
-db /work/GIF/GIF3/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/nt/nt \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-culling_limit 5 \
-num_threads 12 \
-evalue 1e-5 \
-out ${FASTA%.**}.vs.nt.cul5.1e25.megablast.out

#################################################################################

#This finished before the read mapping, so checking for contamination manually

#how many top hits are not to tylenchida species?
less nematode_sp_19Jul2018_IbtP1.vs.nt.cul5.1e25.megablast.out |sort -u -k1,1V |awk '$15!="Heterodera"' |grep -v "Globodera" |grep -v "Meloidogyne" |awk '{print $1}' |wc
     45      45    1074


#manually blasted sections of these scaffolds to NCBI nucleotide 8/1/18
Scaffold_26;HRSCAF=95  -- This one had a nematode hit, but the best hits were all over the board. (This is H. glycines. has hits to TN20 clones)
Scaffold_37;HRSCAF=129 -- this one did not have a nematode hit.  only tick and moth (looks like a satellite, no decent hits to anything)
Scaffold_181;HRSCAF=321 -- Saprolegnia parasitica scaffold --200bp hit.  No nematode hits (has hits to TN20 clones and other nematodes)
Scaffold_336;HRSCAF=511 -- This had hits to plants and fungi.  no nematodes.  appears to be a polyubiquitin (has hits to tn20 clones)
Scaffold_360;HRSCAF=541 -- no nematode hits, only fish (has hits to TN20 clones)
Scaffold_404;HRSCAF=589 -- Closest hit was to a leech, no nematodes (hits to Rhabditophanes, a nematode)
Scaffold_437;HRSCAF=689  --  Magnaporthe oryzae 70-15 signal recognition particle protein, fungus, no nematodes (has hits to Nippostrongylus brasiliensi, a nematode)
```

#### Map pacbio reads to genome
```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/19_Dovetail2Blobtools/02_PacBio_Bam

ln -s ../../../CamTechGenomeComparison/36_deduplicate122016/subreadMapping/SCN.all.subreads.sl.fasta
ln -s ../nematode_sp_19Jul2018_IbtP1.fasta

module load singularity
singularity build blasr.simg docker://dnadocker/blasr

 printf "module load singularity \n singularity exec blasr.simg blasr SCN.all.subreads.sl.fasta nematode_sp_19Jul2018_IbtP1.fasta --nproc 12 --nCandidates 1   -sam --out subreads2genome.out --bestn 1 --unaligned  subreads2genomeUnaligned.fasta\n" >blasr.sh

singularity exec blasr.simg blasr SCN.all.subreads.sl.fasta nematode_sp_19Jul2018_IbtP1.fasta --nproc 12 --nCandidates 1   -sam --out subreads2genome.out --bestn 1 --unaligned  subreads2genomeUnaligned.fasta


#What kind of mapping percentages did we get?
grep -c ">" subreads2genomeUnaligned.fasta
150719
grep -c ">" SCN.all.subreads.sl.fasta
2382864


#convert to sorted bam
samtools view --threads 16 -b -o subreads2genome.bam subreads2genome.out
samtools sort -m 7G -o subreads2genome_sorted.bam -T  subreads2genome_temp --threads 16  subreads2genome.bam
```
### Run blobtools
```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/19_Dovetail2Blobtools/03_blobtools

ln -s ../nematode_sp_19Jul2018_IbtP1.fasta
ln -s  ../01_MegaBlast/nematode_sp_19Jul2018_IbtP1.vs.nt.cul5.1e25.megablast.out
ln -s ../02_PacBio_Bam/subreads2genome_sorted.bam

sh runBlobtools.sh subreads2genome_sorted.bam nematode_sp_19Jul2018_IbtP1.fasta nematode_sp_19Jul2018_IbtP1.vs.nt.cul5.1e25.megablast.out
```
