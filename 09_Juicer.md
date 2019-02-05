# Test juicer on test dataset to see if I can use it on SCN

### Install
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff

wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.8_jcuda.0.8.jar
wget https://github.com/aidenlab/juicer/wiki/data/test.txt.gz


git clone https://github.com/theaidenlab/juicer.git
cd juicer/
ln -s CPU/ scripts
cd scripts/
ln -s juicer_tools.1.7.6_jcuda.0.8.jar juicer_tools.jar
cd common/
ln -s ../juicer_tools.1.7.6_jcuda.0.8.jar juicer_tools.jar
cd ../..
mkdir references
cd references
#download a test reference and the bwa index. great!
wget https://s3.amazonaws.com/juicerawsmirror/opt/juicer/references/Homo_sapiens_assembly19.fasta
wget https://s3.amazonaws.com/juicerawsmirror/opt/juicer/references/Homo_sapiens_assembly19.fasta.amb
wget https://s3.amazonaws.com/juicerawsmirror/opt/juicer/references/Homo_sapiens_assembly19.fasta.ann
wget https://s3.amazonaws.com/juicerawsmirror/opt/juicer/references/Homo_sapiens_assembly19.fasta.bwt
wget https://s3.amazonaws.com/juicerawsmirror/opt/juicer/references/Homo_sapiens_assembly19.fasta.pac
wget https://s3.amazonaws.com/juicerawsmirror/opt/juicer/references/Homo_sapiens_assembly19.fasta.sa
cd ..
mkdir restriction_sites
cd restriction_sites/
wget https://s3.amazonaws.com/juicerawsmirror/opt/juicer/restriction_sites/hg19_MboI.txt
cd ..
mkdir fastq;cd fastq
wget http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/HIC003/fastq/HIC003_S2_L001_R1_001.fastq.gz
wget http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/HIC003/fastq/HIC003_S2_L001_R2_001.fastq.gz
cd ..





Dependencies:

module load bwa
module load gnutls/3.5.13-7a3mvfy

#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr

#Current directory
/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/juicer
bash /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/juicer/scripts/juicer.sh

Need to have GPU to do the last step, HiCCUPs.  
#While I did not get HiCCUPs to run due to the lack of a gpu for testing.  The pipeline appears to have worked for the test data.  I will try this on SCN D2 now.  
```
# Run Juicer on Dovetail2
#SCN Dovetail2 scaffolded with SSPace subreads, gap filled with ccs, polished and gap-filled with pilon
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run

git clone https://github.com/theaidenlab/juicer.git
cd juicer/
ln -s CPU/ scripts
cd scripts/
ln -s juicer_tools.1.7.6_jcuda.0.8.jar juicer_tools.jar
cd common/
ln -s ../juicer_tools.1.7.6_jcuda.0.8.jar juicer_tools.jar
cd ../..
mkdir references
cd references
#this is the ref as well as some previously made bwa indices
for f in ../../../juicer/MisAssFixed.Pilon.fast*; do ln -s $f;done
mkdir restriction_sites
cd restriction_sites


```
### Generate MboI restriction cutting sites
#In retrospect, I don't think this was necessary.
```
module load python/2.7.15-ief5zfp
../misc/generate_site_positions.py MboI MisAssFixed.Pilon.fasta /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/restriction_sites/MisAssFixed.Pilon.fasta
#created
MisAssFixed.Pilon.fasta_MboI.txt
```
### Continue Juicer
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/restriction_sites
cd ..
mkdir fastq;cd fastq
for f in /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/02_DovetailFastq/CP4476_chicago_hiseq/*gz ; do ln -s $f;done


#create script to run juicer
module load bwa
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
cd /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer
bash scripts/juicer.sh  -z references/MisAssFixed.Pilon.fasta -p chrom.sizes

```
### parallelize dedup
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/01_split
split -n 16 ../merged_sort.txt
for f in *; do mkdir $f.dir;done
for f in *; do mv $f $f.dir;done
mv ../merged_sort.txt .

#took xaa and xab, changed the name to merged_sort.txt and then moved the file back to aligned/.  Then ran the default juicer script
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes -t 16
#once done, I moved the merged_nodups.txt back to their respective folder
#This can be done in parallel if I recreate the juicer directory for each file.  

cd ../xac.dir.dir/
mv xac merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &
###Did not finish in 96h
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/01_split/xac.dir.dir/aligned
rm opt_dups.txt
rm merged_nodups.txt
rm dups.txt
less merged_sort.txt |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../../01_RemoveRepeats/Centromeres.bed |cut -f 4- >BlacklistedMerged_sort.txt
mv merged_sort.txt PreCentDupMergedSort.txt
cat header BlacklistedMerged_sort.txt >merged_sort.txt
cd ..
module load bwa
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &


cd ../xad.dir.dir/
mv xad merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &
###Did not finish in 96h
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/01_split/xad.dir.dir/aligned
rm optdups.txt
rm merged_nodups.txt
rm dups.txt
less merged_sort.txt |awk 'NR>2' |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../../01_RemoveRepeats/Centromeres.bed |cut -f 4- >BlacklistedMerged_sort.txt
mv merged_sort.txt PreCentDupMergedSort.txt
cat header BlacklistedMerged_sort.txt >merged_sort.txt
cd ..
module load bwa
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xae.dir.dir/
mv xae merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &


cd ../xaf.dir.dir/
mv xaf merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &


cd ../xag.dir.dir/
mv xag merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xah.dir.dir/
mv xah merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xai.dir.dir/
mv xai merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &
###Did not finish in 96h
rm optdups.txt
rm merged_nodups.txt
rm dups.txt
less merged_sort.txt |awk 'NR>2' |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../../01_RemoveRepeats/Centromeres.bed |cut -f 4- >BlacklistedMerged_sort.txt
mv merged_sort.txt PreCentDupMergedSort.txt
cat header BlacklistedMerged_sort.txt >merged_sort.txt
cd ..
module load bwa
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &


cd ../xaj.dir.dir/
mv xaj merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xak.dir.dir/
mv xak merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &
###Did not finish in 96h
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/01_split/xak.dir.dir/aligned/
rm optdups.txt
rm merged_nodups.txt
rm dups.txt
less merged_sort.txt |awk 'NR>2' |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../../01_RemoveRepeats/Centromeres.bed |cut -f 4- >BlacklistedMerged_sort.txt
mv merged_sort.txt PreCentDupMergedSort.txt
cat header BlacklistedMerged_sort.txt >merged_sort.txt
cd ..
module load bwa
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &


cd ../xal.dir.dir/
mv xal merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xam.dir.dir/
mv xam merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xan.dir.dir/
mv xan merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xao.dir.dir/
mv xao merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xap.dir.dir/
mv xap merged_sort.txt
cp ../../../* .
ln -s ../../../fastq/
ln -s ../../../misc/
ln -s ../../../references/
ln -s ../../../restriction_sites/
ln -s ../../../scripts/
ln -s ../../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

```

### DEDUP never finishes.  Trying to create a blacklist of sequences.  Centromeres, telomeres, and tandem satellites are typically included here.
#### Split approach

Awk is not running at 16 cores, so trying to parallelize each to get rid of most before merging into a merged_sort.txt file.
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/01_split/01_RemoveRepeats



###TELOMERE
module load blast-plus/2.7.1-py2-vvbzyor
ln -s ../../../references/MisAssFixed.Pilon.fasta
ln -s /work/GIF/remkv6/Baum/04_Dovetail2Restart/06_Centromere/CentRepeatRead.fast
vi Telomere.fasta
#####################################################################################################
>Telomere
tttagggtttagggtttagggtttagggtttagggtttagggtttagggtttagggtttaggcttaggcttaggctttagggtttagggtttaggg
#####################################################################################################

makeblastdb -in MisAssFixed.Pilon.fasta -dbtype nucl -out Genome.DB

blastn -db Genome.DB -dust no -num_threads 14 -outfmt 6 -query CentRepeatRead.fasta -evalue 100 -num_alignments 100000  -out Centromere2Genome.blastout

blastn -db Genome.DB -dust no -num_threads 14 -outfmt 6 -query Telomere.fasta -evalue 10000 -num_alignments 10000 -word_size 5 -task blastn-short -out Telomere2Genome.blastout

#Could not find a good set of telomeres, so ignoring them.

###CENTROMERE



less Centromere2Genome.blastout |awk '{if($10>$9){print $2,$9,$10} else {print $2,$10,$9}}' |tr " " "\t" |sort -k1,1V -k2,3n |uniq >CentromereUnmerge.bed ;bedtools merge -d 1000 -i CentromereUnmerge.bed >Centromeres.bed


```


### remove centromere reads from merged_nodups.txt, This is faster in cpu time and wait time, but more labor intensive
```
#These two splits were done separately, as they have a different file structure
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/01_split/xaa.dir.dir
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/01_split/xab.dir.dir
less merged_nodups.txt |awk 'NR>2'|awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../01_RemoveRepeats/Centromeres.bed |cut -f 4- >BlacklistedMerged_nodups.txt
mv merged_nodups.txt PreCentDupmergedNodups.txt
mv BlacklistedMerged_nodups.txt merged_nodups.txt


###For the rest of the merged_nodups.txt files that were ran with the softlinked juicer file structure.  split badly split the reads lines, but salvagable. Only losing a ~14x reads x 2 ends =28 reads.   

#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/01_split/xac.dir.dir/aligned

less merged_sort.txt |awk 'NR>1' |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../../01_RemoveRepeats/Centromeres.bed |cut -f 4- >BlacklistedMerged_sort.txt
mv merged_nodups.txt PreCentDupmergedNodups.txt
mv BlacklistedMerged_nodups.txt merged_nodups.txt


#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/01_split/xae.dir.dir/aligned


```
