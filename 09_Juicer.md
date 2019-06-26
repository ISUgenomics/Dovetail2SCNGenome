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

less merged_nodups.txt |awk 'NR>1' |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../../01_RemoveRepeats/Centromeres.bed |cut -f 4- >BlacklistedMerged_sort.txt
mv merged_nodups.txt PreCentDupmergedNodups.txt
mv BlacklistedMerged_nodups.txt merged_nodups.txt



cat xaa.dir.dir/header x*.dir.dir/aligned/merged_nodups.txt >merged_nodups.txt
less merged_nodups.txt |awk 'NR>4' |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b 01_RemoveRepeats/Centromeres.bed |cut -f 4- >BlacklistedMerged_nodups.txt
mv merged_nodups.txt PreCentDupmergedNodups.txt
mv BlacklistedMerged_nodups.txt merged_sort.txt



```




# HALTED APPROACH 2/5/19.  Restart here from failed dedup
```
This will remove the centromeric reads from the sequencing.  They just cause massive drag from multiple alignment issues.  
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned

less merged_sort.txt |awk 'NR>1' |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -iobuf 100G -v -a - -b 01_split/01_RemoveRepeats/Centromeres.bed |cut -f 4- >BlacklistedMerged_sort.txt &

# then run this
split -n l/16 BlacklistedMerged_sort.txt
for f in x*; do mkdir $f.dir;done
for f in x*; do mv $f $f.dir;done

cd xaa.dir.dir/
mv xaa merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xab.dir.dir/
mv xab merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xac.dir.dir/
mv xac merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xad.dir.dir/
mv xad merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xae.dir.dir/
mv xae merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &


cd ../xaf.dir.dir/
mv xaf merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &


cd ../xag.dir.dir/
mv xag merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xah.dir.dir/
mv xah merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xai.dir.dir/
mv xai merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xaj.dir.dir/
mv xaj merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xak.dir.dir/
mv xak merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xal.dir.dir/
mv xal merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xam.dir.dir/
mv xam merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xan.dir.dir/
mv xan merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xao.dir.dir/
mv xao merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

cd ../xap.dir.dir/
mv xap merged_sort.txt
cp ../../* .
ln -s ../../fastq/
ln -s ../../misc/
ln -s ../../references/
ln -s ../../restriction_sites/
ln -s ../../scripts/
ln -s ../../splits/
mkdir aligned
cd aligned/
mv ../merged_sort.txt .
cd ../
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &
```

### Looking around in duplicates to see what is causing the problem.
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/xam.dir.dir/aligned
 less dups.txt |awk '{print $2,$3,$6,$7}' |sort --parallel 10 -k1,1V -k2,4n >dupCoords.tst


samtools faidx references/MisAssFixed.Pilon.fasta
 samtools faidx ../xam.dir.dir/references/MisAssFixed.Pilon.fasta "scaffold4|size13779824|pilon:5570000-5610000" > scaffold4_5570000-5610000.fasta
 module load blast-plus/2.7.1-py2-vvbzyor
 makeblastdb -in ../xam.dir.dir/references/MisAssFixed.Pilon.fasta -dbtype nucl -out Genome.blastdb
 blastn -db Genome.blastdb -query scaffold4_5570000-5610000.fasta -outfmt 6 -num_threads 10 -out scaffold4_5570000-5610000Genome.blastout


#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/02_InvestigateDupRegions
 mkdir 02_InvestigateDupRegions
 cd 02_InvestigateDupRegions/
 samtools faidx ../xam.dir.dir/references/MisAssFixed.Pilon.fasta "scaffold4|size13779824|pilon:5570000-5610000" > scaffold4_5570000-5610000.fasta

 #this is a small repeatmodeler repeat region
less scaffold4_5570000-5610000Genome.blastout  |awk '$12>100{if ($10>$9) {print $2"\t"$9"\t"$10} else {print $2"\t"$10"\t"$9}}' |awk '{print $1,$2"#"$1,$3}' |tr "#" "\n" |awk '{print $1,$2,$2+150}' |tr " " "\t" |sort -k1,1V -k2,3n |bedtools  merge -d 1000 -i -  >FilterRepeatscaffold4_5570000-5610000.bed

less merged_sort.txt |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../../02_InvestigateDupRegions/FilterRepeatscaffold4_5570000-5610000.bed |cut -f 4- >BlacklistedMerged_sort.txt
###XAM continued, as I found a better method -- Remove all genomic repeats.




#xan
less dups.txt |awk '{print $2,$3/1000}' |sed 's/\./\t/g' |awk '{print $1,$2}' |uniq -c  |awk '$1>6000' |less

21552 scaffold5|size11380747|pilon 1661
 38209 scaffold5|size11380747|pilon 2405
 11320 scaffold5|size11380747|pilon 3626
  9626 scaffold5|size11380747|pilon 4286
  7908 scaffold5|size11380747|pilon 4287
  9033 scaffold5|size11380747|pilon 4288
 13465 scaffold5|size11380747|pilon 4616
 23058 scaffold5|size11380747|pilon 7541
 31932 scaffold5|size11380747|pilon 7542
 11627 scaffold5|size11380747|pilon 7546
 32657 scaffold5|size11380747|pilon 7547
 34223 scaffold5|size11380747|pilon 7544
 39082 scaffold5|size11380747|pilon 7545
 51931 scaffold5|size11380747|pilon 7546
 19728 scaffold5|size11380747|pilon 7547
  6824 scaffold5|size11380747|pilon 7548
  6081 scaffold5|size11380747|pilon 7551
  7467 scaffold5|size11380747|pilon 7629
618355 scaffold5|size11380747|pilon 9162
208989 scaffold5|size11380747|pilon 9163

#this is buried in tandem repeats
 samtools faidx ../../xam.dir.dir/references/MisAssFixed.Pilon.fasta "scaffold5|size11380747|pilon:9162000-9164000" > ../../02_InvestigateDupRegions/scaffold5:9162000-9164000.fasta
 blastn -db Genome.blastdb -query scaffold5:9162000-9164000.fasta -outfmt 6 -num_threads 10 -out scaffold5:9162000-9164000Genome.blastout

 less scaffold5\:9162000-9164000Genome.blastout  |awk '$12>200{if ($10>$9) {print $2"\t"$9"\t"$10} else {print $2"\t"$10"\t"$9}}' |awk '{print $1,$2"#"$1,$3}' |tr "#" "\n" |awk '{print $1,$2,$2+150}' |tr " " "\t" |sort -k1,1V -k2,3n |bedtools  merge -d 1000 -i -  >FilterTDRscaffold59162000-9164000Genome.bed

 less merged_sort.txt |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../../02_InvestigateDupRegions/FilterTDRscaffold59162000-9164000Genome.bed |cut -f 4- >BlacklistedMerged_sort.txt
 mv merged_sort.txt Oldmerged_sort.txt
 mv BlacklistedMerged_sort.txt merged_sort.txt
 bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &


###XAP
/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/xap.dir.dir/aligned
less dups.txt |awk '{print $2,$3/1000}' |sed 's/\./\t/g' |awk '{print $1,$2}' |uniq -c  |awk '$1>10000' |less
11131 scaffold7|size7744779|pilon 6288
17354 scaffold7|size7744779|pilon 36
10619 scaffold85|size264353|pilon 54
11755 scaffold85|size264353|pilon 55
10079 scaffold85|size264353|pilon 57
10317 scaffold85|size264353|pilon 58
22465 scaffold85|size264353|pilon 59
16531 scaffold85|size264353|pilon 60
24959 scaffold85|size264353|pilon 61
45725 scaffold85|size264353|pilon 62
37467 scaffold85|size264353|pilon 63
28481 scaffold85|size264353|pilon 64
23497 scaffold86|size209102|pilon 52
34269 scaffold86|size209102|pilon 52
10388 scaffold86|size209102|pilon 53
17559 scaffold86|size209102|pilon 61
13577 scaffold86|size209102|pilon 63
18282 scaffold86|size209102|pilon 68
31201 scaffold86|size209102|pilon 70
28023 scaffold86|size209102|pilon 71
23396 scaffold86|size209102|pilon 73
12737 scaffold86|size209102|pilon 52
11347 scaffold86|size209102|pilon 53
22096 scaffold87|size186194|pilon 6
15658 scaffold87|size186194|pilon 31
36698 scaffold87|size186194|pilon 32
187994 scaffold87|size186194|pilon 33
49147 scaffold87|size186194|pilon 31
63996 scaffold87|size186194|pilon 32
180672 scaffold87|size186194|pilon 33

module load blast-plus/2.7.1-py2-vvbzyor
cd ../../02_InvestigateDupRegions/
samtools faidx ../../xam.dir.dir/references/MisAssFixed.Pilon.fasta "scaffold87|size186194|pilon:31000-34000" > ../../02_InvestigateDupRegions/scaffold87-31000-34000.fasta
blastn -db Genome.blastdb -query scaffold87-31000-34000.fasta -outfmt 6 -num_threads 10 -out scaffold87-31000-34000.blastout

less scaffold87-31000-34000.blastout |awk '$12>200{if ($10>$9) {print $2"\t"$9"\t"$10} else {print $2
"\t"$10"\t"$9}}' |awk '{print $1,$2"#"$1,$3}' |tr "#" "\n" |awk '{print $1,$2,$2+150}' |tr " " "\t" |sort -k1,1V -k2,3n |awk 'NR>2' |bedtools  me
rge -d 1000 -i -  >Filterscaffold87-31000-34000.bed

less merged_sort.txt |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../../02_InvestigateDupRegions/Filterscaffold87-31000-34000.bed |cut -f 4- >BlacklistedMerged_sort.txt
mv merged_sort.txt Oldmerged_sort.txt
mv BlacklistedMerged_sort.txt merged_sort.txt
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &




### XAG
less merged_sort.txt |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../../02_InvestigateDupRegions/DedupList.bed |cut -f 4- >BlacklistedMerged_sorted.txt
mv merged_sort.txt OldMergedSort.txt
mv BlacklistedMerged_sorted.txt merged_sort.txt
cd ..
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

### XAD -- never completely finished problems with scaffold 1 at 10.9MB in the 31.4Mb scaffold
less merged_sort.txt |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../../02_InvestigateDupRegions/DedupList.bed |c
ut -f 4- >BlacklistedMerged_sorted.txt
mv merged_sort.txt OldMergedSort.txt
mv BlacklistedMerged_sorted.txt merged_sort.txt
cd ..
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &

###XAC
less merged_sort.txt |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../../02_InvestigateDupRegions/DedupList.bed |c
ut -f 4- >BlacklistedMerged_sorted.txt
mv merged_sort.txt OldMergedSort.txt
mv BlacklistedMerged_sorted.txt merged_sort.txt
cd ..
bash scripts/juicer.sh -S dedup -z references/MisAssFixed.Pilon.fasta -p chrom.sizes &
```
### generate hic file manually with mostly deduped reads, minus jettisoning a few Gb of duplicated regions.  
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned
cat x*/aligned/merged_nodups.txt >merged_nodups.txt

cd /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/02_juicerUnmasked/aligned
module load bwa
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
java -Xmx2g -jar ../scripts/juicer_tools.jar pre merged_nodups.txt merged_nodups2.hic ../chrom.sizes

```

### run 3D DNA pipeline
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/02_juicebox/
git clone https://github.com/theaidenlab/3d-dna.git
cd 3d-dna/

module use /work/GIF/software/modules
module load GIF2/lastz/1.03.73
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
module load python
module load parallel/20170322-36gxsog
bash run-asm-pipeline.sh -e -m diploid /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/references/MisAssFixed.Pilon.fasta /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/merged_nodups.txt

runs in abot 8-9hrs 16cpu.  140gb bam merged_nodups.txt file.
```




#  Set up juicer to run with Hi-C reads

```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/

module load bedtools2
module load bwa
git clone https://github.com/theaidenlab/juicer.git
cd juicer/

mkdir references; cd references
ln -s  /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/references/MisAssFixed.Pilon.fasta
cd ..

ln -s CPU/scripts/ scripts
cd scripts/
wget http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/juicer_tools.1.7.6_jcuda.0.8.jar
ln -s juicer_tools.1.7.6_jcuda.0.8.jar juicer_tools.jar
cd ..
mkdir fastq
cd fastq/
for f in /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/02_DovetailFastq/CP4477_hic_hiseq/*gz ; do ln -s $f;done

cd ../references
ln -s /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/02_InvestigateDupRegions/DedupList.bed
bedtools maskfasta -fi MisAssFixed.Pilon.fasta -fo MaskedMisAssFixed.Pilon.fasta -bed DedupList.bed
unlink MisAssFixed.Pilon.fasta

module load bwa
bwa index MaskedMisAssFixed.Pilon.fasta

cd ../
ln -s ../../01_scnD2run/juicer/chrom.sizes
mkdir restriction_sites
cd restriction_sites
python ../misc/generate_site_positions.py MboI MaskedMisAssFixed.Pilon.fasta /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/juicer/references/MaskedMisAssFixed.Pilon.fasta
cd ..

#RUN JUICER
bash scripts/juicer.sh  -y restriction_sites/MaskedMisAssFixed.Pilon.fasta_MboI.txt -z references/MaskedMisAssFixed.Pilon.fasta -p chrom.sizes

#This finished in about 2-3 days

```

### 3DNA pipeline
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/01_JuiceBox

git clone https://github.com/theaidenlab/3d-dna.git
cd 3d-dna/

module use /work/GIF/software/modules
module load GIF2/lastz/1.03.73
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
module load python
module load parallel/20170322-36gxsog
bash run-asm-pipeline.sh -e -m diploid /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/references/MisAssFixed.Pilon.fasta /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/aligned/merged_nodups.txt

runs in abot 7:06 hrs 16cpu.  144gb bam merged_nodups.txt file.

################################################################################
Number of scaffolds        968
Total size of scaffolds  160071602
   Longest scaffold   13754137
  Shortest scaffold          1
Number of scaffolds > 1K nt        925  95.6%
Number of scaffolds > 10K nt        766  79.1%
Number of scaffolds > 100K nt        158  16.3%
Number of scaffolds > 1M nt         17   1.8%
Number of scaffolds > 10M nt          2   0.2%
 Mean scaffold size     165363
Median scaffold size      25000
N50 scaffold length    5161912
 L50 scaffold count         10
 ##############################################################################

#Tried scaffolding with juicebox, and numerous scaffolding errors seemed to be present.  A problem with the masked genome removed coverage from some key areas.  Have to redo juicer using all parts of the genome.
```


#  Juicer with HI-C reads on an unmasked genome
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/02_juicerUnmasked
cp ../juicer/* .
ln -s ../juicer/debug
ln -s ../juicer/scripts
ln -s ../juicer/fastq
ln -s ../juicer/misc
ln -s ../../01_scnD2run/juicer/references/
ln -s ../../01_scnD2run/juicer/restriction_sites/

#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 96:00:00
#SBATCH -J juiceIt_0
#SBATCH -o juiceIt_0.o%j
#SBATCH -e juiceIt_0.e%j
#SBATCH --mail-user=remkv6@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
module use /work/GIF/software/modules
module load bwa
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
cd /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/02_juicerUnmasked
bash scripts/juicer.sh -y restriction_sites/MisAssFixed.Pilon.fasta_MboI.txt -z references/MisAssFixed.Pilon.fasta -p chrom.sizes
scontrol show job $SLURM_JOB_ID

#This ran for 96hrs on 16cpu and timed out during dedup.

Splitting file and running dedup in parallel.  Deleted dups/ optdups, merged_nodups files.  kept header.
```

### Run deduplication in parallel
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/02_juicerUnmasked/aligned

split -d  --additional-suffix=merged_sort.txt -n l/64  merged_sort.txt &
for f in x*txt; do mkdir $f.dir; cd $f.dir; cp ../../* .; ln -s ../../fastq/; ln -s ../../misc/ ;ln -s ../../debug/; ln -s ../../references/;ln -s ../../restriction_sites/;ln -s ../../scripts/;ln -s ../../splits/; mkdir aligned; mv ../$f aligned/merged_sort.txt; cd ../   ;done

for f in x*dir; do echo "cd "$f"; bash scripts/juicer.sh -S dedup -y restriction_sites/MisAssFixed.Pilon.fasta_MboI.txt -z references/MisAssFixed.Pilon.fasta -p chrom.sizes";done >>dedup.sh






All but three ended in less than 2 days:x40 x02 x13 will require further modification.  obtaining unique lines and removing those that were not finished processing (causes unix to think it is a binary)
cd ../../x02merged_sort.txt.dir/aligned/
cat dups.txt merged_nodups.txt merged_sort.txt optdups.txt |sort |uniq -c |awk '$1==1'  >RemainingMergedSort.txt &
tail -n +9 RemainingMergedSort.txt |sed 's/     //g' |sed 's/ /\t/g' |cut -f 3- |tr "\t" " " >RemainingMergedSort2.txt

cd ../../x13merged_sort.txt.dir/aligned/
cat dups.txt merged_nodups.txt merged_sort.txt optdups.txt |sort |uniq -c |awk '$1==1'  >RemainingMergedSort.txt &
tail -n +4 RemainingMergedSort.txt |sed 's/     //g' |sed 's/ /\t/g' |cut -f 3- |tr "\t" " " >RemainingMergedSort2.txt

cd ../../x40merged_sort.txt.dir/aligned/cat
cat dups.txt merged_nodups.txt merged_sort.txt optdups.txt |sort |uniq -c |awk '$1==1'  >RemainingMergedSort.txt &
tail -n +4 RemainingMergedSort.txt |sed 's/     //g' |sed 's/ /\t/g' |cut -f 3- |tr "\t" " " >RemainingMergedSort2.txt


split -d  --additional-suffix=merged_sort.txt -n l/16  RemainingMergedSort2.txt &
for f in x*txt; do mkdir $f.dir; cd $f.dir; cp ../../../../* .; ln -s ../../../../fastq/; ln -s ../../../../misc/ ;ln -s ../../../../debug/; ln -s ../../../../references/;ln -s ../../../../restriction_sites/;ln -s ../../../../scripts/;ln -s ../../../../splits/; mkdir aligned; mv ../$f aligned/merged_sort.txt; cd ../   ;done

for f in x*dir; do echo "cd "$f"; bash scripts/juicer.sh -S dedup -y restriction_sites/MisAssFixed.Pilon.fasta_MboI.txt -z references/MisAssFixed.Pilon.fasta -p chrom.sizes";done >>dedup.sh
cp /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/02_juicerUnmasked/aligned/x02merged_sort.txt.dir/aligned/makeSLURMp.py .
python makeSLURMp.py 16 dedup.sh
sbatch dedup_0.sub

#Concatenation and removal of binary lines from merged_nodups.
cd x02merged_sort.txt.dir/aligned/
mv merged_nodups.txt merged_nodups.txtOLD

cat merged_nodups.txtOLD x*.dir/aligned/merged_nodups.txt >merged_nodups.txt
cd ../../x13merged_sort.txt.dir/aligned/
mv merged_nodups.txt merged_nodups.txtOLD
cat merged_nodups.txtOLD x*.dir/aligned/merged_nodups.txt >merged_nodups.txt
cd ../../x40merged_sort.txt.dir/aligned/
mv merged_nodups.txt merged_nodups.txtOLD
cat merged_nodups.txtOLD x*.dir/aligned/merged_nodups.txt >merged_nodups.txt
cd ../../
cat  x*.dir/aligned/merged_nodups.txt >merged_nodups.txt



#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/02_juicerUnmasked/aligned/
mv merged_sort.txt Round1MergedSort.txt
mv test_nodups.txt merged_sort.txt
sbatch  juiceIt_0.sub

For some reason it could not make the hic file, likely due a broken juicer_tools.jar file in the scripts/ folder. # this is due to some split not getting incorporated in order.  Need to sort before running.  
sort -k2,2d -k6,6d --parallel=16 merged_nodups.txt >Sorted_merged_nodups.txt

module load bwa
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
java -Xmx2g -jar ../scripts/juicer_tools.jar pre merged_nodups.txt merged_nodups2.hic ../chrom.sizes

takes 7 hours to sort and create hic file
```

### run 3D DNA pipeline
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/02_juicebox/
git clone https://github.com/theaidenlab/3d-dna.git
cd 3d-dna/

module use /work/GIF/software/modules
module load GIF2/lastz/1.03.73
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
module load python
module load parallel/20170322-36gxsog
bash run-asm-pipeline.sh -m diploid /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/01_scnD2run/juicer/references/MisAssFixed.Pilon.fasta /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/02_juicerUnmasked/aligned/merged_nodups.txt
```
