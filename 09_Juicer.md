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

#Current directory
/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/juicer
sh /scripts/juicer.sh



Dependencies:

module load bwa
module load gnutls/3.5.13-7a3mvfy

#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
```
