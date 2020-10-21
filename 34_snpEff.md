# Run SNPEff

#### Rationale
```
Need to find proteins with functional coding differences in the genome in other SCN lines
Particularly, we are interested in seeing how coding differences can change a signal peptide.
```
### File and database setup
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/43_FreeBayes/

mkdir 01_SNPEFF
cd 01_SNPEFF/

#softlink snp file
ln -s ../AllSNPS.vcf

#install snpeff
ml miniconda3
conda create -n snpEff
source activate snpEff
conda install -c bioconda snpeff

#copy in the genome file
cp ../SCNgenome.fasta .

#make the data structure
mkdir data/
cd data/
mkdir Hetgly
cd Hetgly


#create chromosome length file in gff format
ml bioawk
bioawk -c fastx '{print $name,"bioawk","chromosome","1",length($seq),".",".",".","ID="$name}' ..//../SCNgenome.fasta |tr " " "\t" >chromosome.gff

#create the gff in the format that snpeff expects
cat <(grep -v "#" ../../../../49_RenameChromosomes/01_Transfer2Box/OrderedSCNGenePredictions.gff3 ) chromosome.gff ../../SCNgenome.fasta >genes.gff

#add these lines between the gff's end and the beginning of the fasta
 vi genes.gff
 ######################################################################################3
 etc.
 contig_10383_pilon_pilon_pilon  AUGUSTUS        gene    897     8132    0.31    +       .       .
 contig_10383_pilon_pilon_pilon  AUGUSTUS        chromosome 897     1159    .    .       0       ID=contig_10383_pilon_pilon_pilon
 ###
 ##FASTA
 >contig_1038_pilon_pilon_pilon
 TGATATCATCAAAGCGCTATGTCCCGTGGAAGCTAGTTTGGGATGATAAATATAGAATACTAAACGAGCTCCGTATAATA
 etc.
 #########################################################################################

 #get in proper directory
 cd ../../

 vi snpEff.config
 ########################
 Hetgly.genome : SCN
 sacCer.genome : Yeast
 ########################

 #this takes time.
 snpEff build -gff3 Hetgly
```

### Run SNPEff
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/43_FreeBayes/01_SNPEFF

echo "ml miniconda3; source activate snpEff;  java -Xmx950G  -jar ~/.conda/pkgs/snpeff-5.0-0/share/snpeff-5.0-0/snpEff  -interval data/Hetgly/genes.gff Hetgly AllSNPS.vcf > snpEffOUT.vcf" >runSNPeff.sh

#note, it looks like this program requires a large node, but may not be necessary for SCN.  Also this needs to be run on each line separately, so need to check that out before running.
```
