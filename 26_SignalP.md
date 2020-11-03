# Signal peptide prediction


### Signalp 5.0 Rerun with all transcripts assembled.  
```
wget https://services.healthtech.dtu.dk/download/66d87146-6d25-45b7-9bd9-62f21520bcbe/signalp-5.0b.Linux.tar.gz

tar -zxvf signalp-5.0b.Linux.tar.gz
cd signalp-5.0b/bin/

./signalp -fasta MergingTestSCNVHEJ_proteins.fasta -gff3


#so there are some duplicate mrna names I guess, so there may be some flaws to using this
less MergingTestSCNOnly_proteins_summary.signalp5|awk '$2!="OTHER"' |cut -f 1 |awk 'NR>2'|wc
   8898    8898  114428
[remkv6@condo045 bin]$ less MergingTestSCNOnly_proteins_summary.signalp5|awk '$2!="OTHER"' |cut -f 1 |awk 'NR>2'|sort|uniq|wc
   7964    7964  100750


less MergingTestSCNgffreadSecretedOnly.gff3|awk  '$2=="Secreted"' |awk '{print $4,"Secreted",$6,$7,$8,$9,$10,$11,"ID="$3 }' |tr " " "\t" >mRNAMergingTestSCNgffreadSecretedOnly.gff3

bedtools intersect -wo -f .7 -a <(awk '$3=="mRNA"' ../../../../25_AnnotateGenes/06_Combine/SCNgenomeFunctionalGeneAndMrnaAnnotations.gff3) -b mRNAMergingTestSCNgffreadSecretedOnly.gff3   | cut -f 1-9 |sort |uniq|sort -k1,1V -k4,5nr >AnnotatedSecretedmRNAMergingTestSCNgffread.gff3

#what are the number of secreted heterodera avenae gland expressed?
bedtools intersect -wo -f .7 -a <(awk '$3=="mRNA"' ../../../../25_AnnotateGenes/06_Combine/SCNgenomeFunctionalGeneAndMrnaAnnotations.gff3) -b mRNAMergingTestSCNgffread
SecretedOnly.gff3  |cut -f 1-9 |sort|uniq|grep -c -i "heterodera avenae "
69
## and the total number of annotations
awk '$3=="mRNA"' ../../../../25_AnnotateGenes/06_Combine/SCNgenomeFunctionalGeneAndMrnaAnnotations.gff3 |grep -c -i "heterodera avenae "
88

#what are the number of secreted gmapped effectors?
bedtools intersect -wo -f .7 -a <(awk '$3=="mRNA"' ../../../../25_AnnotateGenes/06_Combine/SCNgenomeFunctionalGeneAndMrnaAnnotations.gff3) -b ../../../../29_Effectors/SCNgenome.effector_sorted.gff |cut -f 10-18 |awk '$3=="mRNA"' |sort|uniq|wc
     91     819   25744

#and the total number of gmapped effectors?
awk '$3=="mRNA"' ../../../../29_Effectors/SCNgenome.effector_sorted.gff|sort|uniq|wc                                                                                       126    1134   35026


```


### SignalP 3.0 -- new gene names

```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/32_SignalP/01_SignalP3.1/signalp-3.0

#downloaded and extracted signalp 3.0

vi signalp
#modified AWK= , SIGNALP= , and HOW=
# Full path to the SignalP 3.0 software (MUST BE SET!):
SIGNALP=/work/gif/remkv6/Baum/04_DovetailSCNGenome/32_SignalP/01_SignalP3.1/signalp-3.0

# Other programs (assumed to be in the path)
# mandatory
SH=/bin/sh                              # POSIX-compliant shell
AWK=/usr/bin/gawk                               # nawk, gawk, or equivalent
PASTE=paste
# optional
PLOTTER=gnuplot                         # needed only for graphics
PPMTOGIF=ppmtogif                       # needed only for graphics
SIGNALPDISP="ghostview -landscape"      # needed only for graphics

#modified where the program finds HOW by adding ../
HOW=$SIGNALP/../how/how_$SYSTEM

#split the protein fasta into 8 parts
fasta-splitter.pl --n-parts 8 /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/mikado.loci.ancestralVHEJ_proteins.fasta

#create the run scripts
for f in *fasta; do echo "/work/gif/remkv6/Baum/04_DovetailSCNGenome/32_SignalP/01_SignalP3.1/signalp-3.0/signalp -t euk "$f" >"${f%.*}".out";done >signalp3.sh


paste <(cat *out |grep "Prediction:" -B 1 |grep ">") <(cat *out |grep "Max cleavage"  ) <(cat *out |grep "Signal peptide probability" ) |sort -k1,1V |paste - <(sort -k1,1V ProteinLengths.txt) |awk '$15>.499999' |awk '{print "samtools faidx mikado.loci.ancestralVHEJ_proteins.fasta "$1":"$11"-"$17" >>SignalPeptidesSubtracted.fasta"}' |sed 's/>Het/Het/g' >ExtractSignalPContainingProts.sh
sh ExtractSignalPContainingProts.sh

```


### Signalp 4.1 -- new gene names

```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/32_SignalP/03_Signalp4.1

#split the protein fasta into 8 parts
fasta-splitter.pl --n-parts 8 /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/mikado.loci.ancestralVHEJ_proteins.fasta


[remkv6@nova038 03_Signalp4.1]$ for f in *part* ; do echo "ml signalp;signalp -f short "$f" > "${f%.*}".signalP4.out"; done >signalp4.sh

```

### Signalp 5.0 -- final gene names
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/32_SignalP/02_SignalP5/signalp-5.0b/bin

fasta-splitter.pl --n-parts 8 /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/mikado.loci.ancestralVHEJ_proteins.fasta

for f in ../*fasta;do echo "./signalp -fasta "$f" -gff3"; done >signalp5.sh
python ~/common_scripts/makeSLURMs.py 1 signalp5.sh
sed -i 's/=36/=2/g' *sub
for f in *sub; do sbatch $f;done



ln -s /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/mikado.loci.ancestralVHEJ_proteins.fasta
bioawk -c fastx '{print $name,length($seq)}' /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/mikado.loci.ancestralVHEJ_proteins.fasta >ProteinLengths.txt


ml samtools
samtools faidx mikado.loci.ancestralVHEJ_proteins.fasta

cat *signalp5 |grep -v "#" |sort -k1,1nr |paste - <(sort -k1,1nr ProteinLengths.txt)  |awk '$3>.4999999' |awk '{print $1,$7,$12}' |sed 's/-/ /g' |sed 's/\./\t/2' |awk '{print "samtools faidx mikado.loci.ancestralVHEJ_proteins.fasta "$1":"$3"-"$4 " >>SignalPeptidesSubtracted.fasta"}' >ExtractSignalPContainingProts.sh
sh ExtractSignalPContainingProts.sh

cat *out |grep -v "#" |sort -k1,1V |paste - <(sort -k1,1V ProteinLengths.txt ) |awk '$9>.499999 {print "samtools faidx mikado.loci.ancestralVHEJ_proteins.fasta "$1":"$5"-"$14" >>SignalPeptidesSubtracted.fasta"}' >ExtractSignalPContainingProts.sh

sh ExtractSignalPContainingProts.sh



```
