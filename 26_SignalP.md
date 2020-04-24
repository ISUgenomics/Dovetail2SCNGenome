# Signal peptide prediction

### SignalP 5.0
```
wget https://services.healthtech.dtu.dk/download/66d87146-6d25-45b7-9bd9-62f21520bcbe/signalp-5.0b.Linux.tar.gz

tar -zxvf signalp-5.0b.Linux.tar.gz
cd signalp-5.0b/bin/

ln -s /work/GIF/remkv6/Baum/04_Dovetail2Restart/25_AnnotateGenes/07_NewGenes/OrderedSCNGenePredictionsVHEJ_proteins.fasta

./signalp -fasta OrderedSCNGenePredictionsVHEJ_proteins.fasta -gff3
less OrderedSCNGenePredictionsVHEJ_proteins_summary.signalp5 |grep -v "#" |awk '{print $1,$3}' |tr " " "\t" >SignalP5.tab

### how many??
less OrderedSCNGenePredictionsVHEJ_proteins_summary.signalp5 |grep -v "#" |awk '{print $1,$3}' |tr " " "\t" |awk '$2>.6' |wc
   3437    6874   64432
[remkv6@condo042 bin]$ less OrderedSCNGenePredictionsVHEJ_proteins_summary.signalp5 |grep -v "#" |awk '{print $1,$3}' |tr " " "\t" |awk '$2>.5' |wc
   3714    7428   69635
[remkv6@condo042 bin]$ less OrderedSCNGenePredictionsVHEJ_proteins_summary.signalp5 |grep -v "#" |awk '{print $1,$3}' |tr " " "\t" |awk '$2>.8' |wc
   2841    5682   53265
[remkv6@condo042 bin]$ less OrderedSCNGenePredictionsVHEJ_proteins_summary.signalp5 |grep -v "#" |awk '{print $1,$3}' |tr " " "\t" |awk '$2>.9' |wc
   2353    4706   44102

```

#### SignalP 4.1

```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/32_SignalP
ln -s /work/GIF/remkv6/Baum/04_Dovetail2Restart/25_AnnotateGenes/07_NewGenes/OrderedSCNGenePredictionsVHEJ_proteins.fasta

ml signalp/4.1f-py2-binsfdx

fasta-splitter.pl --n-parts 100 OrderedSCNGenePredictionsVHEJ_proteins.fasta
for f in *part* ; do signalp -f short $f > ${f}.signalP4.out; done

 cat *signalP4.out |grep -v "#" |awk '{print $1"\t"$9}' >OrderedSCNGenePredictionsVHEJ_proteinsSignalP4.tab

 awk '$2>.5' OrderedSCNGenePredictionsVHEJ_proteinsSignalP4.tab |wc
   3823    7646   60188

```
