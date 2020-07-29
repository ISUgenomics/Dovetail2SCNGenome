# Figure out which of these proteins is a peptidase or peptidase inhibitor
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/41_Merops

#get the blast database
for f in ../../11_PeptidaseInvestigation/04_Merops/protease.lib.p*; do ln -s $f;done
ln -s ../25_AnnotateGenes/07_NewGenes/OrderedSCNGenePredictionsVHEJ_proteins.fasta


ml blast-plus; blastp -db protease.lib -query mikado.loci.ancestralVHEJ_proteins.fasta -num_threads 16 -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -out SCNProts2Merops.blastout

#filter by evalue less than 0.001, print mrna, evalue, merops annotation
 less SCNProts2Merops.blastout |awk '$14<.001' |sort -k1,1 -u |cut -f 1,14,17 >Merops.tab


 wc Merops.tab
 9303  132007 1549663 Merops.tab


```
