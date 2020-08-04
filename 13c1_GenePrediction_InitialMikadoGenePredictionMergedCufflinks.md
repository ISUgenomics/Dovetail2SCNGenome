# Investigate initial mikado prediction reduced via cufflinks merge
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/47_MikadoFinalize

ln -s ../23_Mikado/mikado.loci.ancestral.gff3
ln -s ../23_Mikado/mikado.loci.ancestralVHEJ_proteins.fasta
ln -s ../23_Mikado/mikado.loci.ancestralVHEJ_transcripts.fasta
```

#How many mRNA's in the gene prediction are clearly transposable elements? i.e. overlap with EDTA track
```
ml bedtools2/2.27.1-s2mtpsu

#60% overlap
bedtools intersect -f .6 -wo -a <(awk '$3=="mRNA"' mikado.loci.ancestral.gff3) -b ../33_EDTA/EDTA/SCNgenome.fasta.EDTA.TEanno.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
   1653    1653   40865

#50% overlap
bedtools intersect -f .5 -wo -a <(awk '$3=="mRNA"' mikado.loci.ancestral.gff3) -b ../33_EDTA/EDTA/SCNgenome.fasta.EDTA.TEanno.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
2026    2026   50085

#40% overlap
bedtools intersect -f .4 -wo -a <(awk '$3=="mRNA"' mikado.loci.ancestral.gff3) -b ../33_EDTA/EDTA/SCNgenome.fasta.EDTA.TEanno.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
2523

#30% overlap
bedtools intersect -f .3 -wo -a <(awk '$3=="mRNA"' mikado.loci.ancestral.gff3) -b ../33_EDTA/EDTA/SCNgenome.fasta.EDTA.TEanno.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
3206

#20% overlap
bedtools intersect -f .2 -wo -a <(awk '$3=="mRNA"' mikado.loci.ancestral.gff3) -b ../33_EDTA/EDTA/SCNgenome.fasta.EDTA.TEanno.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
4431
#10% overlap
bedtools intersect -f .1 -wo -a <(awk '$3=="mRNA"' mikado.loci.ancestral.gff3) -b ../33_EDTA/EDTA/SCNgenome.fasta.EDTA.TEanno.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
7544
```
