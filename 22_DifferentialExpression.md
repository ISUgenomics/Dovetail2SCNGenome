# Generate differential expression for gland and whole worm RNASEQ


# generate bigwig and count files
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate

ln -s ../../25_AnnotateGenes/SCNgenomeFunctionalGeneAn
ln -s ../../10_RepeatModeler/SCNgenome.fasta
for f in  ../../11_AlignRNA/*sorted.bam; do ln -s $f;done



for f in *bam; do echo "ml samtools; samtools index "$f"; ml subread; featureCounts -T 16 -p -t gene -s 2 -g ID -a SCNgenomeFunctionalGeneAnnotationsHighConfidenceAllFeatures.gff3  -o "${f%.*}"_GeneCounts.txt " $f"; ml bedtools2; bedtools genomecov -ibam "$f" -bga -g SCNgenome.fasta > "${f%.*}"_sorted.bam.bdg; ml bioawk; bioawk -c fastx '{print \$name,length(\$seq)}'  SCNgenome.fasta >Chr.sizes ; /home/remkv6/common_scripts/kentUtils/bin/linux.x86_64/bedGraphToBigWig "${f%.*}"_sorted.bam.bdg Chr.sizes "${f%.*}"_sorted.bw"; done
```
