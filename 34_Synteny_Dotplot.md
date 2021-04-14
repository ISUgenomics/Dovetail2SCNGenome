# Make a synteny dotplot for TN10 pseudo, tn10 draft, and X12 genomes



## X12 to TN10 pseudo
### Gather resources
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/59_SyntenyPlotMinimap
cp ../49_RenameChromosomes/01_Transfer2Box/RepeatMaskerFormatted.gff3.gz .
cp ../49_RenameChromosomes/01_Transfer2Box/SCNgenome.fasta .

cp ../01_CondoPorts/28_X12GenomeComparison/SCN_genome.fa X12Genome.fa
cp ../01_CondoPorts/28_X12GenomeComparison/repeat_annotation.gff X12repeats.gff
```


### Repeatmask before minimap alignment
```
ml bedtools2
bedtools maskfasta -fi SCNgenome.fasta -bed RepeatMaskerFormatted.gff3 -fo MaskedTN10.fasta

bedtools maskfasta -fi X12Genome.fa -bed X12repeats.gff -fo MaskedX12.fasta
```

### run minimap to generate PAF file
```

ml minimap2;sh runMinimap.sh MaskedX12.fasta MaskedTN10.fasta

#runMinimap.sh
##############################################################################
#!/bin/bash
query=$1
target=$2
outname="${query%.*}_${target%.*}_minimap2.paf"
module load minimap2
minimap2 -x asm5 -t 36 $target $query > ${outname}
##############################################################################

```

### Create dotplot
```
git clone https://github.com/tpoorten/dotPlotly.git
ml R


ml miniconda3
conda create -n dotplotly
source activate dotplotly
ml r/3.6.3-py3-sxv6dw3
R
install.packages(c("optparse", "ggplot2", "plotly")).
ml gcc/7.3.0-xegsmw4
git clone https://github.com/tpoorten/dotPlotly.git

 ./pafCoordsDotPlotly.R -i ../OnlyPseudoMaskedTN10_MaskedX12_minimap.paf -o SCNgenome_X12Gfuckenome_minimap2.dimbplot -l
```

## TN10 PSEUDO TO TN10 Draft
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/59_SyntenyPlotMinimap


ln -s ../01_CondoPorts/31_Synteny/02_738Assembly/genome738sl.polished.mitofixed.fasta
wget https://scnbase.org/files/download/genome738sl.polished.mitofixed.fa_.out_sorted.gff_.gz
gunzip genome738sl.polished.mitofixed.fa_.out_sorted.gff_.gz


bedtools maskfasta -fi genome738sl.polished.mitofixed.fasta -bed genome738sl.polished.mitofixed.fa_.out_sorted.gff_ -fo Masked738.fasta

echo "sh runMinimap.sh Masked738.fasta MaskedTN10.fasta" >minimap738.sh


```
