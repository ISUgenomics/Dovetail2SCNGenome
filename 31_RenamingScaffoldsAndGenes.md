# Name scaffolds by size, and rename genes/mRNAs something standard and relevant



### Rename Scaffolds
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/49_RenameChromosomes

bioawk -c fastx '{print $name, length($seq)}' SCNgenome.fasta |sort -k2,2nr  |awk '{print $1"\tChromosome_"NR }' >chromosome.map

singularity shell "/opt/rit/singularity/images/maker/2.31.10_3.1/maker.simg"

map_data_ids chromosome.map mRNALocationssorted
map_data_ids chromosome.map RiboCoords.bed
map_data_ids chromosome.map Repeatmasker.gff
map_data_ids chromosome.map SCNgenome.fasta.EDTA.TEanno.gff
map_data_ids chromosome.map SCNgenome.effector_sorted.gff
map_data_ids chromosome.map TRF.bed
map_data_ids chromosome.map Centromere.gff
map_data_ids chromosome.map mikado.loci.ancestral.gff3

map_fasta_ids chromosome.map SCNgenome.fasta


```
### Rename Genes
```
#change locus to gene, or it is not recognized.  
awk '{if($3=="locus") {print $1,$2,"gene",$4,$5,$6,$7,$8,$9} else {print $1,$2,$3,$4,$5,$6,$7,$8,$9}}'  mikado.loci.ancestral.gff3 |tr " " "\t" >mikado.lociGenes.gff3

singularity shell "/opt/rit/singularity/images/maker/2.31.10_3.1/maker.simg"

#Create map of old gene name to new gene name, mRNAs too.
maker_map_ids --prefix Hetgly --suffix .t --iterate 1 --justify 5 mikado.lociGenes.gff3 >GeneNames.map


#have to fix a parent child relationship with the genes and mRNAs
awk '{if($3=="locus") {print $1,$2,"gene",$4,$5,$6,$7,$8,$9} else {print $1,$2,$3,$4,$5,$6,$7,$8,$9}}'  mikado.loci.ancestral.gff3 |tr "
 " "\t" |sed 's/locus=/Parent=/g' |sed 's/;/\t/g' |awk '{if($3=="gene"){print $1,$2,$3,$4,$5,$6,$7,$8,$9";"}else if($3=="mRNA") {print $1,$2,$3,$4,$5,$6,$7,$8,$9";"$1
2";"}else {print $1,$2,$3,$4,$5,$6,$7,$8,$9";"$10}} ' |tr " " "\t" |awk 'NR>1' >PolishedForNameChange.gff3

#rename the gff genes/mRNAs
map_gff_ids  GeneNames.map PolishedForNameChange.gff3

#change name of transcripts and proteins generated
map_fasta_ids GeneNames.map mikado.loci.ancestralVHEJ_proteins.fasta
map_fasta_ids GeneNames.map mikado.loci.ancestralVHEJ_transcripts.fasta

#Get gene names changed to reflect gff changes
map_data_ids GeneNames.map mRNALocations.tab
```
