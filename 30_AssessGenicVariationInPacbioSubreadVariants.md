
# Try to figure out what genes are variable in the subreads.
```
#extract subreads that are variants, only considering those with at least 1kb of variation
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/03_Subread/candidates
cat *bed |awk '($3-$2)>1000' |sed 's/|/\t/1' |sed 's/|/\t/1' |sed 's/|/\t/1'|cut -f 10|sed 's/]/\t/g' |sed 's/|/\t/1' |cut -f 2 |sed '/^$/d' |sort|uniq>SubreadsToExtract.list

#/work/GIF/remkv6/Baum/04_Dovetail2Restart/45_PacbioVariants

#extract variant subreads
ln -s ../30_ReadMapping/AllSubreads.fastq
sed -n '1~4s/^@/>/p;2~4p' AllSubreads.fastq >AllSubreads.fasta
less SubreadsToExtract.list |cdbyank AllSubreads.fasta.cidx >VariantSubreads.fasta

#map genes to subreads
cp ../28_X12GenomeComparison/01_AlignX12Genes2OurGenome/runGmap.sh .
ln -s ../25_AnnotateGenes/07_NewGenes/OrderedSCNGenePredictions_VHEJtranscripts.fasta
ln -s ../25_AnnotateGenes/07_NewGenes/OrderedSCNGenePredictions_VHEJtranscripts.fasta

echo "sh runGmap.sh VariantSubreads /work/GIF/remkv6/Baum/04_Dovetail2Restart/45_PacbioVariants/ VariantSubreads.fasta OrderedSCNGenePredictions_VHEJtranscripts.fasta" >mapGenes.sh


How many subreads with a 1kb variant had a gene mapping to them?
less VariantSubreads.OrderedSCNGenePredictions_VHEJtranscripts.gff3 |awk '$3=="gene"'|cut -f 1 |wc
  22953   22953 1794320

How many full genes could be mapped to these pacbio sureads?
less VariantSubreads.OrderedSCNGenePredictions_VHEJtranscripts.gff3 |awk '$3=="gene"'|cut -f 1 |sort|uniq|wc
   2831    2831  221749

#how many genes had more than 5 subreads with variants affecting genes?
less VariantSubreads.OrderedSCNGenePredictions_VHEJtranscripts.gff3 |awk '$3=="gene"'|cut -f 9 |sed 's/;/\t/g'|cut -f 2 |sort|uniq -c |sort -k1,1nr  |awk '$1>5 {print $2}' |sed 's/Name=//g' |while read line; do grep -w $line ../25_AnnotateGenes/06_Combine/CombeinAnnotIPRs.tab2;done  |wc
     45     435   10446

#What are the functions of the genes that are affected by variation in subreads?
less VariantSubreads.OrderedSCNGenePredictions_VHEJtranscripts.gff3 |awk '$3=="gene"'|cut -f 9 |sed 's/;/\t/g'|cut -f 2 |sort|uniq -c |sort -k1,1nr  |awk '$1>5 {print $2}' |sed 's/Name=//g' |while read line; do grep -w $line ../25_AnnotateGenes/06_Combine/CombeinAnnotIPRs.tab2;done  |cut -f 1 |sort|less
###############################################################################
Top 8 genes enriched in the variant pacbio reads are:
MAGL2_HUMAN MAGE-like protein 2 – a circadian rhythm gene
Mitochondrial RNA pseudouridine synthase RPUSD4
WD40 repeat containing protein
Neurogenic locus notch homolog protein 1
TRPL translocation defect protein 14
histone-lysine N-methyltransferase SETMAR-like
THIF-type NAD FAD binding fold and Ubiquitin-activating enzyme and Ubiquitin-activating enzyme repeat and Ubiquitin-activating enzyme e1 domain containing protein
AP-3 complex subunit delta-1
###############################################################################

```
