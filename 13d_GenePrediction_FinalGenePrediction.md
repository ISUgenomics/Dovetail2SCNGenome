# Need to predict genes for the SCN TN10 pseudomolecule assembly

### Run Mikado on Filtered transcriptome's
```

#list4.txt
#################################################################################
MetAllRNASEQClass2.proteins.gff3        CL      False
Metaugustus.hints.proteins.gff3 BR      False
MetNonRrnaRNASEQ_stringtie.protein.gff3 ST      False
MetSCNgenome.DovetailSCNMaker4.all.maker.transcripts.proteins.gff3      MA      False   -1
MetSCNgenome.Trinity-GG.proteins.gff3   TR      False
MetSCNgenome.TylenchidaESTNotH.glycines.proteins.gff3   TY      False   -1
MetSpadesSCNgenome.transcripts.proteins.gff3    SP      False
FixedBrakerUnmasked.gff3        BS      True

#################################################################################

#MikadoScript_0.sub
################################################################################
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=35
#SBATCH -t 96:00:00
#SBATCH -J MikadoScript_0
#SBATCH -o MikadoScript_0.o%j
#SBATCH -e MikadoScript_0.e%j
#SBATCH --mail-user=remkv6@istate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
ulimit -s unlimited
cd /work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn
ml miniconda3
source activate mikado
ml py-networkx/1.10-py2-2rrvpu6
#!/bin/bash
#setup variables
genome=SCNgenome.fasta
bam="AllRNASEQ_sorted.bam"
list="list4.txt"
#run splice junction prediction
junctions="portcullis_filtered.pass.junctions.bed"
echo "launch configure"
mikado configure \
   --list $list \
   --reference $genome \
   --mode stringent \
   --scoring worm.yaml \
   --copy-scoring worm.yaml \
   --junctions $junctions configuration.yaml
echo "launch prepare"
mikado prepare \
   --json-conf configuration.yaml
   blast db
   makeblastdb \
      -in uniprotCombined.fasta \
      -dbtype prot \
      -parse_seqids
   echo "start blast"
   blastx \
      -max_target_seqs 5 \
      -num_threads 35 \
      -query mikado_prepared.fasta \
      -outfmt 5 \
      -db uniprotCombined.fasta \
      -evalue 0.000001 2> blast.log | sed '/^$/d' > mikado.blast.xml
   blastxml=mikado.blast.xml
   echo "launch transdecoder"
   TransDecoder.LongOrfs \
      -t mikado_prepared.fasta
   TransDecoder.Predict \
      -t mikado_prepared.fasta \
      --cpu 35
   orfs=$(find $(pwd) -name "mikado_prepared.fasta.transdecoder.bed")
   echo "Please do not produce tons of crappy incorrect annotations, I mean launch serialise"
   mikado serialise \
      --start-method spawn \
      --procs 35 \
      --blast_targets ${genome} \
      --json-conf configuration.yaml \
      --xml ${blastxml} \
      --orfs ${orfs} \
      --blast_targets uniprotCombined.fasta \
      -mr 1

   echo "launch pick"
   mikado pick \
      --start-method spawn \
      --procs 35 \
      --json-conf configuration.yaml \
      --subloci_out mikado.subloci.gff3 \
      --pad
      -v
   scontrol show job $SLURM_JOB_ID



Cannot seem to produce mikado output with reliable splicing, start codons, and stop codons.  So filter all proteins that do not meet that criteria.
```
### Filter transcripts going into mikado
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn/03_MethionineProteinsOnly

#braker genes, masked
gffread ../augustus.hints.gff3 -g SCNgenome.fasta -t mRNA -VHEJ -x BrakerMasked_transcriptsVHEJ.fasta -y BrakerMasked_proteinsVHEJ.fasta
grep -c ">" BrakerMasked_proteinsVHEJ.fasta
24113

#maker genes on 368 scaffold genome
gffread ../SCNgenome.DovetailSCNMaker4.all.maker.transcripts.gff3  -g SCNgenome.fasta -t mRNA -VHEJ  -x SCNgenome.DovetailSCNMaker4.all.maker.transcripts_transcriptsVHEJ.fasta -y SCNgenome.DovetailSCNMaker4.all.maker.transcripts_proteinsVHEJ.fasta
grep -c ">" SCNgenome.DovetailSCNMaker4.all.maker.transcripts_proteinsVHEJ.fasta
39896

#spades transcripts
gffread ../SCNgenome.transcripts.gff3  -g SCNgenome.fasta -t mRNA -VHEJ  -x Spades_transcriptsVHEJ.fasta -y Spades_proteinsVHEJ.fasta
grep -c ">" Spades_proteinsVHEJ.fasta
23598

#Trinity transcripts
gffread  ../SCNgenome.Trinity-GG.gff3  -g SCNgenome.fasta -t mRNA -VHEJ  -x Trinity_transcriptsVHEJ.fasta -y Trinity_proteinsVHEJ.fasta
grep -c ">" Trinity_proteinsVHEJ.fasta
41366

# NCBI EST proteins
gffread  ../SCNgenome.TylenchidaESTNotH.glycines.gff3  -g SCNgenome.fasta -t mRNA -VHEJ  -x NCBIEST_transcriptsVHEJ.fasta -y NCBIEST_proteinsVHEJ.fasta
[remkv6@nova005 03_MethionineProteinsOnly]$ grep -c ">" NCBIEST_proteinsVHEJ.fasta
1750

# Class2 transcripts
gt gtf_to_gff3 -tidy <(sort -k1,1V -k4,5nr AllRNASEQClass2.gtf|grep -v "#" |uniq) >AllRNASEQClass2.gff3
gt gff3 -sortlines -tidy -fixregionboundaries  AllRNASEQClass2.gff3 >tidiedAllRNASEQClass2.gff3
gt cds -startcodon -finalstopcodon -seqfile SCNgenome.fasta  -matchdesc tidiedAllRNASEQClass2.gff3 >FixedClass2.gff3
gffread ../FixedClass2.gff3 -g SCNgenome.fasta -t mRNA -VHEJ  -x AllRNASEQClass2_transcriptsVHEJ.fasta -y AllRNASEQClass2_proteinsVHEJ.fasta
grep -c ">" AllRNASEQClass2_proteinsVHEJ.fasta
36212

#stringtie conversion

gt gtf_to_gff3 -tidy <(sort -k1,1V -k4,5nr NonRrnaRNASEQ_stringtie.gtf|grep -v "#" |uniq) >NonRrnaRNASEQ_stringtie.gff3
gt gff3 -sortlines -tidy -fixregionboundaries  NonRrnaRNASEQ_stringtie.gff3 >tidiedNonRrnaRNASEQ_stringtie.gff3
gt cds -startcodon -finalstopcodon -seqfile SCNgenome.fasta  -matchdesc tidiedNonRrnaRNASEQ_stringtie.gff3 >FixedNonRrnaRNASEQ_stringtie.gff3
gffread ../FixedNonRrnaRNASEQ_stringtie.gff3 -g SCNgenome.fasta -VHEJ -t mRNA -x Stringtie_transcriptsVHEJ.fasta -y Stringtie_proteinsVHEJ.fasta
grep -c ">" Stringtie_proteinsVHEJ.fasta
29063

#braker unmasked

gt gtf_to_gff3 -tidy <(sort -k1,1V -k4,5nr BrakerUnmasked.gff|grep -v "#" |uniq) >BrakerUnmasked.gff3
gt gff3 -sortlines -tidy -fixregionboundaries  BrakerUnmasked.gff3 >tidiedBrakerUnmasked.gff3
gt cds -startcodon -finalstopcodon -seqfile SCNgenome.fasta  -matchdesc tidiedBrakerUnmasked.gff3 >FixedBrqakerUnmasked.gff3
gffread ../FixedBrakerUnmasked.gff3 -g SCNgenome.fasta -t mRNA -VHEJ -x BrakerUnmasked_transcriptsVHEJ.fasta -y BrakerUnmasked_proteinsVHEJ.fasta
grep -c ">" BrakerUnmasked_proteinsVHEJ.fasta
37257
```



#### Rerun mikado with these filtered gene sets
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn

less mikado_proteins.fasta |tr "\n" "\t" |sed 's/>/\n>/g' |awk '{print substr($3,1,1)}' |sort|uniq -c |less
#Results
#########################
1
177 .
374 A
297 C
290 D
239 E
324 F
331 G
212 H
440 I
477 K
1877 L
34562 M
525 N
253 P
230 Q
360 R
1414 S
285 T
343 V
197 W
31 X
181 Y
#########################

#So how many are good?
less mikado_proteinsVHEJ.fasta |tr "\n" "\t" |sed 's/>/\n>/g' |awk '{print substr($3,1,1)}' |sort|uniq -c |less
1
29959 M

#Create a format of 'mrna_name\tgene_name'
grep ">" mikado_proteinsVHEJ.fasta|awk '{print $1}' |sed 's/>//g' |cat - AllProteins.list |sort|uniq -c |awk '$1==1{print $2,$2}' |sed 's/\./\t/4' |cut -f 1  |tr " " "\t" >BadGenes.list
#How many were bad?
wc -l ../BadGenes.list
13460 ../BadGenes.list


#remove them from the gff
mikado util grep -v BadGenes.list mikado.loci.gff3 RemoveBadGenesmikado.loci.gff3

#how many genes are there?
less RemoveBadGenesmikado.loci.gff3 |awk '$3=="gene"' |wc
  26016  234144 4691582
#how many mrnas/proteins are there?
less RemoveBadGenesmikado.loci.gff3 |awk '$3=="mRNA"' |wc
  29959  269631 5155985

#get the bad ones for assessment
mikado util grep BadGenes.list mikado.loci.gff3 BadGenesmikado.loci.gff3

#get bad genes and compare to braker unmasked prediction to obtaion appropriate CDS for genes.
#how many will I get in the exchange?

bedtools intersect -wo -a ../BadGenesmikado.loci.gff3 -b ../FixedBrakerUnmasked.gff3 |awk '$3=="gene"' |awk '$12=="gene"' |cut -f 18 |sort|uniq|sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
   8893    8893   86151



#Make a tabular file of mrna name, gene name
awk '$3=="mRNA"' ../FixedBrakerUnmasked.gff3 |cut -f 9 |sed 's/ID=//g' |sed 's/Parent=//g' |sed 's/;/\t/g' |cut -f 1,2 >brakerunmaskedGrepMod.list

bedtools intersect -wo -a ../BadGenesmikado.loci.gff3 -b ../FixedBrakerUnmasked.gff3 |awk '$3=="gene"' |awk '$12=="gene"' |cut -f 18 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|grep -w -f - brakerunmaskedGrepMod.list >List2GrabFromBraker.list

mikado util grep List2GrabFromBraker.list ../FixedBrakerUnmasked.gff3 AugustusGenes4Concat.gff3
less AugustusGenes4Concat.gff3 |grep -v "^#" |awk '{if($3=="gene" || $3=="mRNA") {print $1,"AUGUSTUS",$3,$4,$5,$6,$7,$8,$9,$9} else {print $1,"AUGUSTUS",$3,$4,$5,$6,$7,$8,$9}}' |tr " " "\t" |sed 's/;/\t/2' |cut -f 1-10 |sed 's/\tID=/;Name=/2'  >augustusreformatName.gff3

#256 genes are bad.  lets remove
less augustusreformatName.gff3 |awk '$5-$4<3' |awk '$3=="mRNA"' |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |sed 's/Parent=//g' |cut -f 1,2 >BadAugustusGenes.list
mikado util grep -v BadAugustusGenes.list augustusreformatName.gff3 FixedaugustusreformatName.gff3

gffread FixedaugustusreformatName.gff3 -VHEJ -g ../SCNgenome.fasta -t mRNA -x FixedaugustusreformatName_VHEJtranscripts.fasta -y FixedaugustusreformatName_VHEJproteins.fasta
gffread FixedaugustusreformatName.gff3  -g ../SCNgenome.fasta -t mRNA -x FixedaugustusreformatName_transcripts.fasta
 -y FixedaugustusreformatName_proteins.fasta
 cat <(grep ">" FixedaugustusreformatName_proteins.fasta) <(grep ">" FixedaugustusreformatName_VHEJproteins.fasta) |sort|uniq -c |awk '$1==1 {print $2,$3}' |sed 's/gene=//g' |sed 's/>//g' |tr " " "\t" >>BadAugustusGenes.list

 mikado util grep -v BadAugustusGenes.list augustusreformatName.gff3 FixedaugustusreformatName.gff3

cat FixedaugustusreformatName.gff3 ../RemoveBadGenesmikado.loci.gff3|grep -v "^#" >NeedsRenamedFinalGenes.gff3

gt gff3 -tidy -sortlines -checkids  NeedsRenamedFinalGenes.gff3  >SCNGenePredictions.gff3


#get the proteins  and transcript sequences
gffread SCNGenePredictions.gff3 -VHEJ -g ../SCNgenome.fasta -t mRNA -x SCNGenePredictions_VHEJtranscripts.fasta -y SCNGenePredictionsVHEJ_proteins.fasta

#Verify that there are the right number of proteins and mrnas in the gff
 grep -c ">" SCNGenePredictionsVHEJ_proteins.fasta
39181
awk '$3=="mRNA"' SCNGenePredictions.gff3 |wc
39181

perl gff3sort/gff3sort.pl --precise --chr_order natural SCNGenePredictions.gff3 >OrderedSCNGenePredictions.gff3
bgzip OrderedSCNGenePredictions.gff3
tabix -p gff OrderedSCNGenePredictions.gff3.gz


#add the ncrna track to jbrowse
gt gff3 -tidy -sortlines  Ncrnas.gff3  >TidyNCrna.gff3
perl gff3sort/gff3sort.pl --precise --chr_order natural TidyNCrna.gff3 >OrderedTidyNCrna.gff3
bgzip OrderedTidyNCrna.gff3
tabix -p gff OrderedTidyNCrna.gff3.gz
```

### Investigate genes, How many are clearly transposable elements?
```

#How many mRNA's in the gene prediction overlap with repetitive elements?

#60% overlap
(mikado) [remkv6@nova023 05_FinalGenePrediction]$ bedtools intersect -f .6 -wo -a <(awk '$3=="mRNA"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
   5150    5150   49819
#50% overlap
(mikado) [remkv6@nova023 05_FinalGenePrediction]$ bedtools intersect -f .5 -wo -a <(awk '$3=="mRNA"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
   6133    6133   59294
#30% overlap
(mikado) [remkv6@nova023 05_FinalGenePrediction]$ bedtools intersect -f .3 -wo -a <(awk '$3=="mRNA"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
   8700    8700   84152
#20% overlap   
(mikado) [remkv6@nova023 05_FinalGenePrediction]$ bedtools intersect -f .2 -wo -a <(awk '$3=="mRNA"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
  10813   10813  104663

#Any intersect
(mikado) [remkv6@nova023 05_FinalGenePrediction]$ bedtools intersect  -wo -a <(awk '$3=="mRNA"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
  25929   25929  251739


#What if we consider only CDS overlaps with repeats
#Any overlap between any cds and any repeat
bedtools intersect -f .5 -wo -a <(awk '$3=="CDS"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/Parent=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
14159   14159  137011

# 50% of the CDS length, only one CDS required
bedtools intersect -f .5 -wo -a <(awk '$3=="CDS"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/Parent=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc
  11791   11791  114040

# requires at least 2 CDS to overlap any part of repeat
bedtools intersect -f .5 -wo -a <(awk '$3=="CDS"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/Parent=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq -c |awk '$1>1' |wc
  9775   19550  172686
#requires at least 30% overlap of cds by the repeat, and two CDS's
bedtools intersect -f .3 -wo -a <(awk '$3=="CDS"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/Parent=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq -c |awk '$1>1' |wc
   8313   16626  146776

#Probably the most conservative measure, same as above
bedtools intersect -f .3 -wo -a <(awk '$3=="CDS"' SCNGenePredictions.gff3) -b ../04_Expression/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/Parent=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq -c |awk '$1>1 {print $2}' >RepetitiveBy30PercOverlap2CDS.list

#How many of the mRNAs have an annotation that is transposon, helitron, or transposase?
bedtools intersect -wo -f .5 -r -a <(awk '$3=="mRNA"'  SCNGenePredictions.gff3)  -b <(awk '$3=="mRNA"'  ../OldPredictionSCNgenomeFunctionalGeneAnnotations.gff3 ) |grep -i  "transposon" |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort |uniq >TransposonGenes.list
bedtools intersect -wo -f .5 -r -a <(awk '$3=="mRNA"'  SCNGenePredictions.gff3)  -b <(awk '$3=="mRNA"'  ../OldPredictionSCNgenomeFunctionalGeneAnnotations.gff3 ) |grep -i  "transposase" |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort |uniq >TransposaseGenes.list
bedtools intersect -wo -f .5 -r -a <(awk '$3=="mRNA"'  SCNGenePredictions.gff3)  -b <(awk '$3=="mRNA"'  ../OldPredictionSCNgenomeFunctionalGeneAnnotations.gff3 ) |grep -i  "helitron" |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort |uniq >HelitronGenes.list
cat TransposonGenes.list TransposaseGenes.list HelitronGenes.list |sort|uniq >GenesWithTransposonAnnotations.list
less GenesWithTransposonAnnotations.list|wc
    279     279    2707
```

# Are all of the effectors there?
```

#identify by name only "gland" and "effector"
bedtools intersect -wo -f .5 -r -a <(awk '$3=="mRNA"'  SCNGenePredictions.gff3)  -b <(awk '$3=="mRNA"'  ../OldPredictionSCNgenomeFunctionalGeneAnnotations.gff3 ) |grep -i  "effector" |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort |uniq >effector.list
bedtools intersect -wo -f .5 -r -a <(awk '$3=="mRNA"'  SCNGenePredictions.gff3)  -b <(awk '$3=="mRNA"'  ../OldPredictionSCNgenomeFunctionalGeneAnnotations.gff3 ) |grep -i  "gland" |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort |uniq >>glandGenes.list
cat glandGenes.list effector.list |sort|uniq >EffectorGenes.list
wc EffectorGenes.list
 396  396 3780 EffectorGenes.list


#how many of these effectors are repetitive as called by repeatmasker
cat <(sort EffectorGenes.list |uniq) <(sort RepetitiveBy30PercOverlap2CDS.list |uniq) |sort|uniq -c |awk '$1==2' |wc
    176     352    3031

```

### Examine the expression of genes, leaving repetitive elements
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/38_Expression/01_AllExpressionDatasetsSeparate

#zero expression
less DeseqTable.txt |awk '{print $1,$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16}' |awk '$2=="0"' |wc
  10379   20758  121089
#one read expression  
less DeseqTable.txt |awk '{print $1,$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16}' |awk '$2=="1"' |wc
   1981    3962   23137

# make a list of all mrna and gene names, grep to get them and get new gene count
awk '$3=="mRNA"' OrderedSCNGenePredictions.gff3 |cut -f 9 |sed 's/ID=//g' |sed 's/Parent=//g' |sed 's/;/\t/g' |cut -f 1,2 >mrnaGeneAll.list

#make gene list for mikado to remove from the gff
less DeseqTable.txt |awk '{print $1,$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16}' |awk '$2=="0" {print $1}' |grep -w -f - mrnaGeneAll.list >NoExpressedMrnaGene.list

#remove these from the annotation

```


# List of genes to keep for counts
```
awk '$3=="mRNA" ' SCNGenePredictions.gff3 |cut -f 9 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 >AllmRNAs.list

#Remove the overlap with repeats, and annotated transposon proteins
#Figure in non expressed transcripts later when featurecounts is available
cat AllmRNAs.list GenesWithTransposonAnnotations.list RepetitiveBy30PercOverlap2CDS.list |sort|uniq -c |awk '$1==1{print $2}' |cat - EffectorGenes.list |sort|uniq -c |awk '$1==1{print $2}' >NoRepeatFinalGene.List

wc NoRepeatFinalGene.List
 31233  31233 304022 NoRepeatFinalGene.List

```


### Merge good alignmetns with cufflinks gffread -M ### testing still 4/14/20
 While this does substantially lower the number of genes that are generated, it does combine genes with mrnas on two strands, which is unlikely.  Having them all displayed as a potential mRNAs may be helpful, but also display lots of spurious transcripts toois also
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn/03_MethionineProteinsOnly

#rename all gth alignments from each transcript assembler, so that gffread merge can tell the difference in the transcripts
ls */*aln |while read line; do awk '{print substr($1,1,5),$0}' $f;done |sed 's/AllRN/CLASS2/1' |sed 's/Brake/BRAKER/1' |sed 's/NCBIE/NCBIEST/1' |sed 's/SCNge/MAKER/1' |sed 's/Spade/SPADES/1' |sed 's/Strin/STRINGTIE/1' |sed 's/Trini/TRINITY/1' |while read line; do echo "sed -i 's/ID=/ID= "$line" /g' "$line" ; sed -i 's/Parent=/Parent=" $line "/g' ";done |awk '{print $1,$2,$3$4$6,$8,$9,$10,$11,$12$13$15,$14}' >RenameAllGenomethreaderAlignments.sh
sh RenameAllGenomethreaderAlignments.sh

#/work/gif/remkv6/Baum/04_DovetailSCNGenome/01_mikadoRerurn/05_FinalGenePrediction

#merge with all genomethreader alignments of Met start encoded proteins and to All H.glycines EST alignments
gffread -E -K -M -Q -VHEJ -g ../SCNgenome.fasta  OrderedSCNGenePredictions.gff3 ../03_MethionineProteinsOnly/*/*aln  -o AllProtsMetESTMergeOrderedSCNGenePredictions.gff3

#How many genes are there if we merge with all transcripts starting with M from the genomethreader alignments of all genome-guided transcriptomes?
bedtools intersect -wo -f .3 -a <(awk '$3=="locus"' AllProtsMetMergeESTmergeMergeOrderedSCNGenePredictions.gff3) -b FixedaugustusreformatName.gff3 |awk '$12=="gene"' |cut -f 1-9  |wc
 12518  112662 1999026

#make a gff to look at in jbrowse, genometools throws a fit because the parent ID combos are not unique via the parallel gth alignments.
awk '{if($3=="locus") {print $1,$2,"gene",$4,$5,$6,$7,$8,$9} else {print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' AllProtsMetESTMergeOrderedSCNGenePredictions.gff3 |tr " " "\t" >MergingTestSCN.gff3


ml tabix/2013-12-16-py2-talztjl
bgzip MergingTestSCN.gff3
tabix -p gff MergingTestSCN.gff3.gz


```
