##  Need to create stats for each gene so we know which are likely the most important for studying effectors
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/40_PrepareGeneCharacteristics

#mrna locations
awk '$3=="mRNA"' ../25_AnnotateGenes/07_NewGenes/OrderedSCNGenePredictions.gff3 |sed 's/ID=//g' |sed 's/;/\t/1' |awk '{print $9"\t"$1":"$4"-"$5}' >mRNALocations.tab

#mrna sequences in tab format
awk '{print $1}' ../25_AnnotateGenes/07_NewGenes/OrderedSCNGenePredictions_VHEJtranscripts.fasta |tr "\n" "\t" |sed 's/>/\n/g' |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |sed '/^$/d' |awk '{print $1"\t"$2}' >mRNASeqs.tab


#protein sequences in tab format
awk '{print $1}' ../25_AnnotateGenes/07_NewGenes/OrderedSCNGenePredictionsVHEJ_proteins.fasta |tr "\n" "\t" |sed 's/>/\n/g' |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\t/g' |sed '/^$/d' |awk '{print $1"\t"$2}'>protSeqs.tab

#get all of the expression tabs
for f in ../39_DifferentialExpression/02_mRNA/*Allgenes.tab;do ln -s $f;done

#all functional annotations
ln -s  ../25_AnnotateGenes/01_Interpro/interproAnnot.tab1 interproAnnot.tab
cat ../25_AnnotateGenes/06_Combine/NRBlastp.tab  AllmRNAsNA |sort -u -k1,1 >NRBlastpAllmRNA.tab
cat ../25_AnnotateGenes/06_Combine/NTBlastn.tab  AllmRNAsNA |sort -u -k1,1 >NTBlastnAllmRNA.tab
cat ../25_AnnotateGenes/06_Combine/UPBlastp.tab  AllmRNAsNA |sort -u -k1,1 >UPBlastpllmRNA.tab
cat ../25_AnnotateGenes/06_Combine/UPBlastx.tab  AllmRNAsNA |sort -u -k1,1 >UPBlastxllmRNA.tab


#repeat annotations
cat ../10_RepeatModeler/RepeatModeler.tab  AllmRNAsNA |sort -u -k1,1 >RepeatmodelerAllmRNA.tab
cat ../33_EDTA/EDTA/EDTA.tab  AllmRNAsNA |sort -u -k1,1 >EDTAAllmRNA.tab


# effector annotations
cat ../29_Effectors/GmapEffectorsGenes.tab AllmRNAsNA |sort -u -k1,1 >GmapEffectorsAllmRNA.tab
cat ../29_Effectors/01_Diamond/NamedDiamondEffectors.tab  AllmRNAsNA |sort -u -k1,1 >DiamondEffectorsAllmRNA.tab


#merops annotations
less ../41_Merops/Merops.tab |cut -f 1,2 |cat - AllmRNAsNA |sort -u -k1,1 >Merops1.tab
less ../41_Merops/Merops.tab |cut -f 1,3 |cat - AllmRNAsNA |sort -u -k1,1>Merops2.tab

#NLSTRADAMUS
grep -v "#" ../42_NLSPrediction/NLSsignals|cat - AllmRNAsNA |sort -u -k1,1 >NLSsignals.tab

##Tandem repeat finder intersect for CDS
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/02_ProteinRepeats
bedtools intersect -f .3 -wo -b TRF.bed -a OrderedSCNGenePredictions.gff3 |awk '$3=="CDS"' |cut -f 9 |sed 's/;/\t/g' |sed 's/Parent=//g' |cut -f 1 |sort|uniq>TRmRNAs.list



#Compile all annotations
#sort all to be put in the same order
for f in *tab; do sort -k1,1V $f >${f%.*}sorted;done
ls -1 *sorted |tr "\n" " " |sed 's/^/paste /g' |sed 's/$/ >AllFeatures.tab1/g' >paste.sh
sh paste.sh
less AllFeatures.tab1 |cut -f 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66 >DetrashedAllFeatures.tab1


ls -1 *tab |tr "\n" "\t" |awk '{print "mRNA_Name\t"$0}' |sed 's/\.tab//g' |awk '{print $0}' >AllFeaturesHeader.tab1          
paste AllFeaturesHeader.tab1 DetrashedAllFeatures.tab1 >FinishedAllFeatures.tab




```
