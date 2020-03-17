##  Need to create stats for each gene so we know which are likely the most important for studying effectors
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/40_PrepareGeneCharacteristics

less ../39_DifferentialExpression/AllGlandvsAllWholeWorm |sed 's/,/\t/g' >AllGlandvsAllWholeWorm.tab

less ../39_DifferentialExpression/AllGlandvsAllWholeWorm |sed 's/,/\t/g' |awk 'NR<5281 && NR>1 && $3>0{print $1}' >AllGlandvsAllWholeWormListUp.tab
less ../39_DifferentialExpression/PA3GlandvsMM10Gland |sed 's/,/\t/g' |awk 'NR<6 && NR>1 && $3>0{print $1}' >PA3GlandvsMM10GlandListUp.tab
less ../39_DifferentialExpression/pPA3CvsPA3Gland |sed 's/,/\t/g' |awk 'NR<573 && NR>1 && $3>0{print $1}' >pPA3CvsPA3GlandListUp.tab
less ../39_DifferentialExpression/ppPA3vspPA3C |sed 's/,/\t/g' |awk 'NR<16 && NR>1 && $3>0{print $1}' >ppPA3vspPA3CListUp.tab
less ../39_DifferentialExpression/ppPA3vspPA3IC |sed 's/,/\t/g' |awk 'NR<7461 && NR>1 && $3>0{print $1}' >ppPA3vspPA3ICListUp.tab
less ../39_DifferentialExpression/pPA3CvspPA3IC |sed 's/,/\t/g' |awk 'NR<7284 && NR>1 && $3>0{print $1}' >pPA3CvspPA3ICListUp.tab
less ../39_DifferentialExpression/2PA3vs3MM10GlandGland |sed 's/,/\t/g' |awk 'NR<13 && NR>1 && $3>0{print $1}' >2PA3vs3MM10GlandListUp.tab
less ../39_DifferentialExpression/2PA3vs2MM10GlandGland |sed 's/,/\t/g' |awk 'NR<51 && NR>1 && $3>0{print $1}' >2PA3vs2MM10GlandListUp.tab

less ../39_DifferentialExpression/AllGlandvsAllWholeWorm |sed 's/,/\t/g' |awk 'NR<5281 && NR>1 && $3<0{print $1}' >AllGlandvsAllWholeWormListDown.tab
less ../39_DifferentialExpression/PA3GlandvsMM10Gland |sed 's/,/\t/g' |awk 'NR<6 && NR>1 && $3<0{print $1}' >PA3GlandvsMM10GlandListDown.tab
less ../39_DifferentialExpression/pPA3CvsPA3Gland |sed 's/,/\t/g' |awk 'NR<573 && NR>1 && $3<0{print $1}' >pPA3CvsPA3GlandListDown.tab
less ../39_DifferentialExpression/ppPA3vspPA3C |sed 's/,/\t/g' |awk 'NR<16 && NR>1 && $3<0{print $1}' >ppPA3vspPA3CListDown.tab
less ../39_DifferentialExpression/ppPA3vspPA3IC |sed 's/,/\t/g' |awk 'NR<7461 && NR>1 && $3<0{print $1}' >ppPA3vspPA3ICListDown.tab
less ../39_DifferentialExpression/pPA3CvspPA3IC |sed 's/,/\t/g' |awk 'NR<7284 && NR>1 && $3<0{print $1}' >pPA3CvspPA3ICListDown.tab
less ../39_DifferentialExpression/2PA3vs3MM10GlandGland |sed 's/,/\t/g' |awk 'NR<13 && NR>1 && $3<0{print $1}' >2PA3vs3MM10GlandListDown.tab
less ../39_DifferentialExpression/2PA3vs2MM10GlandGland |sed 's/,/\t/g' |awk 'NR<51 && NR>1 && $3<0{print $1}' >2PA3vs2MM10GlandListDown.tab





ln -s ../29_Effectors/GmappedEffectorsGene.list
ln -s ../25_AnnotateGenes/mikado.loci.gff3
awk '$3=="gene"' mikado.loci.gff3 |awk '{print $9,$1":"$4"-"$5}' |tr " " "\t" >GeneLocationsList.tab

less GeneLocationsList.tab |cut -f 2 |samtools faidx SCNgenome.fasta -r - |tr "\n" "#" |sed 's/>/\n/g' |sed 's/#/\t/1' |sed 's/#//g' >geneSeqslist.tab
less mikado_proteins.fasta |awk '{print $1}' |tr "\n" "#" |sed 's/>/\n/g' |sed 's/#/\t/1' |sed 's/#//g' >proteinSeqslist.tab
less mikado_transcripts.fasta |awk '{print $1}' |tr "\n" "#" |sed 's/>/\n/g' |sed 's/#/\t/1' |sed 's/#//g' >trascriptSeqslist.tab

ln -s ../25_AnnotateGenes/06_Combine/interproAnnot.tab1
ln -s ../25_AnnotateGenes/06_Combine/NRBlastp.tab
ln -s ../25_AnnotateGenes/06_Combine/NTBlastn.tab
ln -s ../25_AnnotateGenes/06_Combine/UPBlastp.tab
ln -s ../25_AnnotateGenes/06_Combine/UPBlastx.tab

```
