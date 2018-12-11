# Need to call repeats on the new dovetail2 genome

### run repeatmodeler
```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/20_Dovetail2RepeatModeler


module use /work/GIF/software/modules
module load GIF2/repeatmodeler
BuildDatabase -name scnDovetail2 -engine ncbi -dir /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/20_Dovetail2RepeatModeler nematode_sp_19Jul2018_IbtP1.fasta
RepeatModeler -database  scnDovetail2 -engine ncbi -pa 12

```

### concatenate repeatmodeler and repeatexplorer fasta, mask with fasta

```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/20_Dovetail2RepeatModeler/01_Repeatmasker

ln -s ../nematode_sp_19Jul2018_IbtP1.fasta
ln -s ../RM_19764.FriAug31522082018/consensi.fa.classified
ln -s  ../../03_RepeatModelerMasker/ModelerExplorerClassified.fasta

module load cdbfasta
cdbfasta ModelerExplorerClassified.fasta
grep ">CL" ModelerExplorerClassified.fasta |sed 's/>//g' |cdbyank ModelerExplorerClassified.fasta.cidx |cat - consensi.fa.classified >RepeatModelerExplorerClassified.fasta

printf "module load repeatmasker/4.0.7\nRepeatMasker -pa 12 -gff -lib RepeatModelerExplorerClassified.fasta nematode_sp_19Jul2018_IbtP1.fasta\n"
```
