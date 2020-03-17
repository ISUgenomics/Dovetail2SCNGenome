# Need a comparison with the draft genome as to how many tandem duplications are maintained in the genome

```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/37_redTandem

ln -s ../10_RepeatModeler/SCNgenome.fasta.masked

#formats fasta for redtandem
sed -e 's/a/A/g' -e 's/c/C/g' -e 's/g/G/g' -e 's/t/T/g' SCNgenome.fasta.masked | bioawk -c fastx '{print ">"$name"_1-"length($seq)" \n"$seq}' SCNgenome.fasta.masked  |fold -79 >SCNGenomeFormatedMasked.fasta


#run but leave the temp files so I can troubleshoot
 for f in *part*; do echo " ml use /shared/software/GIF/modules ; ml GIF2/redtandem;  perl /shared/software/GIF/programs/redtandem/ReDtandem/ReDtandem.pl --species scn --noclean --dnafile "$f;done >redtandem.sh


#5/9 chromosomes did not finish, liekly due to unmasked tandem repeats

#try again with masked tandems

ml bedtools2

bedtools maskfasta  -fi SCNgenome.fasta.masked -bed TRF.bed -fo TRFMaskedSCNgenome.fasta.masked
sed -e 's/a/A/g' -e 's/c/C/g' -e 's/g/G/g' -e 's/t/T/g' TRFMaskedSCNgenome.fasta.masked|bioawk -c fastx '{print ">"$name"_1-"length($seq)" \n"$seq}' TRFMaskedSCNgenome.fasta.masked |fold -79 >TRFMaskedSCNgenomeFormatted.fasta.masked

fasta-splitter.pl --n-parts 9 TRFMaskedSCNgenomeFormatted.fasta.masked

for f in *part*;do echo " ml use /shared/software/GIF/modules ; ml GIF2/redtandem;  perl /shared/software/GIF/programs/redtandem/ReDtandem/ReDtandem.pl --species scn --noclean --dnafile " $f;done >redtandem.sh


```
