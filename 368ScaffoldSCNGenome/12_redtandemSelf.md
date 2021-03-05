# Want to get a handle on a second predictor for tandem duplications in the genome
```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/14_redtandemSelf
bioawk -c fastx '{print "SCN"$name"_1-"length($seq)"\n"$seq}' renamedDovetailSCNMaker4.genome.fasta >SCNGenomeMod.fasta

#had to add a module purge before loading GIF/redtandem
perl /shared/software/GIF/programs/redtandem/ReDtandem/ReDtandem.pl --species SCN --noclean --dnafile DovetailSCNMaker4.genome.fasta

```
