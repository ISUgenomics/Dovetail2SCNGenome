#


### run tandem repeat finder for circos
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/31_Synteny/02_738Assembly/02_iadhore/01_circos/01_TRF

#version 4.09
ml trf
trf ../SCNgenome.fasta 2 7 7 80 10 50 500 -h -d -ngs >trf.out  &

less -S 01_TRF/trf.out |awk -v x=0 '{if (substr($1,1,1)=="@") {x=$1} else {print x,$1,$2}}' |sort|uniq|tr " " "\t" |sed 's/@//g'   >TRF.bed

less TRF.bed |awk '{print $3-$2}' |summary.sh
Total:  16,390,271
Count:  87,791
Mean:   186
Median: 69
Min:    24
Max:    46,311

```

### EDTA repeat prediction
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/33_EDTA/EDTA

ln -s ../../10_RepeatModeler/SCNgenome.fasta

ml miniconda3; source activate EDTA; cd /work/GIF/remkv6/Baum/04_Dovetail2Restart/33_EDTA/EDTA; ./EDTA.pl --genome SCNgenome.fasta --threads 16 --overwrite 1 --anno 1 --sensitive 1

Repeat Classes
==============
Total Sequences: 9
Total Length: 156300699 bp
Class                  Count        bpMasked    %masked
=====                  =====        ========     =======
DNA                    --           --           --
    DTA                66646        13237760     8.47%
    DTC                3717         818023       0.52%
    DTH                845          158435       0.10%
    DTM                38907        8374967      5.36%
    DTT                126          44517        0.03%
    Helitron           1791         511988       0.33%
LTR                    --           --           --
    Copia              57           26694        0.02%
    Gypsy              17137        8156160      5.22%
    unknown            50861        11956505     7.65%
MITE                   --           --           --
    DTA                21689        2715415      1.74%
    DTC                376          42636        0.03%
    DTH                5            731          0.00%
    DTM                5701         716083       0.46%
    DTT                3            413          0.00%
                      ---------------------------------
    total interspersed 207861       46760327     29.92%

---------------------------------------------------------
Total                  207861       46760327     29.92%

```

### Get gene counts that are not repetitive
```
#use EDTA repeats
less ../../47_MikadoFinalize/mikado.loci.ancestral.gff3 |awk '$3=="CDS"' |bedtools intersect -v -wo -a - -b ../../10_RepeatModeler/SCNgenome.fasta.out.gff |cut -f 9 |sed 's/\./\t/2' |awk '{print $1}' |sort|uniq |wc -l
19652

#Use EDTA repeats and TRF repeats
less ../../47_MikadoFinalize/mikado.loci.ancestral.gff3 |awk '$3=="CDS"' |bedtools intersect -v -wo -a - -b ../../10_RepeatModeler/SCNgenome.fasta.out.gff |bedtools intersect -v -wo -a - -b ../../31_Synteny/02_738Assembly/02_iadhore/01_circos/TRF.bed |cut -f 9 |sed 's/\./\t/2' |awk '{print $1}' |sort|uniq |wc -l
19451

#Use repeatmodeler repeats, EDTA, and TRF
less ../../47_MikadoFinalize/mikado.loci.ancestral.gff3 |awk '$3=="CDS"' |bedtools intersect -v -wo -a - -b ../../10_RepeatModeler/SCNgenome.fasta.out.gff |bedtools intersect -v -wo -a - -b ../../31_Synteny/02_738Assembly/02_iadhore/01_circos/TRF.bed |bedtools intersect -v -wo -a - -b ../../10_RepeatModeler/SCNgenome.fasta.out.gff  |cut -f 9 |sed 's/\./\t/2' |awk '{print $1}' |sort|uniq |wc -l
19451


```


## Make a tabular file of repetitive mRNA's for feature table

### EDTA tab file
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/33_EDTA/EDTA


#50% coverage of a mrna by a repeat
bedtools intersect -wo -f .5 -b SCNgenome.fasta.EDTA.intact.gff -a ../../47_MikadoFinalize/mikado.loci.ancestral.gff3 |awk '$3=="mRNA"' |cut -f 9,18 |sed 's/;/\t/1' |cut -f 1,3 |sed 's/ID=//1' |awk 'substr($2,1,6)!="Parent"' |sort -u -k1,1 >EDTA.tab

#50% coverage of a mrna by a repeat
bedtools intersect -wo -f .5 -b SCNgenome.fasta.EDTA.intact.gff -a ../../47_MikadoFinalize/mikado.loci.ancestral.gff3 |awk '$3=="mRNA"' |cut -f 9,18 |sed 's/;/\t/1' |cut -f 1,3 |sed 's/ID=//1' |awk 'substr($2,1,6)!="Parent"' |sort -u -k1,1 |wc
586    1172   90338


#Any overlap with an EDTA repeat
bedtools intersect -wo  -b SCNgenome.fasta.EDTA.intact.gff -a ../../47_MikadoFinalize/mikado.loci.ancestral.gff3 |awk '$3=="mRNA"' |cut -f 9,18 |sed 's/;/\t/1' |cut -f 1,3 |sed 's/ID=//1' |awk 'substr($2,1,6)!="Parent"' |sort -u -k1,1 |wc
2620    5240  397235


```

### Repeatmodeler tab file
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/10_RepeatModeler

#50% coverage of an mRNA
bedtools intersect -f .5 -wo  -b SCNgenome.fasta.out.gff -a ../47_MikadoFinalize/mikado.loci.ancestral.gff3 |awk '$3=="mRNA"' |cut -f 9,18 |sed 's/;/\t/1' |cut -f 1,3 |sed 's/ID=//1' |awk 'substr($2,1,6)!="Parent"' |sort -u -k1,1 |awk '{print $1"\t"$3}' |sed 's/"//g' >RepeatModeler.tab

#50% coverage of an mRNA by a repeat
bedtools intersect -f .5 -wo  -b SCNgenome.fasta.out.gff -a ../47_MikadoFinalize/mikado.loci.ancestral.gff3 |awk '$3=="mRNA"' |cut -f 9,18|sed 's/;/\t/1' |cut -f 1,3 |sed 's/ID=//1' |awk 'substr($2,1,6)!="Parent"' |sort -u -k1,1 |wc
2330   11650  153217


#Any overlap with a repeatmodeler repeat   
bedtools intersect -wo  -b SCNgenome.fasta.out.gff -a ../47_MikadoFinalize/mikado.loci.ancestral.gff3 |awk '$3=="mRNA"' |cut -f 9,18 |sed 's/;/\t/1' |cut -f 1,3 |sed 's/ID=//1' |awk 'substr($2,1,6)!="Parent"' |sort -u -k1,1 |wc
15711   78555 1027886



```
