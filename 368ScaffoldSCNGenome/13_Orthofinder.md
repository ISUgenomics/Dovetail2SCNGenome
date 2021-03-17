### Need to find more distant synteny and selection studies for each of the genomes, so need to establish orthologs

```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder

module load orthofinder
orthofinder -t 12 -f 01_Orthofinder/ -n OrthofinderTylenchida


OrthoFinder assigned 137810 genes (80.5% of total) to 13399 orthogroups. Fifty percent of all genes were in orthogroups with 11 or more genes (G50 was 11) and were contained in the largest 3928 orthogroups (O50 was 3928). There were 3050 orthogroups with all species present and 208 of these consisted entirely of single-copy genes.
# This is kind of a neat stat: 208 are conserved at single copy across all 8 species. So these must be the genes that are dosage sensitive, because they are found at single copy in allopolyploid Meloidogyne incognita.
```
### An in-depth look at the gene familes
```


How large are these gene familes?
less Orthogroups.txt |awk '{print NF}' |less
#Some of these gene families are huge!
################################################################################
1781
1371
573
510
415
403
364
299
277
259
243
222
211
204
200
191
190
187
186
etc.
################################################################################

#How many gene familes are larger than 8 ?
less Orthogroups.txt |awk '{print NF}' |awk '$1>8' |wc
   6921    6921   19780



Which gene familes have the most duplicated genes per species?
less Orthogroups.GeneCount.csv
###############################################################################
GeneFamily DovetailSCNMaker4.all.maker.proteins    G.ellingtonae_protein   bursaphelenchus_xylophilus.PRJEA64437.WBPS10.protein    ditylenchus_destructor.PRJNA312427.WBPS10.protein       globodera_pallida.PRJEB123.WBPS10.protein       globodera_rostochiensis.PRJEB13504.WBPS10.protein       meloidogyne_hapla.PRJNA29083.WBPS10.protein     meloidogyne_incognita.PRJEA28837.WBPS10.protein Total
OG0000000       385     837     1       3       475     73      2       4       1780
OG0000001       276     611     4       7       428     38      2       4       1370
OG0000002       312     64      26      57      41      14      17      41      572
OG0000003       354     46      7       2       24      22      16      38      509
OG0000004       279     66      4       3       15      10      36      1       414
OG0000005       191     108     17      2       23      16      34      11      402
OG0000006       108     64      5       5       26      2       24      129     363
OG0000007       0       0       0       0       0       0       102     196     298
OG0000008       123     65      0       12      54      19      2       1       276
OG0000009       75      64      1       1       109     5       2       1       258
OG0000010       0       0       0       0       0       0       103     139     242
OG0000011       116     12      2       3       14      15      23      36      221
OG0000012       140     31      3       5       6       1       24      0       210
OG0000013       79      37      0       1       2       1       9       74      203
OG0000014       121     4       1       13      10      1       22      27      199
OG0000015       92      29      2       1       13      17      4       32      190
OG0000016       108     23      9       4       4       5       6       30      189
OG0000017       95      13      9       13      15      10      9       22      186
OG0000018       107     30      8       3       13      7       8       9       185
OG0000019       92      34      0       0       13      10      9       27      185
OG0000020       88      13      2       4       1       2       29      45      184
OG0000021       24      66      2       10      5       10      10      39      166
OG0000022       72      34      0       0       5       11      4       31      157
OG0000023       18      70      1       2       59      3       1       2       156
OG0000024       56      18      11      9       12      8       12      28      154
OG0000025       38      60      1       1       24      6       11      12      153
OG0000026       99      10      2       0       12      10      8       11      152
OG0000027       112     2       0       2       3       0       19      11      149
OG0000028       99      14      0       0       9       14      3       7       146
OG0000029       20      48      0       1       14      33      0       28      144
OG0000030       123     8       0       0       8       2       1       0       142
OG0000031       132     5       0       0       2       0       1       1       141
OG0000032       56      15      12      11      10      10      11      10      135
OG0000033       104     11      0       0       9       4       1       6       135
OG0000034       10      50      2       7       15      5       24      22      135
OG0000035       59      3       35      12      1       3       5       14      132
OG0000036       67      6       6       18      8       6       10      6       127
OG0000037       0       0       0       0       0       0       15      112     127
OG0000038       69      16      2       1       3       3       6       26      126
OG0000039       30      51      1       3       37      3       0       0       125
OG0000040       55      16      0       26      12      1       3       10      123
OG0000041       18      34      3       26      21      8       2       4       116
OG0000042       70      15      0       3       1       6       1       17      113
OG0000043       43      21      0       4       23      12      2       7       112
OG0000044       0       106     0       0       5       1       0       0       112
OG0000045       28      36      0       0       29      17      0       0       110
OG0000046       87      5       0       1       8       4       3       2       110
OG0000047       30      11      1       5       10      5       29      18      109
OG0000048       94      4       0       0       6       5       0       0       109
###############################################################################
```
# THE ABOVE ANALYSIS IS FLAWED IN THAT MULTIPLE ISOFORMS WERE USED.
# Restart with single isoforms
```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/12_MakerGenesOrthofinder/02_Orthofinder1isoform
less ../01_Orthofinder/DovetailSCNMaker4.all.maker.proteins.fasta |awk '{print $1}' |grep  "\-RA" |sed 's/>//g' |cdbyank ../01_Orthofinder/DovetailSCNMaker4.all.maker.proteins.fasta.cidx >DovetailSCNMaker4.all.maker.proteins.Isoform1Only.fasta
less ../01_Orthofinder/bursaphelenchus_xylophilus.PRJEA64437.WBPS10.protein.fa|awk 'substr($1,14,14)==1' |sed 's/>//g' |awk '{print $1}' |cdbyank ../01_Orthofinder/bursaphelenchus_xylophilus.PRJEA64437.WBPS10.protein.fa.cidx >bursaphelenchus_xylophilus.PRJEA64437.WBPS10.protein.Isoform1Only.fasta
ln -s ../01_Orthofinder/ditylenchus_destructor.PRJNA312427.WBPS10.protein.fa
less ../01_Orthofinder/G.ellingtonae_protein.fasta |grep "t1" |sed 's/>//g' |cdbyank ../01_Orthofinder/G.ellingtonae_protein.fasta.cidx >G.ellingtonae_protein.Isoform1Only.fasta
 less ../01_Orthofinder/G.ellingtonae_protein.fasta |grep "t1" |sed 's/>//g' |cdbyank ../01_Orthofinder/G.ellingtonae_protein.fasta.cidx >G.ellingtonae_protein.Isoform1Only.fasta
 ln -s ../01_Orthofinder/globodera_pallida.PRJEB123.WBPS10.protein.fa
 less ../01_Orthofinder/globodera_rostochiensis.PRJEB13504.WBPS10.protein.fa |awk 'substr($1,14,15)=="t1"' |awk '{print $1}' |sed 's/>//g' |cdbyank ../01_Orthofinder/globodera_rostochiensis.PRJEB13504.WBPS10.protein.fa.cidx >globodera_rostochiensis.PRJEB13504.WBPS10.protein.Isoform1Only.fasta
 ln -s ../01_Orthofinder/meloidogyne_incognita.PRJEA28837.WBPS10.protein.fa

 module load orthofinder
 orthofinder -t 12 -f 02_Orthofinder1isoform -n OrthofinderTylenchida

 ############################################################################
 OrthoFinder assigned 111455 genes (79.9% of total) to 12456 orthogroups. Fifty percent of all genes were in orthogroups with 9 or more genes (G50 was 9) and were contained in the largest 4261 orthogroups (O50 was 4261). There were 2985 orthogroups with all species present and 356 of these consisted entirely of single-copy genes.
 ############################################################################

#############################################################################
DovetailSCNMaker4.all.maker.proteins.Isoform1Only       G.ellingtonae_protein.Isoform1Only      bursaphelenchus_xylophilus.PRJEA64437.WBPS10.protein.Isoform1Only       ditylenchus_destructor.PRJNA312427.WBPS10.protein       globodera_pallida.PRJEB123.WBPS10.protein       globodera_rostochiensis.PRJEB13504.WBPS10.protein.Isoform1Only  meloidogyne_hapla.PRJNA29083.WBPS10.protein     meloidogyne_incognita.PRJEA28837.WBPS10.protein Total
OG0000000       180     795     2       5       478     72      2       5       1539
OG0000001       135     597     5       3       429     39      2       4       1214
OG0000002       96      108     17      2       21      15      52      9       320
OG0000003       44      63      26      65      41      13      17      42      311
OG0000004       42      68      5       5       27      2       24      131     304
OG0000005       0       0       0       0       0       0       94      186     280
OG0000006       0       0       0       0       0       0       100     138     238
OG0000007       18      63      1       2       108     6       3       1       202
OG0000008       63      62      5       5       13      4       35      0       187
OG0000009       37      60      0       12      54      19      2       1       185
OG0000010       27      51      1       3       15      34      6       44      181
OG0000011       39      38      0       1       2       1       9       74      164
OG0000012       19      61      1       23      26      6       12      15      163
OG0000013       9       1       1       12      1       1       31      105     161
OG0000014       2       78      1       2       70      3       1       2       159
OG0000015       12      66      2       10      5       10      10      40      155
OG0000016       73      28      0       0       12      14      7       18      152
OG0000017       17      46      2       7       15      5       24      22      138
OG0000018       40      12      0       2       13      15      19      35      136
OG0000019       34      29      1       1       12      17      4       36      134
OG0000020       36      14      2       4       1       2       27      45      131
OG0000021       55      11      0       10      14      9       7       22      128
OG0000022       0       0       0       0       0       0       15      112     127
OG0000023       28      17      0       4       16      12      29      20      126
OG0000024       28      49      2       3       40      4       0       0       126
OG0000025       1       109     0       0       5       1       0       0       116
OG0000026       16      18      3       9       12      9       12      28      107
OG0000027       4       9       38      24      11      10      3       6       105
OG0000028       9       33      3       23      21      8       2       4       103
OG0000029       62      11      0       0       8       12      3       7       103
OG0000030       21      9       1       7       1       4       14      43      100
#############################################################################
```
