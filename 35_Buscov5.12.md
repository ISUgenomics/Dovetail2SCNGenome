# Rerun all of the buscos with v5

### TN10 DRAFT genome
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/01_738
ln -s ../../01_CondoPorts/31_Synteny/02_738Assembly/genome738sl.polished.mitofixed.fasta
echo "ml miniconda3; source activate busco5_env ; busco -i genome738sl.polished.mitofixed.fasta -o 738Busco -m geno --auto-lineage-euk -c 35 -f " >busco.sh

        --------------------------------------------------
        |Results from generic domain eukaryota_odb10      |
        --------------------------------------------------
        |C:60.4%[S:51.0%,D:9.4%],F:20.0%,M:19.6%,n:255    |
        |154    Complete BUSCOs (C)                       |
        |130    Complete and single-copy BUSCOs (S)       |
        |24     Complete and duplicated BUSCOs (D)        |
        |51     Fragmented BUSCOs (F)                     |
        |50     Missing BUSCOs (M)                        |
        |255    Total BUSCO groups searched               |
        --------------------------------------------------

#run with --long and --augustus
        --------------------------------------------------
        |Results from dataset nematoda_odb10              |
        --------------------------------------------------
        |C:55.1%[S:47.4%,D:7.7%],F:2.0%,M:42.9%,n:3131    |
        |1725   Complete BUSCOs (C)                       |
        |1485   Complete and single-copy BUSCOs (S)       |
        |240    Complete and duplicated BUSCOs (D)        |
        |62     Fragmented BUSCOs (F)                     |
        |1344   Missing BUSCOs (M)                        |
        |3131   Total BUSCO groups searched               |
        --------------------------------------------------


#How many are actually missing from nematoda odb10 using both protein and genome datasets
cat 738BuscoLong/run_nematoda_odb10/missing_busco_list.tsv 01_proteins/738Buscoprot/run_nematoda_odb10/missing_busco_list.tsv |sort|uniq -c |awk '$1==2' |wc
    940    1902   17993




```

### TN10 Draft protein
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/01_738/01_proteins
wget https://scnbase.org/files/download/finalaugustus.pep_.fasta
echo "ml miniconda3; source activate busco5_env ; busco -i finalaugustus.pep_.fasta -o 738Buscoprot -m prot --auto-lineage-euk -c 35 -f " >busco.sh
--------------------------------------------------
|Results from generic domain eukaryota_odb10      |
--------------------------------------------------
|C:89.0%[S:68.2%,D:20.8%],F:3.1%,M:7.9%,n:255     |
|227    Complete BUSCOs (C)                       |
|174    Complete and single-copy BUSCOs (S)       |
|53     Complete and duplicated BUSCOs (D)        |
|8      Fragmented BUSCOs (F)                     |
|20     Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------

--------------------------------------------------
|Results from dataset nematoda_odb10              |
--------------------------------------------------
|C:66.5%[S:47.7%,D:18.8%],F:1.6%,M:31.9%,n:3131   |
|2084   Complete BUSCOs (C)                       |
|1494   Complete and single-copy BUSCOs (S)       |
|590    Complete and duplicated BUSCOs (D)        |
|49     Fragmented BUSCOs (F)                     |
|998    Missing BUSCOs (M)                        |
|3131   Total BUSCO groups searched               |
--------------------------------------------------

#What percentage would be considered complete if we eliminate the 688 busco genes missing from all species from the total
3131-666= 2465
940-666 = 274

#Total percent complete (2465- 274 = 2191)/2465 = 88.9% not missing
```

### TN10 Pseudo genome
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/02_tn10pseudo
ln -s ../../49_RenameChromosomes/01_Transfer2Box/SCNgenome.fasta
echo "ml miniconda3; source activate busco5_env ; busco -i SCNgenome.fasta -o TN10Buscogeno -m geno --auto-lineage-euk -c 35 -f " >busco.sh

        --------------------------------------------------
        |Results from generic domain eukaryota_odb10      |
        --------------------------------------------------
        |C:60.0%[S:56.5%,D:3.5%],F:22.0%,M:18.0%,n:255    |
        |153    Complete BUSCOs (C)                       |
        |144    Complete and single-copy BUSCOs (S)       |
        |9      Complete and duplicated BUSCOs (D)        |
        |56     Fragmented BUSCOs (F)                     |
        |46     Missing BUSCOs (M)                        |
        |255    Total BUSCO groups searched               |
        --------------------------------------------------
#run with --long and --augustus
        --------------------------------------------------
        |Results from dataset nematoda_odb10              |
        --------------------------------------------------
        |C:55.6%[S:52.5%,D:3.1%],F:2.3%,M:42.1%,n:3131    |
        |1740   Complete BUSCOs (C)                       |
        |1644   Complete and single-copy BUSCOs (S)       |
        |96     Complete and duplicated BUSCOs (D)        |
        |72     Fragmented BUSCOs (F)                     |
        |1319   Missing BUSCOs (M)                        |
        |3131   Total BUSCO groups searched               |
        --------------------------------------------------


#How many are actually missing from nematoda odb10 using both protein and genome datasets

cat TN10BuscogenoLong/run_nematoda_odb10/missing_busco_list.tsv 01_protein/TN10Psuedoprot/run_nematoda_odb10/missing_busco_list.tsv |sort|uniq -c |awk '$1==2' |wc
    965




```
### TN10 Pseudo protein
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/02_tn10pseudo/01_protein
ln -s ../../../49_RenameChromosomes/mikado.loci.ancestralVHEJ_proteins.fasta
echo "ml miniconda3; source activate busco5_env ; busco -i mikado.loci.ancestralVHEJ_proteins.fasta -o TN10Psuedoprot -m prot --auto-lineage-euk -c 35 -f " >busco.sh


        --------------------------------------------------
        |Results from generic domain eukaryota_odb10      |
        --------------------------------------------------
        |C:78.0%[S:69.0%,D:9.0%],F:7.8%,M:14.2%,n:255     |
        |199    Complete BUSCOs (C)                       |
        |176    Complete and single-copy BUSCOs (S)       |
        |23     Complete and duplicated BUSCOs (D)        |
        |20     Fragmented BUSCOs (F)                     |
        |36     Missing BUSCOs (M)                        |
        |255    Total BUSCO groups searched               |
        --------------------------------------------------

        --------------------------------------------------
        |Results from dataset metazoa_odb10               |
        --------------------------------------------------
        |C:62.0%[S:52.8%,D:9.2%],F:6.5%,M:31.5%,n:954     |
        |592    Complete BUSCOs (C)                       |
        |504    Complete and single-copy BUSCOs (S)       |
        |88     Complete and duplicated BUSCOs (D)        |
        |62     Fragmented BUSCOs (F)                     |
        |300    Missing BUSCOs (M)                        |
        |954    Total BUSCO groups searched               |
        --------------------------------------------------
        --------------------------------------------------
|Results from dataset nematoda_odb10              |
--------------------------------------------------
|C:60.2%[S:51.5%,D:8.7%],F:1.5%,M:38.3%,n:3131    |
|1884   Complete BUSCOs (C)                       |
|1612   Complete and single-copy BUSCOs (S)       |
|272    Complete and duplicated BUSCOs (D)        |
|47     Fragmented BUSCOs (F)                     |
|1200   Missing BUSCOs (M)                        |
|3131   Total BUSCO groups searched               |
--------------------------------------------------

#What percentage would be considered complete if we eliminate the 688 busco genes missing from all species from the total
3131-666= 2465
965-666 =299

#Total percent complete (2465 - 299=2166)/2465 = 87.9% not missing


```
### X12 genome
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/03_X12
ln -s ../../01_CondoPorts/31_Synteny/01_X12/X12SCN_genome.fa
echo "ml miniconda3; source activate busco5_env ; busco -i X12SCN_genome.fa -o X12Buscogeno -m geno --auto-lineage-euk -c 35 -f " >busco.sh
--------------------------------------------------
|Results from generic domain eukaryota_odb10      |
--------------------------------------------------
|C:57.2%[S:54.1%,D:3.1%],F:17.3%,M:25.5%,n:255    |
|146    Complete BUSCOs (C)                       |
|138    Complete and single-copy BUSCOs (S)       |
|8      Complete and duplicated BUSCOs (D)        |
|44     Fragmented BUSCOs (F)                     |
|65     Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------
#run with --long and --augustus
--------------------------------------------------
|Results from dataset nematoda_odb10              |
--------------------------------------------------
|C:52.2%[S:49.0%,D:3.2%],F:1.8%,M:46.0%,n:3131    |
|1634   Complete BUSCOs (C)                       |
|1535   Complete and single-copy BUSCOs (S)       |
|99     Complete and duplicated BUSCOs (D)        |
|57     Fragmented BUSCOs (F)                     |
|1440   Missing BUSCOs (M)                        |
|3131   Total BUSCO groups searched               |
--------------------------------------------------


#How many are actually missing from nematoda odb10 using both protein and genome datasets
cat X12BuscogenoLong/run_nematoda_odb10/missing_busco_list.tsv 01_protein/*/run_nematoda_odb10/missing_busco_list.tsv |sort|uniq -c |awk '$1==2' |wc
   1142    2306   21813


```
### X12 protein
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/03_X12/01_protein
ln -s ../../../01_CondoPorts/31_Synteny/01_X12/pasa2.longest.filter.pep
echo "ml miniconda3; source activate busco5_env ; busco -i pasa2.longest.filter.pep -o X12Buscoprot -m prot --auto-lineage-euk -c 35 -f " >busco.sh

        --------------------------------------------------
        |Results from generic domain eukaryota_odb10      |
        --------------------------------------------------
        |C:64.7%[S:58.8%,D:5.9%],F:13.7%,M:21.6%,n:255    |
        |165    Complete BUSCOs (C)                       |
        |150    Complete and single-copy BUSCOs (S)       |
        |15     Complete and duplicated BUSCOs (D)        |
        |35     Fragmented BUSCOs (F)                     |
        |55     Missing BUSCOs (M)                        |
        |255    Total BUSCO groups searched               |
        --------------------------------------------------
        --------------------------------------------------
        |Results from dataset nematoda_odb10              |
        --------------------------------------------------
        |C:53.9%[S:49.4%,D:4.5%],F:2.2%,M:43.9%,n:3131    |
        |1689   Complete BUSCOs (C)                       |
        |1547   Complete and single-copy BUSCOs (S)       |
        |142    Complete and duplicated BUSCOs (D)        |
        |70     Fragmented BUSCOs (F)                     |
        |1372   Missing BUSCOs (M)                        |
        |3131   Total BUSCO groups searched               |
        --------------------------------------------------
#What percentage would be considered complete if we eliminate the 688 busco genes missing from all species from the total
3131-666= 2465
1142-666 = 476

#Total percent complete (2465-476=1989)/2465 = 80.7% not missing
```
### Globodera pallida genome
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/04_GlobPal
ln -s ../../../06_Orthology/01_genomicResources/globodera_pallida.PRJEB123.WBPS15.genomic.fa
echo "ml miniconda3; source activate busco5_env ; busco -i globodera_pallida.PRJEB123.WBPS15.genomic.fa -o GpalBuscogeno -m geno --auto-lineage-euk -c 35 -f " >busco.sh

        --------------------------------------------------
        |Results from dataset eukaryota_odb10             |
        --------------------------------------------------
        |C:48.2%[S:43.9%,D:4.3%],F:22.4%,M:29.4%,n:255    |
        |123    Complete BUSCOs (C)                       |
        |112    Complete and single-copy BUSCOs (S)       |
        |11     Complete and duplicated BUSCOs (D)        |
        |57     Fragmented BUSCOs (F)                     |
        |75     Missing BUSCOs (M)                        |
        |255    Total BUSCO groups searched               |
        --------------------------------------------------
#run with --long and --augustus
        --------------------------------------------------
        |Results from dataset nematoda_odb10              |
        --------------------------------------------------
         |C:46.7%[S:43.5%,D:3.2%],F:2.5%,M:50.8%,n:3131    |
         |1464   Complete BUSCOs (C)                       |
         |1363   Complete and single-copy BUSCOs (S)       |
         |101    Complete and duplicated BUSCOs (D)        |
         |77     Fragmented BUSCOs (F)                     |
         |1590   Missing BUSCOs (M)                        |
         |3131   Total BUSCO groups searched               |
         --------------------------------------------------

        #How many are actually missing from nematoda odb10 using both protein and genome datasets
        cat GpalBuscogenoLong/run_nematoda_odb10/missing_busco_list.tsv 01_protein/*/run_nematoda_odb10/missing_busco_list.tsv |sort|uniq -c |awk '$1==2' |wc
           1377    2776   26259



```
### Globodera pallida protein
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/04_GlobPal/01_protein
ln -s ../../../../06_Orthology/01_genomicResources/globodera_pallida.PRJEB123.WBPS15.protein.fa
echo "ml miniconda3; source activate busco5_env ; busco -i globodera_pallida.PRJEB123.WBPS15.protein.fa -o GpalBuscoprot -m prot --auto-lineage-euk -c 35 -f " >busco.sh


        --------------------------------------------------
        |Results from dataset eukaryota_odb10             |
        --------------------------------------------------
        |C:62.3%[S:58.8%,D:3.5%],F:14.1%,M:23.6%,n:255    |
        |159    Complete BUSCOs (C)                       |
        |150    Complete and single-copy BUSCOs (S)       |
        |9      Complete and duplicated BUSCOs (D)        |
        |36     Fragmented BUSCOs (F)                     |
        |60     Missing BUSCOs (M)                        |
        |255    Total BUSCO groups searched               |
        --------------------------------------------------
        --------------------------------------------------
        |Results from dataset nematoda_odb10              |
        --------------------------------------------------
        |C:45.6%[S:41.9%,D:3.7%],F:2.8%,M:51.6%,n:3131    |
        |1427   Complete BUSCOs (C)                       |
        |1312   Complete and single-copy BUSCOs (S)       |
        |115    Complete and duplicated BUSCOs (D)        |
        |87     Fragmented BUSCOs (F)                     |
        |1617   Missing BUSCOs (M)                        |
        |3131   Total BUSCO groups searched               |
        --------------------------------------------------
#What percentage would be considered complete if we eliminate the 688 busco genes missing from all species from the total
3131-666 =2465
1377-666 = 711



#Total percent complete (2465- 711=1754)/2465 = 71.2% not missing

```
### Globodera rostochiensis genome
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/05_GlobRos
ln -s ../../../06_Orthology/01_genomicResources/globodera_rostochiensis.PRJEB13504.WBPS15.genomic.fa
echo "ml miniconda3; source activate busco5_env ; busco -i globodera_rostochiensis.PRJEB13504.WBPS15.genomic.fa -o GrosBuscogeno -m geno --auto-lineage-euk -c 35 -f " >busco.sh
--------------------------------------------------
|Results from generic domain eukaryota_odb10      |
--------------------------------------------------
|C:63.2%[S:62.4%,D:0.8%],F:18.4%,M:18.4%,n:255    |
|161    Complete BUSCOs (C)                       |
|159    Complete and single-copy BUSCOs (S)       |
|2      Complete and duplicated BUSCOs (D)        |
|47     Fragmented BUSCOs (F)                     |
|47     Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------
#run with --long and --augustus
--------------------------------------------------
 |Results from dataset nematoda_odb10              |
 --------------------------------------------------
 |C:60.0%[S:58.9%,D:1.1%],F:1.9%,M:38.1%,n:3131    |
 |1878   Complete BUSCOs (C)                       |
 |1845   Complete and single-copy BUSCOs (S)       |
 |33     Complete and duplicated BUSCOs (D)        |
 |59     Fragmented BUSCOs (F)                     |
 |1194   Missing BUSCOs (M)                        |
 |3131   Total BUSCO groups searched               |
 --------------------------------------------------


#How many are actually missing from nematoda odb10 using both protein and genome datasets
cat GrosBuscogenoLong/run_nematoda_odb10/missing_busco_list.tsv 01_protein/*/run_nematoda_odb10/missing_busco_list.tsv |sort|uniq -c |awk '$1==2' |wc
    881    1784   16890




```
### Globodera rostochiensis protein
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/05_GlobRos/01_protein
ln -s ../../../../06_Orthology/01_genomicResources/globodera_rostochiensis.PRJEB13504.WBPS15.protein.fa
echo "ml miniconda3; source activate busco5_env ; busco -i globodera_rostochiensis.PRJEB13504.WBPS15.protein.fa -o GrosBuscoprot -m prot --auto-lineage-euk -c 35 -f " >busco.sh

#would not finish properly for nematoda on autolineage
 C:87.5%[S:85.9%,D:1.6%],F:6.7%,M:5.8%,n:255

 --------------------------------------------------
|Results from dataset nematoda_odb10              |
--------------------------------------------------
|C:66.7%[S:64.3%,D:2.4%],F:2.1%,M:31.2%,n:3131    |
|2088   Complete BUSCOs (C)                       |
|2013   Complete and single-copy BUSCOs (S)       |
|75     Complete and duplicated BUSCOs (D)        |
|66     Fragmented BUSCOs (F)                     |
|977    Missing BUSCOs (M)                        |
|3131   Total BUSCO groups searched               |
--------------------------------------------------


#What percentage would be considered complete if we eliminate the 688 busco genes missing from all species from the total
3131-666=2465
881-666 = 215



#Total percent complete (2465- 215=2,250)/2465 = 91.3% not missing
```
### Globodera ellingtonae genome
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/06_Globell
ln -s ../../../06_Orthology/01_genomicResources/globodera_ellingtonae.GCA_001723225.1.genomic.fa
echo "ml miniconda3; source activate busco5_env ; busco -i globodera_ellingtonae.GCA_001723225.1.genomic.fa -o GellBuscogeno -m geno --auto-lineage-euk -c 35 -f " >busco.sh
--------------------------------------------------
|Results from generic domain eukaryota_odb10      |
--------------------------------------------------
|C:59.6%[S:59.6%,D:0.0%],F:18.4%,M:22.0%,n:255    |
|152    Complete BUSCOs (C)                       |
|152    Complete and single-copy BUSCOs (S)       |
|0      Complete and duplicated BUSCOs (D)        |
|47     Fragmented BUSCOs (F)                     |
|56     Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------
#run with --long and --augustus
--------------------------------------------------
|Results from dataset nematoda_odb10              |
--------------------------------------------------
|C:54.4%[S:52.8%,D:1.6%],F:2.4%,M:43.2%,n:3131    |
|1701   Complete BUSCOs (C)                       |
|1652   Complete and single-copy BUSCOs (S)       |
|49     Complete and duplicated BUSCOs (D)        |
|74     Fragmented BUSCOs (F)                     |
|1356   Missing BUSCOs (M)                        |
|3131   Total BUSCO groups searched               |
--------------------------------------------------


#How many are actually missing from nematoda odb10 using both protein and genome datasets
cat GellBuscogenoLong/run_nematoda_odb10/missing_busco_list.tsv 01_protein/*/run_nematoda_odb10/missing_busco_list.tsv |sort|uniq -c |awk '$1==2' |wc
    947    1916   18137



```

### Globodera ellingtonae protein --- my old braker annotation is all I can find
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/06_Globell/01_protein
ln -s ../../../../06_Orthology/01_genomicResources/globodera_ellingtonae.GCA_001723225.1.protein.fa
echo "ml miniconda3; source activate busco5_env ; busco -i globodera_ellingtonae.GCA_001723225.1.protein.fa -o GellBuscoprot -m prot --auto-lineage-euk -c 35 -f " >busco.sh

#doesnt finish properly for nematode for some reason
eukaryota_odb10 -- C:87.9%[S:80.4%,D:7.5%],F:5.9%,M:6.2%,n:255

--------------------------------------------------
|Results from dataset nematoda_odb10              |
--------------------------------------------------
|C:65.4%[S:56.5%,D:8.9%],F:2.0%,M:32.6%,n:3131    |
|2048   Complete BUSCOs (C)                       |
|1770   Complete and single-copy BUSCOs (S)       |
|278    Complete and duplicated BUSCOs (D)        |
|64     Fragmented BUSCOs (F)                     |
|1019   Missing BUSCOs (M)                        |
|3131   Total BUSCO groups searched               |
--------------------------------------------------

#What percentage would be considered complete if we eliminate the 688 busco genes missing from all species from the total
3131-666=2465
947-666 = 281



#Total percent complete (2465- 281=2184)/2465 = 88.6% not missing



```
### Meloidogyne hapla genome
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/08_Melhap
ln -s ../../../06_Orthology/01_genomicResources/meloidogyne_hapla.PRJNA29083.WBPS15.genomic.fa
echo "ml miniconda3; source activate busco5_env ; busco -i meloidogyne_hapla.PRJNA29083.WBPS15.genomic.fa -o MhapBuscogeno -m geno --auto-lineage-euk -c 35 -f " >busco.sh

--------------------------------------------------
|Results from generic domain eukaryota_odb10      |
--------------------------------------------------
|C:65.1%[S:64.7%,D:0.4%],F:16.9%,M:18.0%,n:255    |
|166    Complete BUSCOs (C)                       |
|165    Complete and single-copy BUSCOs (S)       |
|1      Complete and duplicated BUSCOs (D)        |
|43     Fragmented BUSCOs (F)                     |
|46     Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------
#run with --long and --augustus
--------------------------------------------------
|Results from dataset nematoda_odb10              |
--------------------------------------------------
|C:59.6%[S:58.5%,D:1.1%],F:2.2%,M:38.2%,n:3131    |
|1866   Complete BUSCOs (C)                       |
|1831   Complete and single-copy BUSCOs (S)       |
|35     Complete and duplicated BUSCOs (D)        |
|68     Fragmented BUSCOs (F)                     |
|1197   Missing BUSCOs (M)                        |
|3131   Total BUSCO groups searched               |
--------------------------------------------------

#How many are actually missing from nematoda odb10 using both protein and genome datasets

cat MhapBuscogenoLong/run_nematoda_odb10/missing_busco_list.tsv 01_protein/*/run_nematoda_odb10/missing_busco_list.tsv |sort|uniq -c |awk '$1==2' |wc
  1050    2122   20064




```
### Meloidogyne hapla protein
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/08_Melhap
ln -s ../../../../06_Orthology/01_genomicResources/meloidogyne_hapla.PRJNA29083.WBPS15.protein.fa
echo "ml miniconda3; source activate busco5_env ; busco -i meloidogyne_hapla.PRJNA29083.WBPS15.protein.fa -o MhapBuscoprot -m prot --auto-lineage-euk -c 35 -f " >busco.sh

--------------------------------------------------
 |Results from generic domain eukaryota_odb10      |
 --------------------------------------------------
 |C:79.3%[S:76.9%,D:2.4%],F:10.2%,M:10.5%,n:255    |
 |202    Complete BUSCOs (C)                       |
 |196    Complete and single-copy BUSCOs (S)       |
 |6      Complete and duplicated BUSCOs (D)        |
 |26     Fragmented BUSCOs (F)                     |
 |27     Missing BUSCOs (M)                        |
 |255    Total BUSCO groups searched               |
 --------------------------------------------------

 --------------------------------------------------
 |Results from dataset metazoa_odb10               |
 --------------------------------------------------
 |C:62.8%[S:60.7%,D:2.1%],F:7.0%,M:30.2%,n:954     |
 |599    Complete BUSCOs (C)                       |
 |579    Complete and single-copy BUSCOs (S)       |
 |20     Complete and duplicated BUSCOs (D)        |
 |67     Fragmented BUSCOs (F)                     |
 |288    Missing BUSCOs (M)                        |
 |954    Total BUSCO groups searched               |
 --------------------------------------------------
 |Results from dataset nematoda_odb10              |
 --------------------------------------------------
 |C:59.6%[S:58.5%,D:1.1%],F:2.2%,M:38.2%,n:3131    |
 |1866   Complete BUSCOs (C)                       |
 |1831   Complete and single-copy BUSCOs (S)       |
 |35     Complete and duplicated BUSCOs (D)        |
 |68     Fragmented BUSCOs (F)                     |
 |1197   Missing BUSCOs (M)                        |
 |3131   Total BUSCO groups searched               |
 --------------------------------------------------

  #What percentage would be considered complete if we eliminate the 688 busco genes missing from all species from the total
  3131-666=2465
  1050-666 = 384

  #Total percent complete (2465- 384=2081)/2465 = 84.4% not missing
```
### Meloidogyne incognita genome
```
/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/07_Melinc
ln -s ../../../06_Orthology/01_genomicResources/meloidogyne_incognita.PRJEB8714.WBPS15.genomic.fa
echo "ml miniconda3; source activate busco5_env ; busco -i meloidogyne_incognita.PRJEB8714.WBPS15.genomic.fa -o MincBuscogeno -m geno --auto-lineage-euk -c 35 -f " >busco.sh
--------------------------------------------------
|Results from dataset eukaryota_odb10             |
--------------------------------------------------
|C:70.2%[S:18.0%,D:52.2%],F:12.2%,M:17.6%,n:255   |
|179    Complete BUSCOs (C)                       |
|46     Complete and single-copy BUSCOs (S)       |
|133    Complete and duplicated BUSCOs (D)        |
|31     Fragmented BUSCOs (F)                     |
|45     Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------

#run with --long and --augustus
--------------------------------------------------
|Results from dataset nematoda_odb10              |
--------------------------------------------------
|C:63.6%[S:24.0%,D:39.6%],F:1.0%,M:35.4%,n:3131   |
|1992   Complete BUSCOs (C)                       |
|752    Complete and single-copy BUSCOs (S)       |
|1240   Complete and duplicated BUSCOs (D)        |
|32     Fragmented BUSCOs (F)                     |
|1107   Missing BUSCOs (M)                        |
|3131   Total BUSCO groups searched               |
--------------------------------------------------



#How many are actually missing from nematoda odb10 using both protein and genome datasets
cat MincBuscogenoLong/run_nematoda_odb10/missing_busco_list.tsv 01_protein/*/run_nematoda_odb10/missing_busco_list.tsv |sort|uniq -c |awk '$1==2' |wc
   889    1800   17008




```
### Meloidogyne incognita protein
```
/work/gif/remkv6/Baum/04_DovetailSCNGenome/60_BuscoAllRedo/07_Melinc/01_protein
ln -s ../../../../06_Orthology/01_genomicResources/meloidogyne_incognita.PRJEB8714.WBPS15.protein.fa
echo "ml miniconda3; source activate busco5_env ; busco -i meloidogyne_incognita.PRJEB8714.WBPS15.protein.fa -o MincBuscoprot -m prot --auto-lineage-euk -c 35 -f " >busco.sh
--------------------------------------------------
|Results from dataset eukaryota_odb10             |
--------------------------------------------------
|C:89.8%[S:18.8%,D:71.0%],F:5.5%,M:4.7%,n:255     |
|229    Complete BUSCOs (C)                       |
|48     Complete and single-copy BUSCOs (S)       |
|181    Complete and duplicated BUSCOs (D)        |
|14     Fragmented BUSCOs (F)                     |
|12     Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------
--------------------------------------------------
|Results from dataset nematoda_odb10              |
--------------------------------------------------
|C:68.5%[S:13.9%,D:54.6%],F:1.7%,M:29.8%,n:3131   |
|2144   Complete BUSCOs (C)                       |
|435    Complete and single-copy BUSCOs (S)       |
|1709   Complete and duplicated BUSCOs (D)        |
|52     Fragmented BUSCOs (F)                     |
|935    Missing BUSCOs (M)                        |
|3131   Total BUSCO groups searched               |
--------------------------------------------------

#What percentage would be considered complete if we eliminate the 688 busco genes missing from all species from the total
3131-666=2465
889-666 = 223

#Total percent complete (2465- 223=2242)/2465 = 91.0% not missing
```

## Busco genes missing from whole clade.

```
#How many are missing from all species in the genome scans
cat */*Long/run_n*/missing_busco_list.tsv |sort |uniq -c|awk '$1==8'|wc
    763    1548   14645

#How many are missing from all species in the protein scans
cat */01_protei*/*/run_n*/missing_busco_list.tsv |sort |uniq -c |awk '$1==8'|wc
   690    1402   13259



#How many were missing from all protein and all genome scans
cat */*Long/run_n*/missing_busco_list.tsv  */01_protei/*/run_n*/missing_busco_list.tsv |sort |uniq -c |less
cat */*Long/run_n*/missing_busco_list.tsv  */01_protei*/*/run_n*/missing_busco_list.tsv |sort |uniq -c |awk '$1==16'|wc
    666    1354   12808


Busco genes missing from Meloidogyne
cat 0*_Mel*/*Long/run_n*/missing* 0*_Mel*/01_prote*/*/run_n*/missing* |sort|uniq -c |awk '$1==4' |wc
    857    1736   16407


Busco genes missing from cyst nematodes
cat 01_*/*Long/run_n*/missing* 01_*/01*/*/run_n*/missing* 02_*/*Long/run_n*/missing* 02_*/01*/*/run_n*/missing* 03_*/*Long/run_n*/missing* 03_*/01*/*/run_n*/missing* 04_*/*Long/run_n*/missing* 04_*/01*/*/run_n*/missing* 05_*/*Long/run_n*/missing* 05_*/01*/*/run_n*/missing* 06_*/*Long/run_n*/missing* 06_*/01*/*/run_n*/missing* |sort|uniq -c |awk '$1=="12"' |wc
    724    1470   13912

Busco genes missing from globodera
cat 04_*/*Long/run_n*/missing* 04_*/01*/*/run_n*/missing* 05_*/*Long/run_n*/missing* 05_*/01*/*/run_n*/missing* 06_*/*Long/run_n*/missing* 06_*/01*/*/run_n*/missing* |sort|uniq -c |awk '$1=="6"' |wc
    775    1572   14882


busco genes missing from all three SCN assemblies
cat 01_*/*Long/run_n*/missing* 01_*/01*/*/run_n*/missing* 02_*/*Long/run_n*/missing* 02_*/01*/*/run_n*/missing* 03_*/*Long/run_n*/missing* 03_*/01*/*/run_n*/missing* |sort|uniq -c |awk '$1=="6"' |wc
    813    1648   15591

```


### genome stats and gaps

```
Used new_Assemblathon.pl on these genomes.  

script to count gaps (any string of N's).  this total minus 1 is the true gap count.
tr "\n" " " <meloidogyne_incognita.PRJEB8714.WBPS15.genomic.fa |sed 's/ //g' |tr "N" "\n" |uniq |wc
```
