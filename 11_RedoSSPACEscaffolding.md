# Andrew and Sebastian say that the default settings for sspace are rediculously aggressive.  We keep getting a 31Mb scaffold, which is at the upper edge of chromosome sizes possible for this genome (n=9 with approximately equal size).  

### Using CCSreads 100bp min
```
module load sspace-longread/1.1-5ae3ijg
perl SSPACE-LongRead.pl -c ../nematode_sp_19Jul2018_IbtP1.fasta -p ../reads_of_insert.fastq -t 16 -s ../ccs2genome.out -o 100 -k 1

 new_Assemblathon.pl scaffolds100.fasta

---------------- Information for assembly 'scaffolds100.fasta' ----------------


                                         Number of scaffolds        527
                                     Total size of scaffolds  166493275
                                            Longest scaffold   31384897
                                           Shortest scaffold        536
                                 Number of scaffolds > 1K nt        526  99.8%
                                Number of scaffolds > 10K nt        518  98.3%
                               Number of scaffolds > 100K nt        135  25.6%
                                 Number of scaffolds > 1M nt         15   2.8%
                                Number of scaffolds > 10M nt          5   0.9%
                                          Mean scaffold size     315927
                                        Median scaffold size      49377
                                         N50 scaffold length   11380747
                                          L50 scaffold count          5
                                                 scaffold %A      31.57
                                                 scaffold %C      18.38
                                                 scaffold %G      18.39
                                                 scaffold %T      31.54
                                                 scaffold %N       0.12
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      92.4%
              Percentage of assembly in unscaffolded contigs       7.6%
                      Average number of contigs per scaffold        3.6
Average length of break (>25 Ns) between contigs in scaffold        141

                                           Number of contigs       1917
                              Number of contigs in scaffolds       1601
                          Number of contigs not in scaffolds        316
                                       Total size of contigs  166297037
                                              Longest contig    1155511
                                             Shortest contig        173
                                   Number of contigs > 1K nt       1904  99.3%
                                  Number of contigs > 10K nt       1863  97.2%
                                 Number of contigs > 100K nt        459  23.9%
                                   Number of contigs > 1M nt          3   0.2%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size      86749
                                          Median contig size      51246
                                           N50 contig length     143023
                                            L50 contig count        298
                                                   contig %A      31.61
                                                   contig %C      18.40
                                                   contig %G      18.41
                                                   contig %T      31.58
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0

```
### CCS reads at 2k

```
module load sspace-longread/1.1-5ae3ijg
perl SSPACE-LongRead.pl -c ../nematode_sp_19Jul2018_IbtP1.fasta -p ../reads_of_insert.fastq -t 16 -s ../ccs2genome.out -o 2000 -k 1


Number of scaffolds        527
Total size of scaffolds  166493275
   Longest scaffold   31384897
  Shortest scaffold        536
Number of scaffolds > 1K nt        526  99.8%
Number of scaffolds > 10K nt        518  98.3%
Number of scaffolds > 100K nt        135  25.6%
Number of scaffolds > 1M nt         15   2.8%
Number of scaffolds > 10M nt          5   0.9%
 Mean scaffold size     315927
Median scaffold size      49377
N50 scaffold length   11380747
 L50 scaffold count          5
        scaffold %A      31.57
        scaffold %C      18.38
        scaffold %G      18.39
        scaffold %T      31.54
        scaffold %N       0.12
scaffold %non-ACGTN       0.00
Number of scaffold non-ACGTN nt          0

Percentage of assembly in scaffolded contigs      92.4%
Percentage of assembly in unscaffolded contigs       7.6%
Average number of contigs per scaffold        3.6
Average length of break (>25 Ns) between contigs in scaffold        141

  Number of contigs       1917
Number of contigs in scaffolds       1601
Number of contigs not in scaffolds        316
Total size of contigs  166297037
     Longest contig    1155511
    Shortest contig        173
Number of contigs > 1K nt       1904  99.3%
Number of contigs > 10K nt       1863  97.2%
Number of contigs > 100K nt        459  23.9%
Number of contigs > 1M nt          3   0.2%
Number of contigs > 10M nt          0   0.0%
   Mean contig size      86749
 Median contig size      51246
  N50 contig length     143023
   L50 contig count        298
          contig %A      31.61
          contig %C      18.40
          contig %G      18.41
          contig %T      31.58
          contig %N       0.00
  contig %non-ACGTN       0.00
Number of contig non-ACGTN nt          0

```

### Because subreads seemed to have a better effect on collapsing scaffolds, I used those with varying parameters that were much more stringent than default below. 
```
module load sspace-longread/1.1-5ae3ijg
perl SSPACE-LongRead.pl -c nematode_sp_19Jul2018_IbtP1.fasta -p SCN.all.subreads.sl.7k.fasta -t 12 -o 1000 -k 1 -l 10 -b 1000_10Scaffolding -s PacBio_scaffolder_results/intermediate_files/BLASR_results.txt
module load sspace-longread/1.1-5ae3ijg
perl SSPACE-LongRead.pl -c nematode_sp_19Jul2018_IbtP1.fasta -p SCN.all.subreads.sl.7k.fasta -t 12 -o 1000 -k 1 -l 5 -b 1000_5Scaffolding -s PacBio_scaffolder_results/intermediate_files/BLASR_results.txt
module load sspace-longread/1.1-5ae3ijg
perl SSPACE-LongRead.pl -c nematode_sp_19Jul2018_IbtP1.fasta -p SCN.all.subreads.sl.7k.fasta -t 12 -o 1000 -k 1 -l 15 -b 1000_15Scaffolding -s PacBio_scaffolder_results/intermediate_files/BLASR_results.txt
module load sspace-longread/1.1-5ae3ijg
perl SSPACE-LongRead.pl -c nematode_sp_19Jul2018_IbtP1.fasta -p SCN.all.subreads.sl.7k.fasta -t 12 -o 1000 -k 1 -l 20 -b 1000_20Scaffolding -s PacBio_scaffolder_results/intermediate_files/BLASR_results.txt

```
