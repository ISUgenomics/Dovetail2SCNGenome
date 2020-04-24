# I need an assessment of structural variation in the genome via long reads.  

Which ones could be candidates for endoreduplication?

### Run SVIM
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/03_Subread


ml miniconda3
conda create -n svim_env --channel bioconda svim
source activate svim_env
conda install --channel bioconda svim

echo "ml miniconda3; source activate svim_env; svim alignment /work/gif/remkv6/Baum/04_DovetailSCNGenome/03_Subread SCNgenome.AllSubreads_sorted.bam SCNgenome.fasta " >SVIM.sh

```

### upload vcf output to jbrowse
```
#/work/gif/remkv6/Baum/04_DovetailSCNGenome/03_Subread

#could not get this formatted for jbrowse...


#how many variants do I find?
(svim_env) [remkv6@nova006 03_Subread]$ less candidates/candidates_novel_insertions.bed |awk '($3-$2)>100' |wc
  26769  187383 6881248
(svim_env) [remkv6@nova006 03_Subread]$ less candidates/candidates_novel_insertions.bed |awk '($3-$2)>1000' |wc
    905    6335  200214
(svim_env) [remkv6@nova006 03_Subread]$ less candidates/candidates_novel_insertions.bed |awk '($3-$2)>5000' |wc
     20     140    3334

```
