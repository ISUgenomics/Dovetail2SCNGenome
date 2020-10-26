# Identify methylation status in the genome using existing pacbio reads



### Example run of smrttools ipdSummary
```
#explanation of use found here
https://www.pacb.com/wp-content/uploads/SMRT-Tools-Reference-Guide-v8.0.pdf

# this is their example command
ipdSummary aligned.bam --reference ref.fasta  m6A,m4C --gff basemods.gff --csv_h5 kinetics.h5

```


### smrttools setup, attempt to use previously aligned data *FAIL
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/55_Methylation

ln -s ../30_ReadMapping/02_Subread/SCNgenome.AllSubreads_sorted.bam
ln -s ../30_ReadMapping/02_Subread/SCNgenome.AllSubreads_sorted.bam.bai
ln -s ../30_ReadMapping/02_Subread/SCNgenome.fasta

ml dafoam; ml smrttools/4.0; ipdSummary SCNgenome.AllSubreads_sorted.bam --reference SCNgenome.fasta  --gff basemods.gff --debug


#error
#KeyError: 'pb'
#this is likely from the program wanting a pbindex of the bam file.

#pbindex SCNgenome.AllSubreads_sorted.bam
#error
#pbindex ERROR: read group ID not found
#this is likely from using a different aligner.

```



### Blasr alignment
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/55_Methylation

#removed previous softlinks
ln -s ../49_RenameChromosomes/01_Transfer2Box/SCNgenome.fasta


for f in  /work/GIF/archive1/Baum/091615_SCN/RawReads/*/Analysis_Results/subreads.bam; do ls $f;done >bamFile.fofn

echo "ml dafoam;ml blasr/5.1; blasr bamFile.fofn SCNgenome.fasta  --nproc 16 --useQuality --bam --out subreadsPBbam.bam;ml smrttools/4.0; pbindex subreadsPBbam.bam; ipdSummary subreadsPBbam.bam --reference SCNgenome.fasta  --gff basemods.gff"


```
