### Need NLS to predict effectors
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/42_NLSPrediction

wget http://www.moseslab.csb.utoronto.ca/NLStradamus/NLStradamus/NLStradamus.1.8.tar.gz
tar -zxvf NLStradamus.1.8.tar.gz
perl nlstradamus.pl -tab -cpu 16 -i OrderedSCNGenePredictionsVHEJ_proteins.fasta >NLSsignals

less NLSsignals |awk 'NR>1{print $1"\t"$3}'  |sort -k1,1 -u  >NLSsignals.tab

```
