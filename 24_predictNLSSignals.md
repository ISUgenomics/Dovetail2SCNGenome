### Need NLS to predict effectors
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/42_NLSPrediction

wget http://www.moseslab.csb.utoronto.ca/NLStradamus/NLStradamus/NLStradamus.1.8.tar.gz
tar -zxvf NLStradamus.1.8.tar.gz
 perl nlstradamus.pl -tab -cpu 16 -i OrderedSCNGenePredictionsVHEJ_proteins.fasta >NLSsignals

less NLSsignals.tab |sort -u -k1,1V >SingleBestScoreNLS.tab

```
