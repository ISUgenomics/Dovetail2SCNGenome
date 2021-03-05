#  How did dovetail do our scaffolding?
```
/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/10_MummerSelf
module load GIF/mummer/4.0.0beta2
mummer -mum -l 10000 -b -n -qthreads 15  genome738sl.polished.mitoFixed.fa  nematode_sp._22Aug2017_DZkUC.fasta >738vs104.mummer.out
```

#### transform mummer output
```
less 738vs104.mummer.out |awk -v source=0 '{if(substr($1,1,1)==">") {source=$0;} else {print source,$0}}' |awk '{if($3=="Reverse") {print $2,$6,$6+$7"\t"$4,$5,$5+$7} else {print $2,$5,$5+$6"\t"$3,$4,$4+$6}}' |sort -k1,1V -k2,3n >SyntenicRibbons.txt
```

#### make the karyotype
```
 bioawk -c fastx '{print $name,length($seq)}'  |awk '{print "chr","-",$1,$1,"0",$2,"blue"}' >dovetail.karyotype.conf
 bioawk -c fastx '{print $name,length($seq)}' genome738sl.polished.mitoFixed.fa |awk '{print "chr","-",$1,$1,"0",$2,"green"}' >738.karyotype.conf

cat 738.karyotype.conf dovetail.karyotype.conf >karyotype.conf
```

#### necessary circos files
```
cp /work/GIF/remkv6/USDA/05_OpscanIadhoreCircosAbaloneTest/circos/bands.conf .
cp /work/GIF/remkv6/USDA/05_OpscanIadhoreCircosAbaloneTest/circos/circos.conf .
cp /work/GIF/remkv6/USDA/05_OpscanIadhoreCircosAbaloneTest/circos/housekeeping.conf .
cp /work/GIF/remkv6/USDA/05_OpscanIadhoreCircosAbaloneTest/circos/ideogram.conf .
cp /work/GIF/remkv6/USDA/05_OpscanIadhoreCircosAbaloneTest/circos/ticks.conf .


#made some edits to text size in ticks and ideogram files.
```


####Get rid of that annoying ;HRSCAF... stuff from Dovetail
```
less dovetail.karyotype.conf |sed 's/;/\t/'g| awk '{print $1,$2,$3,$5,$7,$8,$9,$10} '|tr " " "\t" >dovetail.revisedkaryotype.conf

cat dovetail.revisedkaryotype.conf 738.karyotype.conf |tr " " "\t" >karyotype.conf

less SyntenicRibbons.txt |sed 's/;/\t/g' |awk '{print $1,$3,$4,$5,$6,$7}' |tr " " "\t" >SyntenicRibbons.revised.txt
```

#### reorganize the bands so they can be seen
```
cp -rf /work/GIF/remkv6/USDA/05_OpscanIadhoreCircosAbaloneTest/circos/circos-tools-0.22 .


# took like 20 mins
circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.revised.txt -karyotype karyotype.conf

#had some issues with copy and paste of the above output, so had to run this to get rid of newlines
tr "\n" " " <reorderfix |sed 's/ //g' |cat circos.conf - >circosreordered.conf
circos -conf circosreordered.conf
```
![Dovetail genome alignment to 738 genome](../assets/circos_738_dovetail.png)
