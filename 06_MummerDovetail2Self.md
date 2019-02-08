# How many duplications are present in the dovetail 2 genomes


### run mummer
```
sed  's/>/>2/g' nematode_sp_19Jul2018_IbtP1.fasta >Genome2.fasta
ln -s ../16_NewDovetailGenome/nematode_sp_19Jul2018_IbtP1.fasta

module load GIF/mummer/4.0.0beta2
mummer -mum -l 1000 -b -n -qthreads 12  nematode_sp_19Jul2018_IbtP1.fasta Genome2.fasta >2vs2.mummer.out
```

### How large are the duplications
```
less 2vs2.mummer.out |sed 's/2Scaffold/Scaffold/g' |awk -v scaff=1 '{if($1==">") {print $0,scaff=$2} else if($1==scaff  && $2!=$3) {print $0}}' |awk '$1!=">" {print $4}' |summary.sh
Total:  3,430,752
Count:  1,909
Mean:   1,797
Median: 1,373
Min:    1,000
Max:    19,730


#here is an actual representation of the duplications in the genome, in a bed format.  However, it will need to be modified to show which duplicated fragment belongs to which.  
less 2vs2.mummer.out |sed 's/2Scaffold/Scaffold/g' |awk -v scaff=1 '{if($1==">") {print $0,scaff=$2} else if($1==scaff  && $2!=$3) {print $0}}' |awk '$1!=">" ' |awk -v scaff=1 '{print $1,$2,$2+$4"\n"$1,$3,$3+$4}' |sort -k1,1V -k2,3n | tr " " "\t" |bedtools merge -d 500 >SelfMummer.bed

less SelfMummer.bed |awk '{print $3-$2}' |summary.sh
Total:  6,560,828
Count:  2,761
Mean:   2,376
Median: 1,565
Min:    1,000
Max:    26,508

```
### Make circos plot

```
#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/22_Dovetail2MummerSelf/01_Circos
less 2vs2.mummer.out |awk -v source=0 '{if(substr($1,1,1)==">") {source=$0;} else {print source,$0}}' |awk '{if($3=="Reverse") {print $2,$6,$6+$7"\t"$4,$5,$5+$7} else {print $2,$5,$5+$6"\t"$3,$4,$4+$6}}' |awk '$2!=$5 && $3!=$6' |sort -k1,1V -k2,3n >SyntenicRibbons.txt

#/work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/22_Dovetail2MummerSelf/01_Circos
cp ../../18_2692MummerDovetail2/bands.conf .
cp ../../18_2692MummerDovetail2/circos.conf .
cp ../../18_2692MummerDovetail2/housekeeping.conf .
cp ../../18_2692MummerDovetail2/ideogram.conf .
cp ../../18_2692MummerDovetail2/ticks.conf .
cp ../../18_2692MummerDovetail2/karyotype.conf .

less ../SyntenicRibbons.txt |sed 's/;/\t/g' |awk '{print $1,$3,$4,$5,$7,$8}' |sed 's/2//1' |awk '($3-$2)>2000' >SyntenicRibbons.conf


circos -conf circos.conf

module use /work/GIF/software/modules
 module load GIF2/circos
[remkv6@condo163 01_Circos]$ /shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons.conf -karyotype karyotype.conf
calculating round 0
report round 0 minimize init 271029 final 165323 change 39.00%
calculating round 1
report round 1 minimize init 165323 final 109707 change 33.64%
calculating round 2
report round 2 minimize init 109707 final 81236 change 25.95%
calculating round 3
report round 3 minimize init 81236 final 73117 change 9.99%
scorereport init 271029 final 73117 change 73.02%
chromosomes_order = Scaffold_390,Scaffold_9,Scaffold_78,Scaffold_146,Scaffold_32,Scaffold_191,Scaffold_197,Scaffold_327,Scaffold_512,Scaffold_252,Scaffold_221,Scaffold_2,Scaffold_36,Scaffold_317,Scaffold_396,Scaffold_555,Scaffold_213,Scaffold_423,Scaffold_531,Scaffold_324,Scaffold_292,Scaffold_192,Scaffold_77,Scaffold_103,Scaffold_573,Scaffold_472,Scaffold_498,Scaffold_351,Scaffold_287,Scaffold_101,Scaffold_273,Scaffold_127,Scaffold_392,Scaffold_227,Scaffold_459,Scaffold_38,Scaffold_288,Scaffold_89,Scaffold_538,Scaffold_244,Scaffold_34,Scaffold_163,Scaffold_308,Scaffold_228,Scaffold_225,Scaffold_311,Scaffold_401,Scaffold_277,Scaffold_42,Scaffold_364,Scaffold_150,Scaffold_69,Scaffold_126,Scaffold_485,Scaffold_565,Scaffold_373,Scaffold_458,Scaffold_91,Scaffold_281,Scaffold_181,Scaffold_159,Scaffold_53,Scaffold_366,Scaffold_264,Scaffold_445,Scaffold_517,Scaffold_399,Scaffold_475,Scaffold_402,Scaffold_148,Scaffold_562,Scaffold_297,Scaffold_79,Scaffold_489,Scaffold_30,Scaffold_315,Scaffold_263,Scaffold_558,Scaffold_570,Scaffold_138,Scaffold_417,Scaffold_420,Scaffold_41,Scaffold_586,Scaffold_343,Scaffold_363,Scaffold_80,Scaffold_17,Scaffold_307,Scaffold_466,Scaffold_189,Scaffold_104,Scaffold_436,Scaffold_76,Scaffold_545,Scaffold_559,Scaffold_474,Scaffold_550,Scaffold_345,Scaffold_425,Scaffold_539,Scaffold_575,Scaffold_5,Scaffold_280,Scaffold_175,Scaffold_551,Scaffold_58,Scaffold_59,Scaffold_172,Scaffold_142,Scaffold_231,Scaffold_572,Scaffold_60,Scaffold_291,Scaffold_220,Scaffold_73,Scaffold_492,Scaffold_560,Scaffold_359,Scaffold_218,Scaffold_182,Scaffold_406,Scaffold_43,Scaffold_232,Scaffold_584,Scaffold_235,Scaffold_415,Scaffold_140,Scaffold_143,Scaffold_461,Scaffold_428,Scaffold_553,Scaffold_430,Scaffold_117,Scaffold_125,Scaffold_18,Scaffold_319,Scaffold_411,Scaffold_380,Scaffold_484,Scaffold_116,Scaffold_211,Scaffold_87,Scaffold_98,Scaffold_409,Scaffold_453,Scaffold_8,Scaffold_188,Scaffold_282,Scaffold_441,Scaffold_266,Scaffold_314,Scaffold_522,Scaffold_81,Scaffold_230,Scaffold_180,Scaffold_452,Scaffold_165,Scaffold_518,Scaffold_121,Scaffold_289,Scaffold_284,Scaffold_147,Scaffold_119,Scaffold_576,Scaffold_457,Scaffold_75,Scaffold_14,Scaffold_238,Scaffold_313,Scaffold_295,Scaffold_66,Scaffold_22,Scaffold_395,Scaffold_63,Scaffold_243,Scaffold_488,Scaffold_525,Scaffold_540,Scaffold_413,Scaffold_6,Scaffold_160,Scaffold_265,Scaffold_328,Scaffold_110,Scaffold_173,Scaffold_355,Scaffold_193,Scaffold_152,Scaffold_347,Scaffold_23,Scaffold_257,Scaffold_27,Scaffold_569,Scaffold_94,Scaffold_190,Scaffold_229,Scaffold_122,Scaffold_310,Scaffold_135,Scaffold_205,Scaffold_233,Scaffold_469,Scaffold_259,Scaffold_320,Scaffold_400,Scaffold_468,Scaffold_40,Scaffold_29,Scaffold_532,Scaffold_133,Scaffold_495,Scaffold_217,Scaffold_20,Scaffold_419,Scaffold_331,Scaffold_248,Scaffold_71,Scaffold_407,Scaffold_464,Scaffold_537,Scaffold_303,Scaffold_440,Scaffold_45,Scaffold_109,Scaffold_203,Scaffold_298,Scaffold_337,Scaffold_151,Scaffold_255,Scaffold_497,Scaffold_356,Scaffold_354,Scaffold_513,Scaffold_276,Scaffold_202,Scaffold_206,Scaffold_196,Scaffold_19,Scaffold_368,Scaffold_486,Scaffold_223,Scaffold_139,Scaffold_201,Scaffold_240,Scaffold_102,Scaffold_432,Scaffold_465,Scaffold_21,Scaffold_323,Scaffold_155,Scaffold_156,Scaffold_306,Scaffold_429,Scaffold_68,Scaffold_526,Scaffold_504,Scaffold_31,Scaffold_132,Scaffold_414,Scaffold_241,Scaffold_62,Scaffold_283,Scaffold_462,Scaffold_269,Scaffold_405,Scaffold_500,Scaffold_179,Scaffold_245,Scaffold_481,Scaffold_483,Scaffold_194,Scaffold_64,Scaffold_24,Scaffold_408,Scaffold_55,Scaffold_370,Scaffold_339,Scaffold_254,Scaffold_131,Scaffold_299,Scaffold_506,Scaffold_360,Scaffold_250,Scaffold_533,Scaffold_114,Scaffold_448,Scaffold_177,Scaffold_85,Scaffold_216,Scaffold_523,Scaffold_566,Scaffold_162,Scaffold_198,Scaffold_13,Scaffold_352,Scaffold_322,Scaffold_530,Scaffold_377,Scaffold_384,Scaffold_214,Scaffold_556,Scaffold_153,Scaffold_136,Scaffold_326,Scaffold_482,Scaffold_178,Scaffold_123,Scaffold_74,Scaffold_528,Scaffold_321,Scaffold_285,Scaffold_25,Scaffold_579,Scaffold_99,Scaffold_195,Scaffold_157,Scaffold_115,Scaffold_82,Scaffold_209,Scaffold_404,Scaffold_581,Scaffold_10,Scaffold_426,Scaffold_61,Scaffold_92,Scaffold_333,Scaffold_96,Scaffold_342,Scaffold_210,Scaffold_318,Scaffold_90,Scaffold_11,Scaffold_305,Scaffold_534,Scaffold_433,Scaffold_48,Scaffold_100,Scaffold_212,Scaffold_106,Scaffold_329,Scaffold_112,Scaffold_128,Scaffold_174,Scaffold_473,Scaffold_334,Scaffold_491,Scaffold_394,Scaffold_451,Scaffold_300,Scaffold_385,Scaffold_50,Scaffold_403,Scaffold_237,Scaffold_582,Scaffold_567,Scaffold_507,Scaffold_358,Scaffold_535,Scaffold_346,Scaffold_28,Scaffold_144,Scaffold_583,Scaffold_249,Scaffold_577,Scaffold_164


### ONLY LARGER THAN 3K
less SyntenicRibbons.conf |awk '($3-$2)>3000' >SyntenicRibbons3k.conf
awk '{print $1"\n"$4}' SyntenicRibbons3k.conf |sort|uniq|grep -w -f - karyotype.conf >karyotype3k

/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links SyntenicRibbons3k.conf -karyotype karyotype3k
calculating round 0
report round 0 minimize init 54813 final 25405 change 53.65%
calculating round 1
report round 1 minimize init 25405 final 16398 change 35.45%
calculating round 2
report round 2 minimize init 16398 final 12749 change 22.25%
calculating round 3
report round 3 minimize init 12749 final 11482 change 9.94%
scorereport init 54813 final 11482 change 79.05%
chromosomes_order = Scaffold_573,Scaffold_351,Scaffold_248,Scaffold_131,Scaffold_13,Scaffold_324,Scaffold_74,Scaffold_115,Scaffold_530,Scaffold_213,Scaffold_555,Scaffold_518,Scaffold_36,Scaffold_2,Scaffold_244,Scaffold_103,Scaffold_497,Scaffold_390,Scaffold_160,Scaffold_29,Scaffold_32,Scaffold_221,Scaffold_127,Scaffold_401,Scaffold_30,Scaffold_288,Scaffold_101,Scaffold_273,Scaffold_163,Scaffold_206,Scaffold_474,Scaffold_5,Scaffold_281,Scaffold_80,Scaffold_523,Scaffold_53,Scaffold_159,Scaffold_535,Scaffold_282,Scaffold_373,Scaffold_311,Scaffold_567,Scaffold_138,Scaffold_307,Scaffold_181,Scaffold_87,Scaffold_41,Scaffold_18,Scaffold_42,Scaffold_22,Scaffold_581,Scaffold_43,Scaffold_411,Scaffold_201,Scaffold_40,Scaffold_485,Scaffold_73,Scaffold_558,Scaffold_359,Scaffold_560,Scaffold_570,Scaffold_551,Scaffold_263,Scaffold_399,Scaffold_295,Scaffold_220,Scaffold_572,Scaffold_291,Scaffold_172,Scaffold_175,Scaffold_117,Scaffold_266,Scaffold_193,Scaffold_310,Scaffold_522,Scaffold_143,Scaffold_68,Scaffold_64,Scaffold_98,Scaffold_104,Scaffold_575,Scaffold_417,Scaffold_333,Scaffold_243,Scaffold_210,Scaffold_59,Scaffold_92,Scaffold_305,Scaffold_385,Scaffold_334,Scaffold_346,Scaffold_152,Scaffold_157,Scaffold_173,Scaffold_289,Scaffold_355,Scaffold_583,Scaffold_11,Scaffold_306,Scaffold_237,Scaffold_195,Scaffold_66,Scaffold_285,Scaffold_69,Scaffold_419,Scaffold_142,Scaffold_188,Scaffold_298,Scaffold_453,Scaffold_337,Scaffold_407,Scaffold_303,Scaffold_55,Scaffold_402,Scaffold_102,Scaffold_428,Scaffold_405,Scaffold_151,Scaffold_354,Scaffold_537,Scaffold_356,Scaffold_202,Scaffold_525,Scaffold_488,Scaffold_465,Scaffold_265,Scaffold_155,Scaffold_45,Scaffold_223,Scaffold_352,Scaffold_19,Scaffold_534,Scaffold_190,Scaffold_569,Scaffold_94,Scaffold_484,Scaffold_368,Scaffold_292,Scaffold_507,Scaffold_245,Scaffold_486,Scaffold_432,Scaffold_429,Scaffold_31,Scaffold_177,Scaffold_240,Scaffold_23,Scaffold_255,Scaffold_139,Scaffold_89,Scaffold_313,Scaffold_132,Scaffold_109,Scaffold_495,Scaffold_156,Scaffold_214,Scaffold_203,Scaffold_441,Scaffold_205,Scaffold_238,Scaffold_235,Scaffold_320,Scaffold_532,Scaffold_133,Scaffold_440,Scaffold_257,Scaffold_217,Scaffold_259,Scaffold_526,Scaffold_377,Scaffold_483,Scaffold_370,Scaffold_481,Scaffold_254,Scaffold_556,Scaffold_225,Scaffold_339,Scaffold_308,Scaffold_153,Scaffold_63,Scaffold_198,Scaffold_284,Scaffold_322,Scaffold_136,Scaffold_360,Scaffold_506,Scaffold_299,Scaffold_404,Scaffold_24,Scaffold_380,Scaffold_75,Scaffold_100,Scaffold_528,Scaffold_48,Scaffold_436,Scaffold_586,Scaffold_413,Scaffold_189,Scaffold_363,Scaffold_21,Scaffold_323,Scaffold_538,Scaffold_231,Scaffold_280,Scaffold_425,Scaffold_319,Scaffold_14,Scaffold_531,Scaffold_513,Scaffold_457,Scaffold_6,Scaffold_384,Scaffold_283,Scaffold_473,Scaffold_452,Scaffold_252,Scaffold_462,Scaffold_406,Scaffold_269,Scaffold_533,Scaffold_423,Scaffold_99,Scaffold_230,Scaffold_584,Scaffold_38,Scaffold_121,Scaffold_165,Scaffold_328,Scaffold_227,Scaffold_126,Scaffold_539,Scaffold_112,Scaffold_174,Scaffold_178,Scaffold_491,Scaffold_96,Scaffold_216,Scaffold_144,Scaffold_81,Scaffold_78,Scaffold_392,Scaffold_211,Scaffold_128,Scaffold_192,Scaffold_77,Scaffold_326,Scaffold_582,Scaffold_61,Scaffold_342,Scaffold_415,Scaffold_553,Scaffold_451,Scaffold_482


```

I was trying to perform a bedtools merge on the overlapping gffs that are separated by 5kb.  not successful.  
### separate out largest scaffolds
```
#make manual list of largest scaffolds with mummer dups
vi grep.list # counted from circos plot
less grep.list |while read line; do grep \"$line\" SyntenicRibbons.conf;done |tr "\t" "\n" |>\""$line"SyntenicRibbons.conf";done >RibbonSubsetter.sh
sh RibbonSubsetter

#individual karyotypes also
for f in S*SyntenicRibbons.conf; do awk '{print $1"\n"$4}' $f|sort|uniq|grep -w -f - karyotype.conf  >${f%.*}Karyotype.conf;done


#Now to iterate through the circos.conf file


##18
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links Scaffold_18SyntenicRibbons.conf -karyotype Scaffold_18SyntenicRibbonsKaryotype.conf
chromosomes_order = Scaffold_116,Scaffold_282,Scaffold_180,Scaffold_195,Scaffold_287,Scaffold_182,Scaffold_284,Scaffold_18,Scaffold_8,Scaffold_125,Scaffold_538,Scaffold_363,Scaffold_584,Scaffold_400,Scaffold_408,Scaffold_188,Scaffold_453,Scaffold_189,Scaffold_572,Scaffold_411,Scaffold_310,Scaffold_181,Scaffold_159,Scaffold_366,Scaffold_87

##32
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links Scaffold_32SyntenicRibbons.conf -karyotype Scaffold_32SyntenicRibbonsKaryotype.conf

chromosomes_order = Scaffold_518,Scaffold_32,Scaffold_2,Scaffold_419,Scaffold_575,Scaffold_572,Scaffold_101,Scaffold_573,Scaffold_324,Scaffold_326,Scaffold_178,Scaffold_482,Scaffold_581,Scaffold_136,Scaffold_322,Scaffold_243,Scaffold_11,Scaffold_321,Scaffold_532,Scaffold_257,Scaffold_21,Scaffold_323,Scaffold_121,Scaffold_328,Scaffold_109,Scaffold_320,Scaffold_198,Scaffold_384,Scaffold_327,Scaffold_127,Scaffold_329,Scaffold_172


##74
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links Scaffold_74SyntenicRibbons.conf -karyotype Scaffold_74SyntenicRibbonsKaryotype.conf

chromosomes_order = Scaffold_74,Scaffold_62,Scaffold_13,Scaffold_426,Scaffold_584,Scaffold_530,Scaffold_115,Scaffold_573,Scaffold_518,Scaffold_248

##101
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links Scaffold_101SyntenicRibbons.conf -karyotype Scaffold_101SyntenicRibbonsKaryotype.conf
chromosomes_order = Scaffold_288,Scaffold_459,Scaffold_532,Scaffold_191,Scaffold_206,Scaffold_573,Scaffold_20,Scaffold_474,Scaffold_142,Scaffold_401,Scaffold_575,Scaffold_468,Scaffold_528,Scaffold_538,Scaffold_101,Scaffold_583,Scaffold_273,Scaffold_38,Scaffold_436,Scaffold_572,Scaffold_32,Scaffold_172,Scaffold_558,Scaffold_127,Scaffold_498,Scaffold_430,Scaffold_280,Scaffold_58

##127
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links Scaffold_127SyntenicRibbons.conf -karyotype Scaffold_127SyntenicRibbonsKaryotype.conf
chromosomes_order = Scaffold_159,Scaffold_308,Scaffold_287,Scaffold_206,Scaffold_281,Scaffold_489,Scaffold_472,Scaffold_420,Scaffold_30,Scaffold_227,Scaffold_53,Scaffold_327,Scaffold_127,Scaffold_32,Scaffold_5,Scaffold_584,Scaffold_461,Scaffold_172,Scaffold_163,Scaffold_321,Scaffold_573,Scaffold_101,Scaffold_565,Scaffold_148,Scaffold_572,Scaffold_538
##172
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links Scaffold_172SyntenicRibbons.conf -karyotype Scaffold_172SyntenicRibbonsKaryotype.conf
chromosomes_order = Scaffold_127,Scaffold_138,Scaffold_485,Scaffold_575,Scaffold_5,Scaffold_175,Scaffold_243,Scaffold_140,Scaffold_101,Scaffold_486,Scaffold_573,Scaffold_329,Scaffold_522,Scaffold_172,Scaffold_117,Scaffold_159,Scaffold_583,Scaffold_104,Scaffold_551,Scaffold_218

##413
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links Scaffold_413SyntenicRibbons.conf -karyotype Scaffold_413SyntenicRibbonsKaryotype.conf
chromosomes_order = Scaffold_193,Scaffold_413,Scaffold_6,Scaffold_310,Scaffold_34,Scaffold_259

##474
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links Scaffold_474SyntenicRibbons.conf -karyotype Scaffold_474SyntenicRibbonsKaryotype.conf
chromosomes_order = Scaffold_474,Scaffold_99,Scaffold_197,Scaffold_206,Scaffold_159,Scaffold_45,Scaffold_17,Scaffold_575,Scaffold_572,Scaffold_210,Scaffold_584,Scaffold_583,Scaffold_343,Scaffold_545,Scaffold_5,Scaffold_573,Scaffold_101,Scaffold_238


##572
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links Scaffold_572SyntenicRibbons.conf -karyotype Scaffold_572SyntenicRibbonsKaryotype.conf
chromosomes_order = Scaffold_572,Scaffold_474,Scaffold_584,Scaffold_406,Scaffold_558,Scaffold_390,Scaffold_518,Scaffold_18,Scaffold_573,Scaffold_159,Scaffold_560,Scaffold_441,Scaffold_127,Scaffold_32,Scaffold_96,Scaffold_295,Scaffold_5,Scaffold_73,Scaffold_291,Scaffold_11,Scaffold_359,Scaffold_101,Scaffold_430,Scaffold_220,Scaffold_59,Scaffold_583

##573
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links Scaffold_573SyntenicRibbons.conf -karyotype Scaffold_573SyntenicRibbonsKaryotype.conf
chromosomes_order = Scaffold_573,Scaffold_14,Scaffold_559,Scaffold_572,Scaffold_235,Scaffold_75,Scaffold_317,Scaffold_11,Scaffold_581,Scaffold_103,Scaffold_534,Scaffold_390,Scaffold_216,Scaffold_518,Scaffold_127,Scaffold_77,Scaffold_214,Scaffold_555,Scaffold_126,Scaffold_192,Scaffold_9,Scaffold_281,Scaffold_351,Scaffold_82,Scaffold_485,Scaffold_38,Scaffold_78,Scaffold_228,Scaffold_531,Scaffold_74,Scaffold_252,Scaffold_287,Scaffold_101,Scaffold_461,Scaffold_221,Scaffold_484,Scaffold_428,Scaffold_324,Scaffold_172,Scaffold_577,Scaffold_423,Scaffold_315,Scaffold_249,Scaffold_284,Scaffold_244,Scaffold_230,Scaffold_583,Scaffold_582,Scaffold_292,Scaffold_159,Scaffold_584,Scaffold_36,Scaffold_474,Scaffold_466,Scaffold_32,Scaffold_329,Scaffold_346,Scaffold_66,Scaffold_433,Scaffold_27

##575
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links Scaffold_575SyntenicRibbons.conf -karyotype Scaffold_575SyntenicRibbonsKaryotype.conf
chromosomes_order = Scaffold_575,Scaffold_142,Scaffold_159,Scaffold_266,Scaffold_417,Scaffold_5,Scaffold_308,Scaffold_475,Scaffold_474,Scaffold_59,Scaffold_101,Scaffold_121,Scaffold_98,Scaffold_550,Scaffold_104,Scaffold_32,Scaffold_172
##583
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links Scaffold_583SyntenicRibbons.conf -karyotype Scaffold_583SyntenicRibbonsKaryotype.conf
chromosomes_order = Scaffold_583,Scaffold_452,Scaffold_573,Scaffold_305,Scaffold_572,Scaffold_318,Scaffold_273,Scaffold_101,Scaffold_334,Scaffold_474,Scaffold_385,Scaffold_237,Scaffold_584,Scaffold_172,Scaffold_346,Scaffold_28

##584
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links Scaffold_584SyntenicRibbons.conf -karyotype Scaffold_584SyntenicRibbonsKaryotype.conf

```
### Circos config Files En masse
```

################################################################################
karyotype =Scaffold_18SyntenicRibbonsKaryotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=Scaffold_18SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
chromosomes_order = Scaffold_116,Scaffold_282,Scaffold_180,Scaffold_195,Scaffold_287,Scaffold_182,Scaffold_284,Scaffold_18,Scaffold_8,Scaffold_125,Scaffold_538,Scaffold_363,Scaffold_584,Scaffold_400,Scaffold_408,Scaffold_188,Scaffold_453,Scaffold_189,Scaffold_572,Scaffold_411,Scaffold_310,Scaffold_181,Scaffold_159,Scaffold_366,Scaffold_87
<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 ################################################################################


 ################################################################################
 karyotype =Scaffold_32SyntenicRibbonsKaryotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=Scaffold_32SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
chromosomes_order = Scaffold_518,Scaffold_32,Scaffold_2,Scaffold_419,Scaffold_575,Scaffold_572,Scaffold_101,Scaffold_573,Scaffold_324,Scaffold_326,Scaffold_178,Scaffold_482,Scaffold_581,Scaffold_136,Scaffold_322,Scaffold_243,Scaffold_11,Scaffold_321,Scaffold_532,Scaffold_257,Scaffold_21,Scaffold_323,Scaffold_121,Scaffold_328,Scaffold_109,Scaffold_320,Scaffold_198,Scaffold_384,Scaffold_327,Scaffold_127,Scaffold_329,Scaffold_172

<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 ################################################################################


 ################################################################################
 karyotype =Scaffold_74SyntenicRibbonsKaryotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=Scaffold_74SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
chromosomes_order = Scaffold_74,Scaffold_62,Scaffold_13,Scaffold_426,Scaffold_584,Scaffold_530,Scaffold_115,Scaffold_573,Scaffold_518,Scaffold_248


<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 ################################################################################


 ################################################################################
 karyotype =Scaffold_101SyntenicRibbonsKaryotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=Scaffold_101SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
Chromosomes_order = Scaffold_288,Scaffold_459,Scaffold_532,Scaffold_191,Scaffold_206,Scaffold_573,Scaffold_20,Scaffold_474,Scaffold_142,Scaffold_401,Scaffold_575,Scaffold_468,Scaffold_528,Scaffold_538,Scaffold_101,Scaffold_583,Scaffold_273,Scaffold_38,Scaffold_436,Scaffold_572,Scaffold_32,Scaffold_172,Scaffold_558,Scaffold_127,Scaffold_498,Scaffold_430,Scaffold_280,Scaffold_58

<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 ################################################################################


 ################################################################################
 karyotype =Scaffold_127SyntenicRibbonsKaryotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=Scaffold_127SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
chromosomes_order = Scaffold_159,Scaffold_308,Scaffold_287,Scaffold_206,Scaffold_281,Scaffold_489,Scaffold_472,Scaffold_420,Scaffold_30,Scaffold_227,Scaffold_53,Scaffold_327,Scaffold_127,Scaffold_32,Scaffold_5,Scaffold_584,Scaffold_461,Scaffold_172,Scaffold_163,Scaffold_321,Scaffold_573,Scaffold_101,Scaffold_565,Scaffold_148,Scaffold_572,Scaffold_538

<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 ################################################################################


 ################################################################################
 karyotype =Scaffold_172SyntenicRibbonsKaryotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=Scaffold_172SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
chromosomes_order = Scaffold_127,Scaffold_138,Scaffold_485,Scaffold_575,Scaffold_5,Scaffold_175,Scaffold_243,Scaffold_140,Scaffold_101,Scaffold_486,Scaffold_573,Scaffold_329,Scaffold_522,Scaffold_172,Scaffold_117,Scaffold_159,Scaffold_583,Scaffold_104,Scaffold_551,Scaffold_218


<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 ################################################################################


 ################################################################################
 karyotype =Scaffold_413SyntenicRibbonsKaryotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=Scaffold_413SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
chromosomes_order = Scaffold_193,Scaffold_413,Scaffold_6,Scaffold_310,Scaffold_34,Scaffold_259


<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 ################################################################################

 ################################################################################
 karyotype =Scaffold_474SyntenicRibbonsKaryotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=Scaffold_474SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
chromosomes_order = Scaffold_474,Scaffold_99,Scaffold_197,Scaffold_206,Scaffold_159,Scaffold_45,Scaffold_17,Scaffold_575,Scaffold_572,Scaffold_210,Scaffold_584,Scaffold_583,Scaffold_343,Scaffold_545,Scaffold_5,Scaffold_573,Scaffold_101,Scaffold_238


<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 ################################################################################


 ################################################################################
 karyotype =Scaffold_572SyntenicRibbonsKaryotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=Scaffold_572SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
chromosomes_order = Scaffold_572,Scaffold_474,Scaffold_584,Scaffold_406,Scaffold_558,Scaffold_390,Scaffold_518,Scaffold_18,Scaffold_573,Scaffold_159,Scaffold_560,Scaffold_441,Scaffold_127,Scaffold_32,Scaffold_96,Scaffold_295,Scaffold_5,Scaffold_73,Scaffold_291,Scaffold_11,Scaffold_359,Scaffold_101,Scaffold_430,Scaffold_220,Scaffold_59,Scaffold_583


<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 ################################################################################


 ################################################################################
 karyotype =Scaffold_573SyntenicRibbonsKaryotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file= Scaffold_573SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
chromosomes_order = Scaffold_573,Scaffold_14,Scaffold_559,Scaffold_572,Scaffold_235,Scaffold_75,Scaffold_317,Scaffold_11,Scaffold_581,Scaffold_103,Scaffold_534,Scaffold_390,Scaffold_216,Scaffold_518,Scaffold_127,Scaffold_77,Scaffold_214,Scaffold_555,Scaffold_126,Scaffold_192,Scaffold_9,Scaffold_281,Scaffold_351,Scaffold_82,Scaffold_485,Scaffold_38,Scaffold_78,Scaffold_228,Scaffold_531,Scaffold_74,Scaffold_252,Scaffold_287,Scaffold_101,Scaffold_461,Scaffold_221,Scaffold_484,Scaffold_428,Scaffold_324,Scaffold_172,Scaffold_577,Scaffold_423,Scaffold_315,Scaffold_249,Scaffold_284,Scaffold_244,Scaffold_230,Scaffold_583,Scaffold_582,Scaffold_292,Scaffold_159,Scaffold_584,Scaffold_36,Scaffold_474,Scaffold_466,Scaffold_32,Scaffold_329,Scaffold_346,Scaffold_66,Scaffold_433,Scaffold_27

<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 ################################################################################


 ################################################################################
 karyotype = Scaffold_575SyntenicRibbonsKaryotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=Scaffold_575SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
chromosomes_order = Scaffold_575,Scaffold_142,Scaffold_159,Scaffold_266,Scaffold_417,Scaffold_5,Scaffold_308,Scaffold_475,Scaffold_474,Scaffold_59,Scaffold_101,Scaffold_121,Scaffold_98,Scaffold_550,Scaffold_104,Scaffold_32,Scaffold_172


<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 karyotype =Scaffold_583SyntenicRibbonsKaryotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=Scaffold_583SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>



<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
chromosomes_order = Scaffold_583,Scaffold_452,Scaffold_573,Scaffold_305,Scaffold_572,Scaffold_318,Scaffold_273,Scaffold_101,Scaffold_334,Scaffold_474,Scaffold_385,Scaffold_237,Scaffold_584,Scaffold_172,Scaffold_346,Scaffold_28


<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 ################################################################################


 ################################################################################
 karyotype =Scaffold_584SyntenicRibbonsKaryotype.conf
chromosomes_units = 100000
  <<include ideogram.conf>>
  <<include ticks.conf>>
  <<include bands.conf>>

  <links>
  <link>
    file=Scaffold_584SyntenicRibbons.conf
    radius = 0.94r
    bezier_radius = 0.1r
    thickness = 1
    ribbon = yes
  </link>
  </links>
<image>
  <<include /shared/software/GIF/programs/circos/0.69.2/etc/image.conf>>
angle_offset* = -46
</image>
chromosomes_order = Scaffold_539,Scaffold_583,Scaffold_461,Scaffold_188,Scaffold_415,Scaffold_142,Scaffold_284,Scaffold_573,Scaffold_484,Scaffold_474,Scaffold_159,Scaffold_310,Scaffold_425,Scaffold_584,Scaffold_90,Scaffold_165,Scaffold_230,Scaffold_492,Scaffold_559,Scaffold_572,Scaffold_453,Scaffold_74,Scaffold_273,Scaffold_127,Scaffold_582
<<include /shared/software/GIF/programs/circos/0.69.2/etc/colors_fonts_patterns.conf>>
 <<include ./housekeeping.conf>>
 ################################################################################

```
