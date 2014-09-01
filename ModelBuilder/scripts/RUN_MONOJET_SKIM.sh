#!/bin/bash

# SUBMIT SKIMMING/Selection jobs for the whole of the Mono-X samples
# arguments for skimmer are 
#	0 -> index
# 	1 -> total events in sample ( pass from baconBatch with --passSumEntries "1:Tree" )
#	2 -> x-section
#	3 -> path to input file 

RUNLEP=0
RUN2L=1
RUNSIG=1

EXEC=skimMetFilesExecDir
OUTP=Output_tree.root
OPTS='--passSumEntries 1:Tree --blacklist Output'
OUTDIR=BaconMonoX

  # Single Lepton Region ....
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:24" -a "2:3354."   -n 5 -o BaconMonoX/DYLLL -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/DYL"

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:25" -a "2:30400.0" -n 5 -o BaconMonoX/WL    -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/WL/"  

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:26" -a "2:48.01"   -n 2 -o BaconMonoX/WHTL_250300 -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/W250300L/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:26" -a "2:38.3"    -n 2 -o BaconMonoX/WHTL_300400 -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/W300400L/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:26" -a "2:25.22 "  -n 2 -o BaconMonoX/WHTL_400inf -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/W400infL"

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:27" -a "2:109.281" -n 2 -o BaconMonoX/TTL_sem     -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/TTSemL/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:27" -a "2:26.1975" -n 2 -o BaconMonoX/TTL_lep     -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/TTLepL/"

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:28" -a "2:3.89"    -n 2 -o BaconMonoX/TL_s -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/TsL/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:28" -a "2:55.531"  -n 2 -o BaconMonoX/TL_t -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/TtL/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:28" -a "2:11.1"    -n 2 -o BaconMonoX/TL_W -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/TWL/"

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:22" -a "2:9.03"    -n 3 -o BaconMonoX/ZZL -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/ZZL/" 

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:29" -a "2:57.2"    -n 3 -o BaconMonoX/WWL -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/WWL/" 

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:23" -a "2:22.44"   -n 3 -o BaconMonoX/WZL -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/WZL/" 

  python baconBatch.py $EXEC $OUTP    --blacklist "Output"    -a "0:30" -a "1:-1" -a "2:-1" -n 5  -o BaconMonoX/Mu2012L -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/Mu2012AL/,/store/group/phys_jetmet/pharris/production/01/ntuples/Mu2012BL/,/store/group/phys_jetmet/pharris/production/01/ntuples/Mu2012CL/,/store/group/phys_jetmet/pharris/production/01/ntuples/Mu2012DL/"

 # Di-lepton Region
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:14" -a "2:3354."   -n 5 -o BaconMonoX/DYLLLL -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/DYL"

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:17" -a "2:109.281" -n 2 -o BaconMonoX/TTLL_sem     -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/TTSemL/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:17" -a "2:26.1975" -n 2 -o BaconMonoX/TTLL_lep     -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/TTLepL/"

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:16" -a "2:3.89"    -n 2 -o BaconMonoX/TLL_s -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/TsL/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:16" -a "2:55.531"  -n 2 -o BaconMonoX/TLL_t -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/TtL/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:16" -a "2:11.1"    -n 2 -o BaconMonoX/TLL_W -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/TWL/"

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:20" -a "2:9.03"    -n 3 -o BaconMonoX/ZZLL -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/ZZL/" 

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:18" -a "2:57.2"    -n 3 -o BaconMonoX/WWLL -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/WWL/" 

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:19" -a "2:22.44"   -n 3 -o BaconMonoX/WZLL -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/WZL/" 

 python baconBatch.py $EXEC $OUTP   --blacklist "Output"     -a "0:21" -a "1:-1" -a "2:-1" -n 5  -o BaconMonoX/Mu2012LL -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/Mu2012AL/,/store/group/phys_jetmet/pharris/production/01/ntuples/Mu2012BL/,/store/group/phys_jetmet/pharris/production/01/ntuples/Mu2012CL/,/store/group/phys_jetmet/pharris/production/01/ntuples/Mu2012DL/"

 # MET (Signal) Region

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:0" -a "2:381.2"  -n 3 -o BaconMonoX/DYNN_50100  -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/DYNN050100/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:0" -a "2:160.3"  -n 3 -o BaconMonoX/DYNN_100200 -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/DYNN100200/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:0" -a "2:41.49"  -n 3 -o BaconMonoX/DYNN_200400 -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/DYNN200400/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:0" -a "2:5.274"  -n 3 -o BaconMonoX/DYNN_400inf -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/DYNN400inf/"

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:2" -a "2:30400.0" -n 5 -o BaconMonoX/W 	     -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/W/"

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:3"  -a "2:48.01" -n 2 -o BaconMonoX/WHT_250300 -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/W250300/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:3"  -a "2:38.3"  -n 2 -o BaconMonoX/WHT_300400 -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/W300400/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:3"  -a "2:25.22" -n 2 -o BaconMonoX/WHT_400inf -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/W400inf/"

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:4" -a "2:109.281" -n 2 -o BaconMonoX/TT_sem     -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/TTSem/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:4" -a "2:26.1975" -n 2 -o BaconMonoX/TT_lep     -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/TTLep/"

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:5" -a "2:3.89"    -n 2 -o BaconMonoX/T_s -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/Ts/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:5" -a "2:55.531"  -n 2 -o BaconMonoX/T_t -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/Tt/"
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:5" -a "2:11.1"    -n 2 -o BaconMonoX/T_W -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/TW/"

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:6" -a "2:9.03"    -n 3 -o BaconMonoX/ZZ -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/ZZ/" 

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:7" -a "2:57.2"    -n 3 -o BaconMonoX/WW -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/WW/" 

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:8" -a "2:22.44"   -n 3 -o BaconMonoX/WZ -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/WZ/" 

  python baconBatch.py $EXEC $OUTP $OPTS -a "0:11" -a "2:19.27 "   -n 1 -o BaconMonoX/GGH -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/GGH/" 
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:12" -a "2:1.578 "   -n 1 -o BaconMonoX/VBF -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/VBF/" 
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:9"  -a "2:0.7046"   -n 1 -o BaconMonoX/WH -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/WH/" 
  python baconBatch.py $EXEC $OUTP $OPTS -a "0:10" -a "2:0.4153"   -n 1 -o BaconMonoX/ZH -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/ZH/" 

  python baconBatch.py $EXEC $OUTP --blacklist "Output"   -a "0:13" -a "1:-1" -a "2:-1" -n 5  -o BaconMonoX/MET -d "3:/store/group/phys_jetmet/pharris/production/01/ntuples/META/,/store/group/phys_jetmet/pharris/production/01/ntuples/METB/,/store/group/phys_jetmet/pharris/production/01/ntuples/METC/,/store/group/phys_jetmet/pharris/production/01/ntuples/METD/"

  #python baconBatch.py $EXEC $OUTP $OPTS -a "1:1" -n 5 -o BaconMonoX/DYLL -d "0:/store/group/phys_jetmet/pharris/production/01/ntuples/DY"   -> Do we not need this guy?
