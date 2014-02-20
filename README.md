Bacon
=====

Bacon makes it Better

scramv1 project CMSSW CMSSW_5_3_14_patch2
cd CMSSW_5_3_14_patch2/src/
eval `scramv1 runtime -sh

git clone https://github.com/violatingcp/Bacon.git
mv Bacon/* .
rm -rf Bacon

Modify $CMSSW_BASE/src/BaconProd/scripts/setup_cvs.sh (to have your user id)
     => Take a look in the script and make sure it can do what you want

source $CMSSW_BASE/src/BaconProd/scripts/setup_cvs.sh
scramv1 b -j20

Now to run a test

cd $CMSSW_BASE/src/BaconProd/Ntupler/python/

you will see 3 example configs

cmsRun makingBacon_MC.py   
cmsRun makingBacon_Data.py   
cmsRun makingExpert_MC.py

Bacon == Student ntuples
Expert == Expert ntuples
