#!/bin/bash
#export CVSROOT=:ext:<cern-user-account>@lxplus.cern.ch:/afs/cern.ch/user/c/cvscmssw/public/CMSSW
export CVSROOT=:ext:pharris@lxplus.cern.ch:/afs/cern.ch/user/c/cvscmssw/public/CMSSW

CURRDIR=$PWD
PATCHDIR=$CMSSW_BASE/src/BaconProd/patch
cd $CMSSW_BASE/src

echo "Checking out MuScleFit package..."
cp -r /afs/cern.ch/work/p/pharris/public/tmp/CMSSW_5_3_13/src/MuScleFit $CMSSW_BASE/src

echo "Checking out Electron MVA ID package..."
cp -r /afs/cern.ch/work/p/pharris/public/tmp/CMSSW_5_3_13/src/EGamma  $CMSSW_BASE/src
mv EGamma/EGammaAnalysisTools/test/BuildFile.xml EGamma/EGammaAnalysisTools/test/BuildFile.xmlSilent
mv EGamma/EGammaAnalysisTools/plugins/BuildFile.xml EGamma/EGammaAnalysisTools/plugins/BuildFile.xmlSilent
cp $PATCHDIR/EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h.53Xpatch EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h

echo "Checking out packages for MET filters..."
cvs co -r V00-00-13-01 RecoMET/METFilters
cvs co -r V00-00-08    RecoMET/METAnalyzers
cvs co -r V00-03-23    CommonTools/RecoAlgos
#cvs co -r V01-00-11-01 DPGAnalysis/Skims
#cvs co -r V00-11-17    DPGAnalysis/SiStripTools
#cvs co -r V00-00-08    DataFormats/TrackerCommon
#cvs co -r V01-09-05    RecoLocalTracker/SubCollectionProducers

git clone https://github.com/nhanvtran/JetTools.git
cp /afs/cern.ch/work/p/pharris/public/tmp/CMSSW_5_3_13/src/JetTools/AnalyzerToolbox/python/njettinessadder_cfi.py JetTools/AnalyzerToolbox/python/

git cms-merge-topic -u TaiSakuma:53X-met-131120-01
git cms-merge-topic -u cms-analysis-tools:5_3_14-updateSelectorUtils
git cms-merge-topic -u cms-analysis-tools:5_3_13_patch2-testNewTau
git cms-merge-topic -u cms-met:53X-MVaNoPuMET-20131217-01
# Remove a bunch of directories bcecause somebody made a poor Met Recipe
rm -rf PhysicsTools
rm -rf GeneratorInterface
rm -rf FWCore
rm -rf RecoBTag

git clone https://github.com/violatingcp/Jets_Short.git
mv Jets_Short/RecoJets/JetProducers/data/*.xml RecoJets/JetProducers/data/
rm -rf Jets_Short

cp /afs/cern.ch/work/p/pharris/public/tmp/CMSSW_5_3_13/src/RecoMET/METPUSubtraction/python/mvaPFMET_leptons_cff.py $CMSSW_BASE/src/RecoMET/METPUSubtraction/python/
cp /afs/cern.ch/cern.ch/p/pharris/public/MVAMetUpdate/*Sep*.root                                                   $CMSSW_BASE/src/RecoMET/METPUSubtraction/data/


#echo "Checking out packages for MET..."
#cvs co -r ??? JetMETCorrections/METPUSubtraction
#cvs co -r ??? DataFormats/JetReco
#cvs co -r ??? DataFormats/METReco
#cvs co -r ??? RecoJets/JetProducers

cd $CURRDIR
