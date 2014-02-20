{    
  if(gSystem->Getenv("CMSSW_VERSION")) {    
    gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
    gSystem->AddIncludePath("-I$CMSSW_RELEASE_BASE/src/");
    gSystem->AddIncludePath("-I$CMSSW_BASE/src/BaconAna/DataFormats/interface");
    
    gInterpreter->AddIncludePath(TString(gSystem->Getenv("CMSSW_BASE"))+"/src/");
    gInterpreter->AddIncludePath(TString(gSystem->Getenv("CMSSW_RELEASE_BASE"))+"/src/");
    gInterpreter->AddIncludePath(TString(gSystem->Getenv("CMSSW_BASE"))+"/src/BaconAna/DataFormats/interface");
    
    gSystem->Load(Form("%s/lib/%s/%s",
                       gSystem->Getenv("CMSSW_BASE"),
                       gSystem->Getenv("SCRAM_ARCH"),
		       "libBaconAnaDataFormats.so"));
  }
               
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
