#include "../include/JetLoader.hh"
#include <cmath>
#include <iostream>

using namespace baconhep;

JetLoader::JetLoader(TTree *iTree,std::string iHLTFile) { 
  fJets     = new TClonesArray("baconhep::TJet");
  fAddJets  = new TClonesArray("baconhep::TAddJet");
  iTree->SetBranchAddress("Jet05",       &fJets);
  iTree->SetBranchAddress("AddJet05",    &fAddJets);
  fJetBr      = iTree->GetBranch("Jet05");
  fAddJetBr   = iTree->GetBranch("AddJet05");

  fTrigger = new TTrigger(iHLTFile);

  fPtr1    = new TLorentzVector();
  fPtr2    = new TLorentzVector();
  fPtr3    = new TLorentzVector();
  fPtr4    = new TLorentzVector();
  fDiJet   = new TLorentzVector();
}
JetLoader::~JetLoader() { 
  delete fJets;
  delete fJetBr;
  delete fAddJets;
  delete fAddJetBr;

  delete fPtr1;
  delete fPtr2; 
  delete fPtr3;
  delete fPtr4;
  delete fDiJet;
}
void JetLoader::reset(TJet &iJet,TAddJet &iAddJet) { 
  //Jet Properties => Probably a better way (but who cares for now)
  iJet.pt   = 0; 
  iJet.eta  = 0; 
  iJet.phi  = 0; 
  iJet.csv  = 0; 
  iJet.csv1 = 0; 
  iJet.csv2 = 0; 
  iJet.mva  = 0; 
  iJet.qgid = 0; 
  iJet.qg1  = 0; 
  iJet.qg2  = 0; 
  iJet.tau1 = 0; 
  iJet.tau2 = 0;  
  iJet.tau3 = 0; 
  iJet.tau4 = 0; 
  iJet.prunedm    = 0; 
  iJet.nCharged   = 0; 
  iJet.nNeutrals  = 0; 
  iJet.nParticles = 0; 
  iJet.beta       = 0; 
  iJet.betaStar   = 0; 
  iJet.dR2Mean    = 0; 
  iJet.ptD        = 0; 
  iJet.q          = 0; 
  iJet.pull       = 0; 
  iJet.pullAngle  = 0; 
  iJet.chEmFrac   = 0; 
  iJet.neuEmFrac  = 0; 
  iJet.chHadFrac  = 0; 
  iJet.neuHadFrac = 0; 
  iJet.mcFlavor   = 0; 
  iJet.mcFlavorPhys = 0; 
  iJet.genpt      = 0; 
  iJet.geneta     = 0; 
  iJet.genphi     = 0; 
  iJet.genm       = 0; 
  iAddJet.pt_p1     = 0; 
  iAddJet.phi_p1    = 0; 
  iAddJet.eta_p1    = 0; 
  iAddJet.mass_p1   = 0; 
  iAddJet.pt_t1     = 0; 
  iAddJet.phi_t1    = 0; 
  iAddJet.eta_t1    = 0; 
  iAddJet.mass_t1   = 0; 
}
void JetLoader::reset() { 
  fNJets      = 0;
  fNoiseClean = 0; 
  fNBTags     = 0; 
  fNBTags10   = 0; 
  fNQTags     = 0; 
  fJDEta      = 0; 
  fJDPhi      = 0; 
  
  fPtr1->SetPtEtaPhiM(1e-9,0,0,0);
  fPtr2->SetPtEtaPhiM(1e-9,0,0,0);
  fPtr3->SetPtEtaPhiM(1e-9,0,0,0);
  fPtr4->SetPtEtaPhiM(1e-9,0,0,0);
  fDiJet->SetPtEtaPhiM(1e-9,0,0,0);
  reset(fJet1,fAJet1); 
  reset(fJet2,fAJet2);
  reset(fJet3,fAJet3);
  reset(fJet4,fAJet4);
}
void JetLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("njets"        , &fNJets     , "fNJets/i");
  fTree->Branch("noiseCleaning", &fNoiseClean, "fNoiseClean/i");
  fTree->Branch("nbtags"       , &fNBTags    , "fNBTags/i");
  fTree->Branch("nbtags10"     , &fNBTags10  , "fNBTags10/i");
  fTree->Branch("nqtags"       , &fNQTags    , "fNQTags/i");
  fTree->Branch("jdeta"        , &fJDEta     , "fJDEta/F");
  fTree->Branch("jdphi"        , &fJDPhi     , "fJDPhi/F");
    
  fTree->Branch("jet1"         , "TLorentzVector", &fPtr1);
  fTree->Branch("jet2"         , "TLorentzVector", &fPtr2);
  fTree->Branch("jet3"         , "TLorentzVector", &fPtr3);
  fTree->Branch("jet4"         , "TLorentzVector", &fPtr4);
  fTree->Branch("dijet"        , "TLorentzVector", &fDiJet);
  //Jet Properites
  fTree->Branch("jet1CHF"     , &fJet1.chEmFrac     , "fJet1CHF/F");
  fTree->Branch("jet1NHF"     , &fJet1.neuHadFrac   , "fJet1NHF/F");
  fTree->Branch("jet1NEMF"    , &fJet1.neuEmFrac    , "fJet1NEMF/F");
  fTree->Branch("jet1Unc"     , &fJet1.unc          , "fJet1Unc/F");
  fTree->Branch("jet1Btag"    , &fJet1.csv          , "fJet1Btag/F");
  fTree->Branch("jet1QGtag"   , &fJet1.qgid         , "fJet1QGtag/F");
  fTree->Branch("jet1PartonId", &fJet1.mcFlavorPhys , "fJet1PartonId/i");
  fTree->Branch("jet1tau1"    , &fJet1.tau1         , "fJet1Tau1/F");
  fTree->Branch("jet1tau2"    , &fJet1.tau2         , "fJet1Tau2/F");
  fTree->Branch("jet1tau3"    , &fJet1.tau3         , "fJet1Tau3/F");
  fTree->Branch("jet1tau4"    , &fJet1.tau4         , "fJet1Tau4/F");
  fTree->Branch("jet1mtrim"   , &fAJet1.mass_t1     , "fJet1MTrim/F");
  fTree->Branch("jet1mprune"  , &fJet1.prunedm      , "fJet1MPrune/F");
  fTree->Branch("jet1npart"   , &fJet1.nParticles   , "fJet1NPart/F");
  fTree->Branch("jet1csv1"    , &fJet1.csv1         , "fJet1csv1/F");
  fTree->Branch("jet1csv2"    , &fJet1.csv2         , "fJet1csv2/F");
  fTree->Branch("jet1qg1"     , &fJet1.qg1          , "fJet1qg1/F");
  fTree->Branch("jet1qg2"     , &fJet1.qg2          , "fJet1qg2/F");
  fTree->Branch("jet1q"       , &fJet1.q            , "fJet1q/F");
  fTree->Branch("jet1pull"    , &fJet1.pull         , "fJet1pull/F");

  fTree->Branch("jet2CHF"     , &fJet2.chEmFrac     , "fJet2CHF/F");
  fTree->Branch("jet2NHF"     , &fJet2.neuHadFrac   , "fJet2NHF/F");
  fTree->Branch("jet2NEMF"    , &fJet2.neuEmFrac    , "fJet2NEMF/F");
  fTree->Branch("jet2Unc"     , &fJet2.unc          , "fJet2Unc/F");
  fTree->Branch("jet2Btag"    , &fJet2.csv          , "fJet2Btag/F");
  fTree->Branch("jet2QGtag"   , &fJet2.qgid         , "fJet2QGtag/F");
  fTree->Branch("jet2PartonId", &fJet2.mcFlavorPhys , "fJet2PartonId/i");
  fTree->Branch("jet2tau1"    , &fJet2.tau1         , "fJet2Tau1/F");
  fTree->Branch("jet2tau2"    , &fJet2.tau2         , "fJet2Tau2/F");
  fTree->Branch("jet2tau3"    , &fJet2.tau3         , "fJet2Tau3/F");
  fTree->Branch("jet2tau4"    , &fJet2.tau4         , "fJet2Tau4/F");
  fTree->Branch("jet2mtrim"   , &fAJet2.mass_t1     , "fJet2MTrim/F");
  fTree->Branch("jet2mprune"  , &fJet2.prunedm      , "fJet2MPrune/F");
  fTree->Branch("jet2npart"   , &fJet2.nParticles   , "fJet2NPart/F");
  fTree->Branch("jet2csv1"    , &fJet2.csv1         , "fJet2csv1/F");
  fTree->Branch("jet2csv2"    , &fJet2.csv2         , "fJet2csv2/F");
  fTree->Branch("jet2qg1"     , &fJet2.qg1          , "fJet2qg1/F");
  fTree->Branch("jet2qg2"     , &fJet2.qg2          , "fJet2qg2/F");
  fTree->Branch("jet2q"       , &fJet2.q            , "fJet2q/F");
  fTree->Branch("jet2pull"    , &fJet2.pull         , "fJet2pull/F");

  fTree->Branch("jet3CHF"     , &fJet3.chEmFrac     , "fJet3CHF/F");
  fTree->Branch("jet3NHF"     , &fJet3.neuHadFrac   , "fJet3NHF/F");
  fTree->Branch("jet3NEMF"    , &fJet3.neuEmFrac    , "fJet3NEMF/F");
  fTree->Branch("jet3Unc"     , &fJet3.unc          , "fJet3Unc/F");
  fTree->Branch("jet3Btag"    , &fJet3.csv          , "fJet3Btag/F");
  fTree->Branch("jet3QGtag"   , &fJet3.qgid         , "fJet3QGtag/F");
  fTree->Branch("jet3PartonId", &fJet3.mcFlavorPhys , "fJet3PartonId/i");
  fTree->Branch("jet3tau1"    , &fJet3.tau1         , "fJet3Tau1/F");
  fTree->Branch("jet3tau2"    , &fJet3.tau2         , "fJet3Tau2/F");
  fTree->Branch("jet3tau3"    , &fJet3.tau3         , "fJet3Tau3/F");
  fTree->Branch("jet3tau4"    , &fJet3.tau4         , "fJet3Tau4/F");
  fTree->Branch("jet3mtrim"   , &fAJet3.mass_t1     , "fJet3MTrim/F");
  fTree->Branch("jet3mprune"  , &fJet3.prunedm      , "fJet3MPrune/F");
  fTree->Branch("jet3npart"   , &fJet3.nParticles   , "fJet3NPart/F");
  fTree->Branch("jet3csv1"    , &fJet3.csv1         , "fJet3csv1/F");
  fTree->Branch("jet3csv2"    , &fJet3.csv2         , "fJet3csv2/F");
  fTree->Branch("jet3qg1"     , &fJet3.qg1          , "fJet3qg1/F");
  fTree->Branch("jet3qg2"     , &fJet3.qg2          , "fJet3qg2/F");
  fTree->Branch("jet3q"       , &fJet3.q            , "fJet3q/F");
  fTree->Branch("jet3pull"    , &fJet3.pull         , "fJet3pull/F");

  fTree->Branch("jet4Btag"    , &fJet4.csv          , "fJet4Btag/F");
  fTree->Branch("jet4QGtag"   , &fJet4.qgid         , "fJet4QGtag/F");

  //Redundant Info (for TTree draw => maybe I get rid of this)
  /*
  fTree->Branch("jpt_1"  ,&fJet1.pt  ,"fPt1/F");
  fTree->Branch("jeta_1" ,&fJet1.eta ,"fEta1/F");
  fTree->Branch("jphi_1" ,&fJet1.phi ,"fPhi1/F");
  fTree->Branch("jm_1"   ,&fJet1.mass,"fM1/F");

  fTree->Branch("jpt_2"  ,&fJet2.pt  ,"fPt2/F");
  fTree->Branch("jeta_2" ,&fJet2.eta ,"fEta2/F");
  fTree->Branch("jphi_2" ,&fJet2.phi ,"fPhi2/F");
  fTree->Branch("jm_2"   ,&fJet2.mass,"fM2/F");

  fTree->Branch("jpt_3"  ,&fJet3.pt  ,"fPt3/F");
  fTree->Branch("jeta_3" ,&fJet3.eta ,"fEta3/F");
  fTree->Branch("jphi_3" ,&fJet3.phi ,"fPhi3/F");
  fTree->Branch("jm_3"   ,&fJet3.mass,"fM3/F");

  fTree->Branch("jpt_4"  ,&fJet4.pt  ,"fPt4/F");
  fTree->Branch("jeta_4" ,&fJet4.eta ,"fEta4/F");
  fTree->Branch("jphi_4" ,&fJet4.phi ,"fPhi4/F");
  fTree->Branch("jm_4"   ,&fJet4.mass,"fM4/F");
  */
}
void JetLoader::load(int iEvent) { 
  fJets     ->Clear();
  fAddJets  ->Clear();
  fJetBr    ->GetEntry(iEvent);
  fAddJetBr ->GetEntry(iEvent);
}
bool JetLoader::selectJets(std::vector<TLorentzVector> &iVetoes) {
  reset(); 
  std::vector<TJet*> lJets;
  for  (int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    TJet *pJet = (TJet*)((*fJets)[i0]);
    //Veto
    bool pMatch = false;
    for(unsigned int i1 = 0; i1 < iVetoes.size(); i1++) {
      double pDEta = pJet->eta      - iVetoes[i1].Eta();
      double pDPhi = fabs(pJet->phi - iVetoes[i1].Phi());
      if(fabs(pDPhi) > 2.*TMath::Pi()-fabs(pDPhi)) pDPhi =  2.*TMath::Pi()-fabs(pDPhi);
      if(sqrt(pDPhi*pDPhi+pDEta*pDEta) > 0.5) continue;
      pMatch = true;
    }
    if(pMatch) continue;
    if(passLoose(pJet)    && pJet->pt > 30 && passPUId(pJet)) fNJets++;
    if(pJet->csv  > 0.679 && pJet->pt > 20 && passPUId(pJet)) fNBTags++;    
    if(pJet->csv  > 0.679 && pJet->pt > 10 && passPUId(pJet)) fNBTags10++;    
    if(pJet->qgid > 0.4   && pJet->pt > 20 && passPUId(pJet)) fNQTags++;
    //Start filling the collection
    if(lJets.size() == 0) {lJets.push_back(pJet); continue;}
    //Vector Manipulation to get sort the jets in corrected pT (collection is sorted in raw pt)
    bool pFill = false;
    for( std::vector<TJet*>::iterator pJetIter = lJets.begin(); pJetIter != lJets.end(); pJetIter++) { 
      if((*pJetIter)->pt > pJet->pt) continue;
      lJets.insert(pJetIter,pJet);
      pFill = true;
      break;
    } 
    if(!pFill) lJets.push_back(pJet);
    //Limit this to the top 4 Jets
    //if(lJets.size() > 4) lJets.pop_back();
  }
  if(lJets.size() > 0)   fillVars(lJets[0],fPtr1,fJet1,fAJet1);
  if(lJets.size() > 1)   fillVars(lJets[1],fPtr2,fJet2,fAJet2);
  if(lJets.size() > 2)   fillVars(lJets[2],fPtr3,fJet3,fAJet3);
  if(lJets.size() > 3)   fillVars(lJets[3],fPtr4,fJet4,fAJet4);
  if(lJets.size() > 1)   {
    TLorentzVector lDiJet = *fPtr1 + *fPtr2;
    fDiJet->SetPtEtaPhiM(lDiJet.Pt(),lDiJet.Eta(),lDiJet.Phi(),lDiJet.M());
    fJDPhi = fPtr1->DeltaPhi(*fPtr2);
    fJDEta = fabs(fPtr1->Eta() - fPtr2->Eta());
  }
  if(lJets.size() > 0) fNoiseClean |= int(fJet1.chHadFrac >0.2) << 0;
  if(lJets.size() > 0) fNoiseClean |= int(fJet1.neuHadFrac>0.7) << 1;
  if(lJets.size() > 0) fNoiseClean |= int(fJet1.neuEmFrac >0.7) << 2;
  if(lJets.size() > 1) fNoiseClean |= int(fJet2.chHadFrac >0.2) << 3;
  if(lJets.size() > 1) fNoiseClean |= int(fJet2.neuHadFrac>0.7) << 4;
  if(lJets.size() > 1) fNoiseClean |= int(fJet2.neuEmFrac >0.7) << 5;
  fHLTMatch   = passTrigObj(&fJet1,0) << 0;
  if(lJets.size() > 1) fHLTMatch   = ((passTrigObj(&fJet1,1) && passTrigObj(&fJet2,1))) << 1;// ||  (passTrigObj(&fJet1,2) && passTrigObj(&fJet2,1)))  << 1;
  return true;
}
void JetLoader::fillVars(TJet *iJet,TLorentzVector *iPtr,TJet &iSaveJet,TAddJet &iASaveJet) { 
  iSaveJet.pt   = iJet->pt;
  iSaveJet.eta  = iJet->eta;
  iSaveJet.phi  = iJet->phi;
  iSaveJet.mass = iJet->mass;
  if(iJet->pt > 0) iPtr->SetPtEtaPhiM(iJet->pt,iJet->eta,iJet->phi,iJet->mass);
  TAddJet *lAJet = addJet(iJet);

  iSaveJet.chEmFrac     = iJet->chEmFrac;
  iSaveJet.neuHadFrac   = iJet->neuHadFrac;
  iSaveJet.neuEmFrac    = iJet->neuEmFrac;
  iSaveJet.unc          = iJet->unc;
  iSaveJet.csv          = iJet->csv;
  iSaveJet.qgid         = iJet->qgid;
  iSaveJet.mcFlavorPhys = iJet->mcFlavorPhys;
  iSaveJet.tau1         = iJet->tau1;
  iSaveJet.tau2         = iJet->tau2;
  iSaveJet.tau3         = iJet->tau3;
  iSaveJet.tau4         = iJet->tau4;
  iSaveJet.mass         = iJet->mass;
  iSaveJet.prunedm      = iJet->prunedm;
  iSaveJet.nParticles   = iJet->nParticles;
  iSaveJet.csv1         = iJet->csv1;
  iSaveJet.csv2         = iJet->csv2;
  iSaveJet.qg1          = iJet->qg1;
  iSaveJet.qg2          = iJet->qg2;
  iSaveJet.q            = iJet->q;
  iSaveJet.pull         = iJet->pull;
  iSaveJet.hltMatchBits = iJet->hltMatchBits;
  iASaveJet.pt_p1     = lAJet->pt_p1    ;
  iASaveJet.phi_p1    = lAJet->phi_p1   ;
  iASaveJet.eta_p1    = lAJet->eta_p1   ;
  iASaveJet.mass_p1   = lAJet->mass_p1  ;
  iASaveJet.pt_t1     = lAJet->pt_t1    ;
  iASaveJet.phi_t1    = lAJet->phi_t1   ;
  iASaveJet.eta_t1    = lAJet->eta_t1   ;
  iASaveJet.mass_t1   = lAJet->mass_t1  ;
}
bool JetLoader::vetoJet() {
  for  (int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    TJet *pJet = (TJet*)((*fJets)[i0]);
    if(passVeto(pJet)) return true;
  }
  return false;
}
//PFJet Ids
bool JetLoader::passLoose(TJet *iJet) { 
  if(iJet->neuEmFrac        >  0.99)                         return false;
  if(iJet->neuHadFrac       >  0.99)                         return false;
  if(iJet->nParticles       <  2)                            return false;
  if(iJet->chHadFrac        <= 0     && fabs(iJet->eta) < 2.4 )    return false;
  if(iJet->chEmFrac         >  0.99  && fabs(iJet->eta) < 2.4 )    return false;
  if(iJet->nCharged         < 1      && fabs(iJet->eta) < 2.4 )    return false;
  return true;
}
bool JetLoader::passTight(TJet *iJet) { 
  if(iJet->neuEmFrac        >  0.9)                          return false;
  if(iJet->neuHadFrac       >  0.9)                          return false;
  if(iJet->nParticles       <  2)                            return false;
  if(iJet->chHadFrac        <= 0     && fabs(iJet->eta) < 2.4 )    return false;
  if(iJet->chEmFrac         >  0.9   && fabs(iJet->eta) < 2.4 )    return false;
  if(iJet->nCharged         < 1      && fabs(iJet->eta) < 2.4 )    return false;
  return true;
}
bool JetLoader::passVeto(TJet *iJet) { 
  return passLoose(iJet);
}
bool JetLoader::passPUId(TJet *iJet) { 
  if(iJet->pt >  0 && iJet->pt < 20    && fabs(iJet->eta) < 2.5)                            return (iJet->mva > -0.95); 
  if(iJet->pt >  0 && iJet->pt < 20    && fabs(iJet->eta) > 2.5  && fabs(iJet->eta) < 2.75) return (iJet->mva > -0.96);
  if(iJet->pt >  0 && iJet->pt < 20    && fabs(iJet->eta) > 2.75 && fabs(iJet->eta) < 3.0 ) return (iJet->mva > -0.94);
  if(iJet->pt >  0 && iJet->pt < 20    && fabs(iJet->eta) > 3.0  )                          return (iJet->mva > -0.95);
  if(iJet->pt > 20 && iJet->pt < 20000 && fabs(iJet->eta) < 2.5)                            return (iJet->mva > -0.63); 
  if(iJet->pt > 20 && iJet->pt < 20000 && fabs(iJet->eta) > 2.5  && fabs(iJet->eta) < 2.75) return (iJet->mva > -0.60);
  if(iJet->pt > 20 && iJet->pt < 20000 && fabs(iJet->eta) > 2.75 && fabs(iJet->eta) < 3.0 ) return (iJet->mva > -0.55);
  if(iJet->pt > 20 && iJet->pt < 20000 && fabs(iJet->eta) > 3.0  )                          return (iJet->mva > -0.45);
  return false;
}
TAddJet *JetLoader::addJet(TJet *iJet) { 
  int lIndex = -1;
  TAddJet *lJet = 0; 
  for(int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    if((*fJets)[i0] == iJet) { lIndex = i0; break;}
  }
  if(lIndex == -1) return 0;
  for  (int i0 = 0; i0 < fAddJets->GetEntriesFast(); i0++) { 
    TAddJet *pJet = (TAddJet*)((*fAddJets)[i0]);
    if(pJet->index == fabs(lIndex)) { lJet = pJet; break;}
  }
  return lJet;
}
//Two triggers to be added
//MonoCentralPFJet80_PFMETnoMu
//HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v"
void JetLoader::addTrigger(std::string iName) { 
  fTrigString.push_back(iName);
}
bool JetLoader::passTrigObj(TJet *iJet,int iId) {
  bool lPass = false;
  if(fTrigger->passObj(fTrigString[iId],iJet->hltMatchBits))  lPass = true;
  return lPass;
}
