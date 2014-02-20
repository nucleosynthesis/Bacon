#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <utility>
#include <iomanip>

#include "Lepton.hh"

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#endif

#define ELE_REFERENCE_IDMVA_CUT_BIN0  0.470   // pT<10, |eta|<0.8
#define ELE_REFERENCE_IDMVA_CUT_BIN1  0.004   // pT<10, 0.8<|eta|<1.479
#define ELE_REFERENCE_IDMVA_CUT_BIN2  0.295   // pT<10, |eta|>1.479
#define ELE_REFERENCE_IDMVA_CUT_BIN3 -0.340   // pT>10, |eta|<0.8
#define ELE_REFERENCE_IDMVA_CUT_BIN4 -0.650   // pT>10, 0.8<|eta|<1.479
#define ELE_REFERENCE_IDMVA_CUT_BIN5  0.600   // pT>10, |eta|>1.479

#define JET_MVA_CUT0 -0.63
#define JET_MVA_CUT1 -0.60
#define JET_MVA_CUT2 -0.55
#define JET_MVA_CUT3 -0.45

#define ELECTRON_PDGID 11
#define MUON_PDGID     13

#define MUON_MASS      0.105658369
#define ELECTRON_MASS  0.000511
#define Z_MASS        91.1876

#define Z1_MASS_MIN  40
#define Z1_MASS_MAX 120
#define Z2_MASS_MIN   4
#define Z2_MASS_MAX 120

#define ISO_CUT 0.4

using namespace baconhep;
using namespace std;

bool passMuonID(const TMuon *muon);
bool passEleID(const TElectron *ele);
bool passJetID(const TJet *jet);
bool isBetterMuon(const TMuon *mu1, const TMuon *mu2);
double getEffArea(const double eta);
double computeIso(const Lepton &lep, const TPhoton *fsrCand, const double rho);
const TPhoton* recoverFsr(const Lepton &lep, const vector<Lepton> &lepvec, const pair<Lepton,Lepton> &zcand, 
                          const TClonesArray *photonArr);

//--------------------------------------------------------------------------------------------------
void testHZZ4l()
{
  TFile *infile = new TFile("../python/ntuple.root");
  TTree *intree = (TTree*)infile->Get("Events");
  
  TEventInfo   *info      = new TEventInfo();
  TClonesArray *eleArr    = new TClonesArray("baconhep::TElectron");
  TClonesArray *muonArr   = new TClonesArray("baconhep::TMuon");
  TClonesArray *photonArr = new TClonesArray("baconhep::TPhoton");
  TClonesArray *jetArr    = new TClonesArray("baconhep::TJet");
  
  intree->SetBranchAddress("Info",    &info);      TBranch *infoBr   = intree->GetBranch("Info");
  intree->SetBranchAddress("Electron",&eleArr);    TBranch *eleBr    = intree->GetBranch("Electron");
  intree->SetBranchAddress("Muon",    &muonArr);   TBranch *muonBr   = intree->GetBranch("Muon");
  intree->SetBranchAddress("Photon",  &photonArr); TBranch *photonBr = intree->GetBranch("Photon");
  intree->SetBranchAddress("Jet",     &jetArr);    TBranch *jetBr    = intree->GetBranch("Jet");
  
  unsigned int n4l=0, n4e=0, n4mu=0, n2e2mu=0;
  
  for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {
    infoBr->GetEntry(ientry);


    //===========================================================
    // Trigger Selection  
    //-----------------------------------------------------------
    ULong64_t triggerBits = kHLT_Mu17_Mu8 | kHLT_Mu17_TkMu8 | 
                            kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL |
			    kHLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL |
			    kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL |
			    kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
    if(!(info->triggerBits & triggerBits)) continue;
   
    
    //===========================================================
    // Primary Vertex Requirement 
    //-----------------------------------------------------------  
    if(!info->hasGoodPV) continue;
  
  
    //===========================================================
    // Lepton Selection 
    //-----------------------------------------------------------
    vector<Lepton> goodLeptons;  // leptons passing ID and IP cuts (before iso cut)
    vector<Lepton> looseMuons;   // muons for lepton cross cleaning
  
    //
    // Loop over muons
    //
    muonArr->Clear();
    muonBr->GetEntry(ientry);
    for(int i=0; i<muonArr->GetEntriesFast(); ++i) {
      const TMuon *muon = (TMuon*)muonArr->At(i);  
    
      Lepton mulep(muon, (muon->q<0 ? MUON_PDGID : -MUON_PDGID));
      mulep.p4.SetPtEtaPhiM(muon->ptHZZ4l, muon->eta ,muon->phi, MUON_MASS);

      if(mulep.p4.Pt()        < 5)   continue;  // pT cut
      if(fabs(mulep.p4.Eta()) > 2.4) continue;  // |eta| cut

      if((muon->selectorBits & kGlobal) || (muon->typeBits & kPFMuon))
        looseMuons.push_back(mulep);

      if(passMuonID(muon) && fabs(muon->sip3d)<4) {    
        // Ghost muon removal
        // NOTE: goodLeptons array only contains muons at this point...
        int ghostindex = -1;
        for(unsigned int j=0; j<goodLeptons.size(); j++) {
          if(mulep.p4.DeltaR(goodLeptons[j].p4)<0.02)
            ghostindex = j;
        }
        if(ghostindex>=0) {
          if(isBetterMuon((TMuon*)(mulep.baconObj), (TMuon*)(goodLeptons[ghostindex].baconObj)))
            goodLeptons[ghostindex] = mulep;
    
        } else {
          goodLeptons.push_back(mulep);
        }    
      }
    }
    
    //
    // Loop over electrons
    //
    eleArr->Clear();
    eleBr->GetEntry(ientry);                
    for(int i=0; i<eleArr->GetEntriesFast(); ++i) {
      const TElectron *ele = (TElectron*)eleArr->At(i);  
    
      Lepton elelep(ele, (ele->q<0 ? ELECTRON_PDGID : -ELECTRON_PDGID));
      elelep.p4.SetPtEtaPhiM(ele->ptHZZ4l, ele->eta, ele->phi, ELECTRON_MASS);
    
      if(elelep.p4.Pt()        <  7)   continue;  // pT cut
      if(fabs(elelep.p4.Eta()) >= 2.5) continue;  // |eta| cut
  
      // Lepton Cross-Cleaning: remove electrons close to loose muons 
      bool isMuon=false;
      for(unsigned int imu=0; imu<looseMuons.size(); imu++) {
        if(elelep.p4.DeltaR(looseMuons[imu].p4)<0.05) {
          isMuon = true;
	  break;
        }
      }
      if(isMuon) continue;
   
      if(passEleID(ele) && fabs(ele->sip3d)<4) {
        goodLeptons.push_back(elelep);
      }
    }


    //========================================================
    // Z candidates preselection
    //--------------------------------------------------------
    photonArr->Clear();
    photonBr->GetEntry(ientry);
    
    vector< pair<Lepton,Lepton> > dilep;
    for(unsigned int i=0; i<goodLeptons.size(); i++) {
      for(unsigned int j=i+1; j<goodLeptons.size(); j++) {
        if(abs(goodLeptons[i].pdgId) != abs(goodLeptons[j].pdgId)) continue;
        if(goodLeptons[i].q          == goodLeptons[j].q)          continue;
      
        pair<Lepton,Lepton> zcand( (goodLeptons[i].pdgId<0) ? goodLeptons[i] : goodLeptons[j],
                                   (goodLeptons[i].pdgId<0) ? goodLeptons[j] : goodLeptons[i] );

        const TPhoton *fsrCand1 = recoverFsr(zcand.first,  goodLeptons, zcand, photonArr);
	TLorentzVector fsrvec1;
	if(fsrCand1) fsrvec1.SetPtEtaPhiM(fsrCand1->pfPt, fsrCand1->pfEta, fsrCand1->pfPhi, 0);

        const TPhoton *fsrCand2 = recoverFsr(zcand.second, goodLeptons, zcand, photonArr);
	TLorentzVector fsrvec2;
	if(fsrCand2) fsrvec2.SetPtEtaPhiM(fsrCand2->pfPt, fsrCand2->pfEta, fsrCand2->pfPhi, 0);	

	bool use1=false, use2=false;
        if(fsrCand1) {
          if(!fsrCand2) {
	    use1 = true;
	  } else {
	    // Both leptons have associated FSR photons
	    // Pick the one photon that brings the mass closer to Z mass
	    double newMass1 = (zcand.first.p4 + zcand.second.p4 + fsrvec1).M();
	    double newMass2 = (zcand.first.p4 + zcand.second.p4 + fsrvec2).M();
	    if(fabs(newMass1 - Z_MASS) < fabs(newMass2 - Z_MASS)) {
	      use1 = true;
	    } else {
	      use2 = true;
	    }
  	  }
        } else if(fsrCand2) {
	  use2 = true;
        }
		
	// Correct the isolation by removing FSR footprint and apply isolation cut
	if(use1) zcand.first.fsrp4  = fsrvec1;
	zcand.first.iso = computeIso(zcand.first, fsrCand1, info->rhoIso);
	if(zcand.first.iso > ISO_CUT*(zcand.first.p4.Pt())) continue;
        
	if(use2) zcand.second.fsrp4 = fsrvec2;
	zcand.second.iso = computeIso(zcand.second, fsrCand2, info->rhoIso);
	if(zcand.second.iso > ISO_CUT*(zcand.second.p4.Pt())) continue;

        dilep.push_back(zcand);
      }
    }


    //==============================================================================
    // Require Z1 with  40<m<120
    //------------------------------------------------------------------------------
    double best_dm = 999;
    pair<Lepton,Lepton> z1dilep;
    TLorentzVector z1vec(0,0,0,0);
    for(unsigned int k=0; k<dilep.size(); k++) {
      double mass = (dilep[k].first.p4 + dilep[k].second.p4 + dilep[k].first.fsrp4 + dilep[k].second.fsrp4).M();
      double dm   = fabs(mass - Z_MASS);
      if(dm < best_dm) {
        best_dm = dm;
        z1dilep = dilep[k];
        z1vec   = dilep[k].first.p4 + dilep[k].second.p4 + dilep[k].first.fsrp4 + dilep[k].second.fsrp4;
      }
    }
    if(z1vec.M()<=Z1_MASS_MIN || z1vec.M()>=Z1_MASS_MAX) continue;

  
    //==============================================================================
    // Require Z2 with  4<m<120
    //------------------------------------------------------------------------------
    double best_sumpt = 0;
    pair<Lepton,Lepton> z2dilep;
    TLorentzVector z2vec(0,0,0,0);
    for(unsigned int k=0; k<dilep.size(); k++) {
      // Z2 cannot share leptons with Z1
      if(z1dilep.first.baconObj  == dilep[k].first.baconObj)  continue;
      if(z1dilep.second.baconObj == dilep[k].second.baconObj) continue;
      //double mass = (dilep[k].first.p4 + dilep[k].second.p4 + dilep[k].first.fsrp4 + dilep[k].second.fsrp4).M();
      //if(mass<=Z2_MASS_MIN || mass>=Z2_MASS_MAX) continue;
      double sumpt = dilep[k].first.p4.Pt() + dilep[k].second.p4.Pt();
      if(sumpt > best_sumpt) {
        best_sumpt = sumpt;
        z2dilep    = dilep[k];
        z2vec      = dilep[k].first.p4 + dilep[k].second.p4 + dilep[k].first.fsrp4 + dilep[k].second.fsrp4;
      }
    }
    if(z2vec.M()<=Z2_MASS_MIN || z2vec.M()>=Z2_MASS_MAX) continue; 


    //******************************************************************************
    // Any two leptons 20/10
    //******************************************************************************
    int nlepPt10=0, nlepPt20=0;
    if     (z1dilep.first.p4.Pt()  > 20) { nlepPt10++; nlepPt20++; }
    else if(z1dilep.first.p4.Pt()  > 10) { nlepPt10++; }
    if     (z1dilep.second.p4.Pt() > 20) { nlepPt10++; nlepPt20++; }
    else if(z1dilep.second.p4.Pt() > 10) { nlepPt10++; }
    if     (z2dilep.first.p4.Pt()  > 20) { nlepPt10++; nlepPt20++; }
    else if(z2dilep.first.p4.Pt()  > 10) { nlepPt10++; }
    if     (z2dilep.second.p4.Pt() > 20) { nlepPt10++; nlepPt20++; }
    else if(z2dilep.second.p4.Pt() > 10) { nlepPt10++; }
    if(nlepPt20<1) continue;
    if(nlepPt10<2) continue;

    //******************************************************************************
    // Low mass resonance veto: any opp. sign pairs with mass <= 4 GeV
    //******************************************************************************  
    if( (z1dilep.first.p4  + z2dilep.second.p4).M() <= 4) continue;
    if( (z1dilep.second.p4 + z2dilep.first.p4 ).M() <= 4) continue;
  
    //******************************************************************************
    // Sync cuts for m(Z2), m(4l)
    //******************************************************************************
    if(z2vec.M() <= 12) continue;
  
    TLorentzVector zzvec = z1vec + z2vec;  
    if(zzvec.M() <= 70) continue;

    //
    // loop through jets
    //
    unsigned int njets=0;
    TLorentzVector jet1, jet2;
    jetArr->Clear();
    jetBr->GetEntry(ientry);
    for(int i=0; i<jetArr->GetEntriesFast(); i++) {
      const TJet *jet = (TJet*)jetArr->At(i);
      
      if(jet->pt        <= 30)  continue;
      if(fabs(jet->eta) >= 4.7) continue;

      TLorentzVector jetvec;
      jetvec.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
      
      if(jetvec.DeltaR(z1dilep.first.p4) <=0.5) continue;
      if(jetvec.DeltaR(z1dilep.second.p4)<=0.5) continue;
      if(jetvec.DeltaR(z2dilep.first.p4) <=0.5) continue;
      if(jetvec.DeltaR(z2dilep.second.p4)<=0.5) continue;

      if(z1dilep.first.fsrp4.Pt() >0 && jetvec.DeltaR(z1dilep.first.fsrp4) <=0.5) continue;
      if(z1dilep.second.fsrp4.Pt()>0 && jetvec.DeltaR(z1dilep.second.fsrp4)<=0.5) continue;
      if(z2dilep.first.fsrp4.Pt() >0 && jetvec.DeltaR(z2dilep.first.fsrp4) <=0.5) continue;
      if(z2dilep.second.fsrp4.Pt()>0 && jetvec.DeltaR(z2dilep.second.fsrp4)<=0.5) continue;
      
      if(!passJetID(jet)) continue;
      
      njets++;

      if(jetvec.Pt() > jet1.Pt()) {
        jet2 = jet1;
        jet1 = jetvec;
      
      } else if(jetvec.Pt() > jet2.Pt()) {
        jet2 = jetvec;
      }
    }
    double mjj    = (njets>1) ? (jet1 + jet2).M() : -1;
    double detajj = (njets>1) ? fabs(jet1.Eta() - jet2.Eta()) : -1;
    double fishjj = (njets>1) ? (0.18*detajj + (1.92e-4)*mjj) : -1;

    if(zzvec.M() <= 100) continue;
    //if(njets < 1) continue;
    if(njets != 2) continue;
    if(fishjj <= 0.4) continue;

    n4l++;
    if(abs(z1dilep.first.pdgId) == MUON_PDGID) {
      if(abs(z2dilep.first.pdgId) == MUON_PDGID) {
        n4mu++;
      } else {
        n2e2mu++;
      }
    } else {
      if(abs(z2dilep.first.pdgId) == MUON_PDGID) {
        n2e2mu++;
      } else {
        n4e++;
      }
    }
/*    
    double sip11 = (abs(z1dilep.first.pdgId)==ELECTRON_PDGID)  ? ((TElectron*)z1dilep.first.baconObj)->sip3d  : ((TMuon*)z1dilep.first.baconObj)->sip3d;
    double sip12 = (abs(z1dilep.second.pdgId)==ELECTRON_PDGID) ? ((TElectron*)z1dilep.second.baconObj)->sip3d : ((TMuon*)z1dilep.second.baconObj)->sip3d;
    double sip21 = (abs(z2dilep.first.pdgId)==ELECTRON_PDGID)  ? ((TElectron*)z2dilep.first.baconObj)->sip3d  : ((TMuon*)z2dilep.first.baconObj)->sip3d;
    double sip22 = (abs(z2dilep.second.pdgId)==ELECTRON_PDGID) ? ((TElectron*)z2dilep.second.baconObj)->sip3d : ((TMuon*)z2dilep.second.baconObj)->sip3d;    
*/        
    cout << "run= "  << info->runNum;
    cout << " evt= " << info->evtNum;
    cout << " ls= "  << info->lumiSec;
//    cout << fixed << setprecision(2) << " mass4l= " << zzvec.M();
//    cout << fixed << setprecision(2) << " mZ1= "    << z1vec.M();
//    cout << fixed << setprecision(2) << " mZ2= "    << z2vec.M();
    cout << endl;
/*    
    cout << fixed << setprecision(3);
    cout << setw(4) << "" << setw(4) << "Id" << setw(9) << "pt" << setw(9) << "eta" << setw(9) << "phi";
    cout << setw(9) << "px" << setw(9) << "py" << setw(9) << "pz" << setw(9) << "y";
    cout << setw(9) << "iso" << setw(9) << "SIP" << endl;
    
    cout << setw(4) << "ZZ" << setw(4) << 0 << setw(9) << zzvec.Pt() << setw(9) << zzvec.Eta() << setw(9) << zzvec.Phi();
    cout << setw(9) << zzvec.Px() << setw(9) << zzvec.Py() << setw(9) << zzvec.Pz() << setw(9) << zzvec.Rapidity();
    cout << setw(9) << "-" << setw(9) << "-" << endl;
    
    cout << setw(4) << "Z1" << setw(4) << 0 << setw(9) << z1vec.Pt() << setw(9) << z1vec.Eta() << setw(9) << z1vec.Phi();
    cout << setw(9) << z1vec.Px() << setw(9) << z1vec.Py() << setw(9) << z1vec.Pz() << setw(9) << z1vec.Rapidity();
    cout << setw(9) << "-" << setw(9) << "-" << endl;
    
    cout << setw(4) << "Z2" << setw(4) << 0 << setw(9) << z2vec.Pt() << setw(9) << z2vec.Eta() << setw(9) << z2vec.Phi();
    cout << setw(9) << z2vec.Px() << setw(9) << z2vec.Py() << setw(9) << z2vec.Pz() << setw(9) << z2vec.Rapidity();
    cout << setw(9) << "-" << setw(9) << "-" << endl;
    
    cout << setw(4) << "L11" << setw(4) << z1dilep.first.pdgId;
    cout << setw(9) << z1dilep.first.p4.Pt() << setw(9) << z1dilep.first.p4.Eta() << setw(9) << z1dilep.first.p4.Phi();
    cout << setw(9) << z1dilep.first.p4.Px() << setw(9) << z1dilep.first.p4.Py() << setw(9) << z1dilep.first.p4.Pz();
    cout << setw(9) << z1dilep.first.p4.Eta() << setw(9) << z1dilep.first.iso/z1dilep.first.p4.Pt() << setw(9) << fabs(sip11) << endl;
    
    cout << setw(4) << "L12" << setw(4) << z1dilep.second.pdgId;
    cout << setw(9) << z1dilep.second.p4.Pt() << setw(9) << z1dilep.second.p4.Eta() << setw(9) << z1dilep.second.p4.Phi();
    cout << setw(9) << z1dilep.second.p4.Px() << setw(9) << z1dilep.second.p4.Py() << setw(9) << z1dilep.second.p4.Pz();
    cout << setw(9) << z1dilep.second.p4.Eta() << setw(9) << z1dilep.second.iso/z1dilep.second.p4.Pt() << setw(9) << fabs(sip12) << endl;
    
    cout << setw(4) << "L21" << setw(4) << z2dilep.first.pdgId;
    cout << setw(9) << z2dilep.first.p4.Pt() << setw(9) << z2dilep.first.p4.Eta() << setw(9) << z2dilep.first.p4.Phi();
    cout << setw(9) << z2dilep.first.p4.Px() << setw(9) << z2dilep.first.p4.Py() << setw(9) << z2dilep.first.p4.Pz();
    cout << setw(9) << z2dilep.first.p4.Eta() << setw(9) << z2dilep.first.iso/z2dilep.first.p4.Pt() << setw(9) << fabs(sip21) << endl;
    
    cout << setw(4) << "L22" << setw(4) << z2dilep.second.pdgId;
    cout << setw(9) << z2dilep.second.p4.Pt() << setw(9) << z2dilep.second.p4.Eta() << setw(9) << z2dilep.second.p4.Phi();
    cout << setw(9) << z2dilep.second.p4.Px() << setw(9) << z2dilep.second.p4.Py() << setw(9) << z2dilep.second.p4.Pz();
    cout << setw(9) << z2dilep.second.p4.Eta() << setw(9) << z2dilep.second.iso/z2dilep.second.p4.Pt() << setw(9) << fabs(sip22) << endl;

    cout << endl;
*/            
  }
  
  cout << endl;
  cout << "All/4e/4mu/2e2mu : " << n4l << "/" << n4e << "/" << n4mu << "/" << n2e2mu << endl;
  
  cout << endl;
  cout << "Finished with " << intree->GetEntries() << " events processed!" << endl;

}


//--------------------------------------------------------------------------------------------------
bool passMuonID(const TMuon *muon)
{
  // Require to be Global muon OR Tracker muon
  if( !(muon->typeBits & kGlobal) && 
      !((muon->typeBits & kTracker) && (muon->selectorBits & kAllArbitrated)) ) {
    return false;
  }
  
  // impact parameter cuts
  if(fabs(muon->sip3d) >= 100) return false;
  if(fabs(muon->d0)    >= 0.5) return false;
  if(fabs(muon->dz)    >= 1.0) return false;
  
  // Require to be PF muon
  if(!(muon->typeBits & kPFMuon)) return false;
  
  return true;
}


//--------------------------------------------------------------------------------------------------
bool passEleID(const TElectron *ele)
{
  // missing hits cut for conversion rejection
  if(ele->nMissingHits > 1) return false;
  
  // impact parameters cuts
  if(fabs(ele->sip3d) >= 100) return false;
  if(fabs(ele->d0)    >= 0.5) return false;
  if(fabs(ele->dz)    >= 1.0) return false;
  
  int ptBin = (ele->ptHZZ4l > 10) ? 1 : 0;
  int etaBin = -1;
  if     (fabs(ele->scEta) < 0.8)   etaBin = 0;
  else if(fabs(ele->scEta) < 1.479) etaBin = 1;
  else                              etaBin = 2;
  
  if(ptBin == 0 && etaBin == 0) return (ele->mva > ELE_REFERENCE_IDMVA_CUT_BIN0);
  if(ptBin == 0 && etaBin == 1) return (ele->mva > ELE_REFERENCE_IDMVA_CUT_BIN1);
  if(ptBin == 0 && etaBin == 2) return (ele->mva > ELE_REFERENCE_IDMVA_CUT_BIN2);
  if(ptBin == 1 && etaBin == 0) return (ele->mva > ELE_REFERENCE_IDMVA_CUT_BIN3);
  if(ptBin == 1 && etaBin == 1) return (ele->mva > ELE_REFERENCE_IDMVA_CUT_BIN4);
  if(ptBin == 1 && etaBin == 2) return (ele->mva > ELE_REFERENCE_IDMVA_CUT_BIN5);
  
  return false;
}

//--------------------------------------------------------------------------------------------------
bool passJetID(const TJet *jet)
{
  if(jet->neuHadFrac >= 0.99) return false;
  if(jet->neuEmFrac  >= 0.99) return false;
  if(jet->nParticles <= 1)    return false;
  if(fabs(jet->eta)<2.4) {
    if(jet->chHadFrac == 0)    return false;
    if(jet->nCharged  == 0)    return false;
    if(jet->chEmFrac  >= 0.99) return false;
  }
  
  if	 (0    <= fabs(jet->eta) && fabs(jet->eta) < 2.5  && jet->mva < JET_MVA_CUT0) return false;
  else if(2.5  <= fabs(jet->eta) && fabs(jet->eta) < 2.75 && jet->mva < JET_MVA_CUT1) return false;
  else if(2.75 <= fabs(jet->eta) && fabs(jet->eta) < 3	  && jet->mva < JET_MVA_CUT2) return false;
  else if(3    <= fabs(jet->eta) && fabs(jet->eta) < 5	  && jet->mva < JET_MVA_CUT3) return false;
  
  return true;
}

//--------------------------------------------------------------------------------------------------
bool isBetterMuon(const TMuon *mu1, const TMuon *mu2)
{
  // based on function found in UserCode/Mangano/WWAnalysis/AnalysisStep/plugins/MuonCleanerBySegments.cc
  
  if(mu2->typeBits==kStandalone) return true;
  if(mu1->typeBits==kStandalone) return false;
  if((mu1->typeBits & kPFMuon) != (mu2->typeBits & kPFMuon)) return (mu1->typeBits & kPFMuon);
  if((mu1->typeBits & kGlobal) != (mu2->typeBits & kGlobal)) return (mu1->typeBits & kGlobal);
  
  TLorentzVector vMu1; vMu1.SetPtEtaPhiM(mu1->pt, mu1->eta, mu1->phi, MUON_MASS);
  TLorentzVector vMu2; vMu2.SetPtEtaPhiM(mu2->pt, mu2->eta, mu2->phi, MUON_MASS);
  if(mu1->q == mu2->q && vMu1.DeltaR(vMu2)<0.03) {
    return (mu1->ptErr/mu1->pt < mu2->ptErr/mu2->pt);
  } else {
    int nm1 = mu1->nMatchStn;
    int nm2 = mu2->nMatchStn;
    return ( nm1 != nm2 ? nm1 > nm2 : mu1->pt > mu2->pt);
  }
}

//--------------------------------------------------------------------------------------------------
double getEffArea(const double eta) {
  
  if     (fabs(eta) < 1.0)   return 0.19;
  else if(fabs(eta) < 1.479) return 0.25;
  else if(fabs(eta) < 2.0)   return 0.12;
  else if(fabs(eta) < 2.2)   return 0.21;
  else if(fabs(eta) < 2.3)   return 0.27;
  else if(fabs(eta) < 2.4)   return 0.44;
  else                       return 0.52;
}

//--------------------------------------------------------------------------------------------------
double computeIso(const Lepton &lep, const TPhoton *fsrCand, const double rho)
{
  TLorentzVector fsrvec;
  if(fsrCand) fsrvec.SetPtEtaPhiM(fsrCand->pfPt, fsrCand->pfEta, fsrCand->pfPhi, 0);
  
  double fsrCorr = 0;
  double dr = (fsrCand) ? lep.p4.DeltaR(fsrvec) : 0;
  if(abs(lep.pdgId)==ELECTRON_PDGID) {
    const TElectron *ele = (TElectron*)lep.baconObj;
    const double extRadius = 0.4;
    const double intRadius = 0.08;
    
    bool matchEle = (fsrCand && (fsrCand->typeBits & kPFPhoton) && fsrCand->mvaNothingGamma>0.99 && ele->nMissingHits>0 &&
                     fsrCand->scID>-1 && ele->scID>-1 && fsrCand->scID == ele->scID);
     
    if(ele->fiducialBits & kIsEB) {
      if(fsrCand && dr<extRadius && !(matchEle)) fsrCorr = fsrCand->pfPt;
    } else {
      if(fsrCand && dr>=intRadius && dr<extRadius && !(matchEle)) fsrCorr = fsrCand->pfPt;
    }
    return ele->chHadIso04 + TMath::Max(ele->gammaIso04 - fsrCorr + ele->neuHadIso04 - rho*getEffArea(ele->scEta), 0.);
  
  } else if(abs(lep.pdgId)==MUON_PDGID) {
    const TMuon *mu = (TMuon*)lep.baconObj;
    const double extRadius = 0.4;
    const double intRadius = 0;//0.01;
    
    if(fsrCand && dr>=intRadius && dr<extRadius) fsrCorr = fsrCand->pfPt;
    
    return mu->chHadIso04 + TMath::Max(mu->gammaIso04 - fsrCorr + mu->neuHadIso04 - 0.5*(mu->puIso04), 0.);
  
  } else {
    return -1;
  }
}

//--------------------------------------------------------------------------------------------------
const TPhoton* recoverFsr(const Lepton              &lep,
                          const vector<Lepton>      &lepvec,
                          const pair<Lepton,Lepton> &zcand,
			  const TClonesArray	    *photonArr)
{
  vector<TLorentzVector> photonMoms;
  vector<const TPhoton*> photonObjs;
  TLorentzVector phovec;
  
  for(int i=0; i<photonArr->GetEntriesFast(); i++) {
    const TPhoton *photon = (TPhoton*)photonArr->At(i);
    
    phovec.SetPtEtaPhiM(0,0,0,0);
    if((photon->typeBits & kPFMuonPhoton) || (photon->typeBits & kPFPhoton)) {
      // case of muon FSR where photon object is not reconstructed, but photon energy
      // is included in the ECAL energy of the muon object
      phovec.SetPtEtaPhiM(photon->pfPt, photon->pfEta, photon->pfPhi, 0);
    
    } else {
      continue;
    }
    
    if(phovec.Pt()        <= 2)   continue;
    if(fabs(phovec.Eta()) >= 2.4) continue;

    double dR = phovec.DeltaR(lep.p4);
     
    // veto if close to an electron supercluster
    bool flagEleSC = false;
    for(unsigned int j=0; j<lepvec.size(); j++) {
      if(abs(lepvec[j].pdgId)!=ELECTRON_PDGID) continue;
      
      double dPhi  = fabs(phovec.DeltaPhi(lepvec[j].p4));
      double dEta  = fabs(phovec.Eta() - lepvec[j].p4.Eta());
      if((dPhi<2. && dEta<0.05) || sqrt(dPhi*dPhi + dEta*dEta)<0.15) {
        flagEleSC = true;
        break;
      }
    }
    if(flagEleSC) continue;
 
    // check input lepton is the closest lepton to this photon
    bool found_closer_lepton = false;
    for(unsigned int j=0; j<lepvec.size(); j++) {
      if(lep.baconObj == lepvec[j].baconObj) continue;
      
      double tmp_dR = phovec.DeltaR(lepvec[j].p4);
      if(tmp_dR < dR) {
        found_closer_lepton = true;
        break;
      }
    }
    if(found_closer_lepton) continue;  

    // Z mass OK?
    double oldMass = (zcand.first.p4 + zcand.second.p4).M();
    double newMass = (zcand.first.p4 + zcand.second.p4 + phovec).M();
    if( newMass <= 4.   || 
        newMass >= 100. ||
        fabs(newMass - Z_MASS) >= fabs(oldMass - Z_MASS) )
      continue;
    
    // "keep all photons close to one of the 4L leptons..."
    bool use = false;
    if(dR < 0.07) {
      use = true;
    } else if(dR<0.5 && phovec.Pt()>4 && photon->isoForFsr03<1.0*phovec.Pt()) {  // "need tighter cuts for other photons..."
      use = true;
    }
    if(use) {
      photonMoms.push_back(phovec);
      photonObjs.push_back(photon);
    }
  }

  // choose the best one
  double maxPt  =  0, minDR  = 999;
  int    iMaxPt = -1, iMinDR = -1;
  for(unsigned int i=0; i<photonMoms.size(); i++) {
    if(photonMoms[i].Pt() > maxPt) {
      maxPt  = photonMoms[i].Pt();
      iMaxPt = i;
    }
    double tmp_dR = photonMoms[i].DeltaR(lep.p4);
    if(tmp_dR < minDR) {
      minDR  = tmp_dR;
      iMinDR = i;
    }
  }
 
  if(maxPt > 4) {
    return photonObjs[iMaxPt];
  } else if(minDR < 999) {
    return photonObjs[iMinDR];
  } else {
    return 0;
  }
}
