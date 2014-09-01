
#include "fewz_weights.C"
#include "TH2F.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TEventList.h"
#include "TTree.h"
#include "TMath.h"
#include "TROOT.h"
#include <vector>
#include <iomanip>

#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/PxPyPzE4D.h"
using namespace std;

double globalXSECSUM=0;
double globalXSECSUMMG=0;

double deltaphi(double phi1, double phi2){
        // do somthing to sort out delta R??
        double mpi = TMath::Pi();
        double res = phi1 - phi2;
        while (res > mpi) res-= 2*mpi;
        while (res <= -mpi) res+= 2*mpi;
        return TMath::Abs(res);
}


TH2F createMadgraphPowhegReweight(TTree *tree, TH2F *target, double mLow, double mHigh){

  //const int nptbins = 21;
  //const int nrapbins = 8;
  //
  //double rap_bin[nrapbins] = {0.,0.2,0.4,0.7,1.1,1.9,2.4,1000.0};
  //double pt_bin[nptbins] = {0.0,20.0,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,100.,120.,150.,300.,600.,1000.};

  std::cout << "Generate weights for mass in range " << mLow << " -> " << mHigh << std::endl;
  TH2F *template_base = (TH2F*) target->Clone();
  template_base->SetName(Form("template_base_%s",target->GetName()));
  // clear !
  for (int i=1;i<template_base->GetNbinsX()+1;i++){
    for (int j=1;j<template_base->GetNbinsY()+1;j++){
      template_base->SetBinContent(i,j,0);
    }
  }
  std::cout << "Draw command, can it be so slow?" << std::endl;
//  tree->Draw(Form("TMath::Abs(genZ.Rapidity()):genZ.Pt()>>template_base_%s",target->GetName()),Form("weight*(genZ.M()>%g && genZ.M()<%g)",mLow,mHigh));  // weight here is /pb so get x-section in pb
  TH2F *rethist = (TH2F*) target->Clone();
  // The histograms for POWHEG are absolute cross-section in each bin (not differential/GeV etc) judging by the shape of the spectrum so no need to scale by bin volume
  rethist->Divide(template_base);
  rethist->SetName(Form("madgraph_to_powheg_%s",target->GetName()));
  std::cout << ".... madgraph x-sec = " << template_base->Integral()   << ", powheg x-sec = " << target->Integral()  << std::endl; 
  globalXSECSUM+= target->Integral() ;
  globalXSECSUMMG+= template_base->Integral() ;
  return *rethist;
}

void fillZvvHistogram(TH2F *hist, TH1F* hist_1d, TTree *tree, std::vector<TH2F> &h_mtop, int *masses, int nm, bool doweight){

  std::string cutstr = "mvamet>200 && jet1CHF>0.2 && jet1NHF<0.7 && jet1NEMF<0.7 && jet1.Pt()>150 && TMath::Abs(jet1.Eta())<2.0 && (metFiltersWord==1023||metFiltersWord==511) && ((trigger&1)==1 || (trigger&2)==2) && njets<3 && ( ntaus == 0 && nphotons == 0 && nlep == 0)";

  float weight_in, mvametPhi, weight_pu, mvamet;

  // remember this weight in isn't scaled to lumi so do it !
  double lumi = 19700.;

  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > *jet  = new ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >();
  TBranch *brjet ;
  brjet = tree->GetBranch("jet1");
  brjet->SetAddress(&jet);
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > *jet2  = new ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >();
  TBranch *brjet2 ;
  brjet2 = tree->GetBranch("jet2");
  brjet2->SetAddress(&jet2);
  
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > *genZ  = new ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >();
  TBranch *brgenZ ;
  brgenZ = tree->GetBranch("genZ");
  brgenZ->SetAddress(&genZ);

  std::cout << "Drawing cutlist" << std::endl;
  // simple selection for the jet  
  tree->Draw(">>cutlist",cutstr.c_str());
  TEventList *keep_points = (TEventList*)gDirectory->Get("cutlist");

  tree->SetBranchAddress("weight",&weight_in);
  tree->SetBranchAddress("mvamet",&mvamet);
  tree->SetBranchAddress("puWeight",&weight_pu);
  tree->SetBranchAddress("mvametPhi",&mvametPhi);
  // Cut M>20!
  std::cout << "Begin loop in Filling Zvv hist" << std::endl;

  double weight=1.;
  for (int ev=0;ev<tree->GetEntries();ev++){
    tree->GetEntry(ev);
    if ( !(keep_points->Contains(ev)) ) continue;

    // Additional cuts which we cannot do with the string
    if (jet2->Pt() > 30){
      if  (deltaphi(jet->Phi(),jet2->Phi()) >= 2) continue;
    }
    if ( deltaphi(jet->Phi(),mvametPhi)<2 ) continue; 

    double M = genZ->M();
    if (M<20) continue;
    //if (genZ->Pt() < 100) continue; 
    double pt = jet->Pt(); 
    double y  = TMath::Abs(jet->Eta()); 
    // Cut on jet pt and y (since photon sample effectively has this applied)
    //if (pt< 150 || y > 2)  continue; 

    weight = weight_in*lumi*weight_pu;

    if (doweight){
     int thindex = -1;
     for (int mi=0; mi<nm-1;mi++){
      double ml = (double) masses[mi];
      double mh = (double) masses[mi+1];
      if ( M>ml && M<mh) { 
		thindex=mi;
		break;
      }
     }
     if (thindex<0 )continue;
     double zpt = genZ->Pt();
     double zy  = TMath::Abs(genZ->Rapidity());
     int b_mtop = h_mtop[thindex].FindBin(zpt,zy);
     double w_mtop=h_mtop[thindex].GetBinContent(b_mtop);
     weight*=w_mtop*weight_FEWZ8TeV( zpt, zy, M);   
    }

    hist->Fill(pt,y,weight);
    hist_1d->Fill(mvamet,weight);
  }
}

void reweight_photonjet() {

  // Generator level reweighting derivation to be applied to Photon + Jet Data
  gROOT->SetBatch(1);
  gROOT->ProcessLine(".L fewz_weights.C"); 

  // Where the Photon+Jet lives ...
  std::string PhotonSample = "Output_Photon.root";

  // Where the Zvv lives 
  std::string ZvvSample    = "root://eoscms//eos/cms/store/cmst3/user/pharris/Marco/s12-zjets-ptz100-v7a.root";
  TFile *f_ZvvSample = TFile::Open(ZvvSample.c_str());
  TTree *tree_Zvv = (TTree*) f_ZvvSample->Get("tree");

  // Where the Zvv lives ! post selection !
 // std::string ZvvSample_selection    = "Output_Marco.root";
//  TFile *f_ZvvSample_selection = TFile::Open(ZvvSample.c_str());
//  TTree *tree_Zvv_selection = (TTree*) f_ZvvSample->Get("DY");

  TFile *f_PhoSample = TFile::Open(PhotonSample.c_str());
  TTree *tree_Pho = (TTree*) f_PhoSample->Get("GJet");

  std::cout << "Generating weights MAGRAPH->POWHEG" << std::endl;
  // Step 1 is to correct the Zvv sample using NNL0 corrections. 
  // We started with Madgraph Zvv sample BEFORE cuts so ...
  // 1. Fill several histograms (one per mass bin of PT-Y), use "weight" to scale will give differential cross-section for each bin
  // 2. Provided are reweights to go from POWHEG->FEWZ as arrays. Also have histograms which give MADGRAPH->POWHEG!  
  
  TFile *f_powheg_xsection = TFile::Open("/afs/cern.ch/work/n/nckw/private/monojet/DYM_40mass8rap.root");
  
  const int nmbins = 41;
  int  massedges[nmbins] = {15, 20 , 25 , 30 , 35 , 40 , 45 , 50 , 55 , 60 , 64 , 68 , 72 , 76 , 81 , 86 , 91 , 96 , 101 , 106 , 110 , 115 , 120 , 126 , 133 , 141 , 150 , 160 , 171 , 185 , 200 , 220 , 243 , 273 , 320 , 380 , 440 , 510 , 600 , 1000 , 1500};

  std::vector<TH2F> th2f_v_madgraph_to_powheg;

  for (int mi=0; mi<nmbins-1;mi++){
    double ml = (double) massedges[mi];
    double mh = (double) massedges[mi+1];

    TH2F *powheg = (TH2F*) f_powheg_xsection->Get(Form("xsec_%d%dbin",massedges[mi],massedges[mi+1]));   
    th2f_v_madgraph_to_powheg.push_back(createMadgraphPowhegReweight(tree_Zvv,powheg,ml,mh));	
  } 

  std::cout << "Total Madgraph/Powheg Z cross-section = " << globalXSECSUMMG << "/" << globalXSECSUM << std::endl;
  std::cout << "Reweigting Zvv MADGRAPH->POWHEG->FEWZ" << std::endl;
  // Step 2 is to produce PT-Y ratio AFTER kinematic gen cuts for photon/corrected_Zvv
  

  // 2D reweighting histogram will be PT-Y between 0->1000 and, 0->10
  const int nptbins = 28;
  const int nrapbins = 11;
  double pt_bin[nptbins] = {100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,450,500,550,600,800,1000,1200};
  double rap_bin[nrapbins] = {0.,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0};
 
  TH2F * th2f_zvv_hist_preweighting = new TH2F("Zvv_madgraph","Zvv_madgraph",nptbins-1,pt_bin,nrapbins-1,rap_bin);
  TH2F * th2f_zvv_hist_weighting    = new TH2F("Zvv_fewz","Zvv_fewz",nptbins-1,pt_bin,nrapbins-1,rap_bin);

  TH1F * th1f_mvamethist_preweighting = new TH1F("Zvv_madgraph_met","Zvv_madgraph_met",50,200,1000);
  TH1F * th1f_mvamethist_weighting    = new TH1F("Zvv_fewz_met","Zvv_fewz_met",50,200,1000);

  fillZvvHistogram(th2f_zvv_hist_preweighting,th1f_mvamethist_preweighting,tree_Zvv,th2f_v_madgraph_to_powheg,massedges,nmbins-1,0);
  fillZvvHistogram(th2f_zvv_hist_weighting,th1f_mvamethist_weighting,tree_Zvv,th2f_v_madgraph_to_powheg,massedges,nmbins-1,1);
  std::cout << "Produced Weighed Zvv generator template" << std::endl;
  
  // repeat for photons
  TH2F * th2f_hist_photons    = new TH2F("Pho_mg","Pho_mg",nptbins-1,pt_bin,nrapbins-1,rap_bin);
  tree_Pho->Draw("TMath::Abs(jeta):jpt>>Pho_mg","weight*( jpt>150 && TMath::Abs(jeta)<2 && jmdphi>2 && jdphi < 2 && mvamet>200 && pho_idmva>0.1)");  // most cuts already appliued  (in fact also these ones are applied?!
  //th2f_hist_photons->Scale(1./19700);   // Account for the fact that luminosity is in the weights, Zvv is already accounting for this part!
  TH1F * th1f_hist_photons_mvamet    = new TH1F("Pho_mg_met","Pho_mg_met",50,200,1000);
  tree_Pho->Draw("mvamet>>Pho_mg_met","weight*( jpt>150 && TMath::Abs(jeta)<2 && jmdphi>2 && jdphi < 2 && mvamet>200 && pho_idmva>0.1)");  // most cuts already appliued  (in fact also these ones are applied?!

  std::cout << "Event Yields --- " << std::endl;
  std::cout << "   Zvv (m-graph)" << th2f_zvv_hist_preweighting->Integral() << std::endl;
  std::cout << "   Zvv (fewz)   " << th2f_zvv_hist_weighting->Integral() << std::endl;
  std::cout << "   Photon       " << th2f_hist_photons->Integral() << std::endl;

  // Finally make weights photon->Zvv /******************** WARNING, FOR NOW I USED PRE FEWZ-WEIGTED SAMPLES!!!!!!!! **********/
  TH2F * hist_photon_reweighting = (TH2F*)th2f_zvv_hist_preweighting->Clone();   // acceptance effect is there but no photon-efficiency correction (assume its high wrt trigger?)
  hist_photon_reweighting->SetName("Photon_to_Zvv_weights");
  hist_photon_reweighting->Divide(th2f_hist_photons);

  TH2F * hist_photon_reweighting_fewz = (TH2F*)th2f_zvv_hist_weighting->Clone();   // acceptance effect is there but no photon-efficiency correction (assume its high wrt trigger?)
  hist_photon_reweighting_fewz->SetName("Photon_to_Zvv_FEWZ_weights");
  hist_photon_reweighting_fewz->Divide(th2f_hist_photons);

  TH1F * hist_photon_reweighting_met = (TH1F*)th1f_mvamethist_preweighting->Clone();   // acceptance effect is there but no photon-efficiency correction (assume its high wrt trigger?)
  hist_photon_reweighting_met->SetName("Photon_to_Zvv_weights_met");
  hist_photon_reweighting_met->Divide(th1f_hist_photons_mvamet);

  TH1F * hist_photon_reweighting_fewz_met = (TH1F*)th1f_mvamethist_weighting->Clone();   // acceptance effect is there but no photon-efficiency correction (assume its high wrt trigger?)
  hist_photon_reweighting_fewz_met->SetName("Photon_to_Zvv_FEWZ_weights_met");
  hist_photon_reweighting_fewz_met->Divide(th1f_hist_photons_mvamet);

  // Save the work
  TFile *fout = new TFile("weight_pho_to_zvv.root","RECREATE");
  th2f_zvv_hist_preweighting->Write();
  th2f_zvv_hist_weighting->Write();
  th2f_hist_photons->Write();
  hist_photon_reweighting->Write();
  hist_photon_reweighting_fewz->Write();
  hist_photon_reweighting_met->Write();
  hist_photon_reweighting_fewz_met->Write();
  TDirectory *dir = (TDirectory*) fout->mkdir("magraph_to_powheg_weights");
  dir->cd();
  for (std::vector<TH2F>::iterator it = th2f_v_madgraph_to_powheg.begin(); it!= th2f_v_madgraph_to_powheg.end();it++){
    (*it).Write();
  }
  std::cout << "Saved reweighting histograms in " << fout->GetName() << std::endl;
  fout->Close();
}
