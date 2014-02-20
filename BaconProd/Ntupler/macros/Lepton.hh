#ifndef BACONPROD_NTUPLER_LEPTON_HH
#define BACONPROD_NTUPLER_LEPTON_HH
#include <TObject.h>
#include <TLorentzVector.h>

class Lepton
{
  public:
    Lepton():baconObj(0),pdgId(0),q(0),iso(-1) {
      p4.SetPtEtaPhiE(0,0,0,0);
      fsrp4.SetPtEtaPhiE(0,0,0,0);
    }
    Lepton(const TObject *lepton, const int id):iso(-1) { 
      baconObj = lepton;
      pdgId    = id;
      q        = (id>0) ? -1 : 1;
      p4.SetPtEtaPhiE(0,0,0,0);
      fsrp4.SetPtEtaPhiE(0,0,0,0);
    }
    ~Lepton(){}
    const TObject *baconObj;
    int pdgId;
    int q;
    double iso;
    TLorentzVector p4;
    TLorentzVector fsrp4;
};
#endif
