//--------------------------------------------------------------------------------------------------
// $Id $
//
// Tau Iso MVA
//
// Authors: M. Chan (implemented by P. Harris) 
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_TAUISOMVA_H
#define MITPHYSICS_UTILS_TAUISOMVA_H

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include <vector>
#include "MitAna/DataTree/interface/PFTau.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"

class TRandom3;
class GBRForest;

namespace TMVA {
  class Reader;
}

namespace mithep {
  class TauIsoMVA {
  public:
    struct IsoRings
    {
      std::vector<int> niso;
      std::vector<std::vector<float> > rings;
      std::vector<std::vector<float> > shapes;
      
      std::vector<float> getVector()
      {
	std::vector<float> all;
	all.reserve(33);
      
	for(uint i = 0; i < niso.size(); i++)
	  all.push_back(niso[i]);
	
	for(uint i = 0; i < rings.size(); i++)
	  all.insert(all.end(), rings[i].begin(), rings[i].end());
      
	for(uint i = 0; i < shapes.size(); i++)
	  all.insert(all.end(), shapes[i].begin(), shapes[i].end());
	
	return all;
      }
    };
    TauIsoMVA();
    ~TauIsoMVA(); 
    void     Initialize     (TString iWeightFile);
    void     InitializeGBR  (TString iWeightFile);
    double   MVAValue       (const PFTau *iTau,double iRho);
    IsoRings computeIsoRings(const PFTau *iTau);
    
  protected:      
    bool                     fGBR;
    GBRForest               *fGBRReader;
    TMVA::Reader            *fReader;
    ClassDef(TauIsoMVA,1)
      };
}
#endif
