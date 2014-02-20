#ifndef BACONPROD_NTUPLER_FILLERPF_HH
#define BACONPROD_NTUPLER_FILLERPF_HH

#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TClonesArray;


namespace baconhep
{
  class FillerPF
  {
    public:
      FillerPF();
      ~FillerPF();
      
       void fill(TClonesArray       *array,    // output array to be filled
		 TClonesArray       *iVtxCol,
		 const edm::Event   &iEvent);  // event info

      
      // EDM object collection names
      std::string fPFName;
      std::string fPVName;

  };
}
#endif
