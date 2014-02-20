#ifndef BACONPROD_UTILS_QGLIKELIHOODCALCULATOR_HH
#define BACONPROD_UTILS_QGLIKELIHOODCALCULATOR_HH

#include <string>

// forward class declarations
namespace TMVA {
  class Reader;
}

namespace baconhep {

  class QGLikelihoodCalculator
  {
    public:
      QGLikelihoodCalculator();
      ~QGLikelihoodCalculator();
      
      void initialize(const std::string methodTag, const std::string weightFile);
      
      bool isInitialized() const {return fIsInitialized;}
      
      void mvaValues(float *qgvals,
                     const float nvtx,
                     const float jetPt, const float jetEta, const float jetPhi,
		     const float beta, const float betaStar,
		     const float nParticles, const float nCharged,
		     const float dRMean, const float ptD,
		     const float frac01, const float frac02, const float frac03, const float frac04, const float frac05,
		     const bool printDebug=false);
      
    
    private:
      bool fIsInitialized;
      
      TMVA::Reader *fReader;
      std::string fMethodTag;
      
      // input variables to compute the likelihood
      float _nvtx,
            _jetPt, _jetEta, _jetPhi,
            _beta, _betaStar,
            _nParticles, _nCharged,
            _dRMean, _ptD,
            _frac01, _frac02, _frac03, _frac04, _frac05;    
  };
}
#endif
