#include "BaconProd/Utils/interface/QGLikelihoodCalculator.hh"
#include "TMVA/Reader.h"
#include <vector>
#include <iostream>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
QGLikelihoodCalculator::QGLikelihoodCalculator():
  fIsInitialized(false),
  fReader(0),
  fMethodTag("")
{}

//--------------------------------------------------------------------------------------------------
QGLikelihoodCalculator::~QGLikelihoodCalculator()
{
  delete fReader;
}

//--------------------------------------------------------------------------------------------------
void QGLikelihoodCalculator::initialize(const std::string methodTag, const std::string weightFile)
{
  fReader = new TMVA::Reader();
  
  fReader->AddVariable("nvtx",       &_nvtx);
  fReader->AddVariable("jetPt",      &_jetPt);
  fReader->AddVariable("jetEta",     &_jetEta);
  fReader->AddVariable("jetPhi",     &_jetPhi);
  fReader->AddVariable("beta",       &_beta);
  fReader->AddVariable("betaStar",   &_betaStar);
  fReader->AddVariable("nParticles", &_nParticles);
  fReader->AddVariable("nCharged",   &_nCharged);
  fReader->AddVariable("dRMean",     &_dRMean);
  fReader->AddVariable("ptD",        &_ptD);
  fReader->AddVariable("frac01",     &_frac01);
  fReader->AddVariable("frac02",     &_frac02);
  fReader->AddVariable("frac03",     &_frac03);
  fReader->AddVariable("frac04",     &_frac04);
  fReader->AddVariable("frac05",     &_frac05);
  
  fMethodTag = methodTag;
  fReader->BookMVA(fMethodTag, weightFile);
  
  fIsInitialized = true;
}

//--------------------------------------------------------------------------------------------------
void QGLikelihoodCalculator::mvaValues(float *qgvals, const float nvtx,
                                       const float jetPt, const float jetEta, const float jetPhi, const float beta, const float betaStar, 
				       const float nParticles, const float nCharged, const float dRMean, const float ptD,
		                       const float frac01, const float frac02, const float frac03, const float frac04, const float frac05,
				       const bool printDebug)
{
  _nvtx       = nvtx;
  _jetPt      = jetPt;
  _jetEta     = jetEta;
  _jetPhi     = jetPhi;
  _beta       = beta;
  _betaStar   = betaStar;
  _nParticles = nParticles;
  _nCharged   = nCharged;
  _dRMean     = dRMean;
  _ptD        = ptD;
  _frac01     = frac01;
  _frac02     = frac02;
  _frac03     = frac03;
  _frac04     = frac04;
  _frac05     = frac05;
  
  const std::vector<float> &vals = fReader->EvaluateMulticlass(fMethodTag);
  
  if(printDebug) {
    std::cout << "[QGLikelihoodCalculator]" << std::endl;
    std::cout << "Method Tag: " << fMethodTag << std::endl;
    std::cout << "Inputs: jetPt= " << _jetPt << "  jetEta= " << _jetEta << "  jetPhi= " << _jetPhi;
    std::cout << "  beta= " << _beta << "  betaStar= " << _betaStar;
    std::cout << "  nParticles= " << _nParticles << "  nCharged= " << _nCharged;
    std::cout << "  dRMean= " << _dRMean << "  ptD= " << _ptD;
    std::cout << "  frac01= " << _frac01 << "  frac02= " << _frac02 << "  frac03= " << _frac03 << "  frac04= " << _frac04 << "  frac05= " << _frac05;
    std::cout << std::endl;
    std::cout << " > Quark LL   = " << vals[0] << std::endl;
    std::cout << " > Gluon LL   = " << vals[1] << std::endl;
    std::cout << " > Pile-up LL = " << vals[2] << std::endl;
  }
  
  qgvals[0] = vals[0];
  qgvals[1] = vals[1];
  qgvals[2] = vals[2];
}
