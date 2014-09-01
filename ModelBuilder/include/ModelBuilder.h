#ifndef MODELBUILDER
#define MODELBUILDER

#include "include/diagonalize.h"
#include "TDirectory.h" 
#include "TFile.h" 
#include "TEventList.h" 
#include "TH2F.h"
#include "RooAbsReal.h"
#include <cstring>

using namespace std;

class ModelBuilder {
  public:

  ModelBuilder( int id ){
  	sel_i=id;
  };
  ~ModelBuilder(){};

  // functions building the model
  const char * doubleexp(RooWorkspace *,RooRealVar &,std::string ext="") ;
   
  void buildModel(TFile *, TDirectory *);

  // These are the defitions of one model (category)
  double luminosity;
  int 	 run_correction;
  int    run_bkgsys;
  int    runphotonmodel;
  std::string cutstring ;//"jet1mprune > 60";

  std::string varstring;
  double xmin;
  double xmax;
  TH1F *lMet;
  
  private:
  // For now, these are pretty inflexible but we could define control samples properly also
  // very hard coded fitting 
  void buildAndFitModels(TDirectory *,RooRealVar &,std::string proc="Zvv") ;
  // very hard coded control regions
  void makeAndImportDataSets(TFile *, RooRealVar & );
  void createPhotonModel(TDirectory *);
  void createWeightedPhotonDataset(TDirectory *);
  // Constants defined 
  //const std::string outfilename="workspace.root";
  //const double brscaleFactorZvv=5.942; 
  //const double brscaleFactorWlv=1.017; 
  //double luminosity(19.7);
  //
  int sel_i;
  RooWorkspace *wspace;
};
#endif
