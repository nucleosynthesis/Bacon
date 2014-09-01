#ifndef DIAGONALIZE
#define DIAGONALIZE

// ROOT includes
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TF1.h"
#include "TAxis.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TIterator.h"

// RooFit includes
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooGlobalFunc.h"
#include "RooCmdArg.h"

// RooStats includes
#include "RooWorkspace.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>

using namespace std;

void freezeParameters(RooArgSet *args, bool freeze=1){
  TIterator* iter(args->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
      rrv->setConstant(freeze);
  }
}
TH1F *generateTemplate(TH1F *base, TTree *tree, std::string varname, std::string weightname, std::string cut="1>0"){

  assert(tree);
  assert(base);
  // Correction point to point
  TH1F *histNew = (TH1F*) base->Clone("hnew");
  histNew->Sumw2();
  histNew->SetName(Form("th1f_%s",tree->GetName()));
  if (weightname!="" ) tree->Draw(Form("%s >> %s",varname.c_str(),histNew->GetName()),Form("%s*(%s)",weightname.c_str(),cut.c_str()));
  else tree->Draw(Form("%s >> %s",varname.c_str(),histNew->GetName()),cut.c_str());
  histNew->GetXaxis()->SetTitle(varname.c_str());
  return histNew;
}

TH1F *generateTemplate(TH1F *base, RooFormulaVar *pdf_num, RooRealVar &var, RooDataSet *data, int dir=1, double weightsf=1.){
  //assert(pdf_num);
  assert(base);
  // Correction point to point
  TH1F *histNew = (TH1F*) base->Clone("hnew");
  histNew->Sumw2();
  if (pdf_num) histNew->SetName(Form("th1f_corrected_%s",data->GetName()));
  else histNew->SetName(Form("th1f_uncorrected_%s",data->GetName()));
  int nevents = data->numEntries();
  const char *varname = var.GetName();
  for (int ev=0;ev<nevents;ev++){
    const RooArgSet *vw = data->get(ev);
    double val    = vw->getRealValue(varname);
    double weight = data->weight();
    //std::cout << val << ", " << weight <<std::endl;
    var.setVal(val);
    double cweight = weight*weightsf;
    if (pdf_num) { 
      if (dir>=0){
    	cweight *= pdf_num->getVal(var);
//	std::cout << "EVT - " << ev << " Corrected " <<var.getVal() << ", orig = " << weight*weightsf << ", new = " << cweight << std::endl;  
      } else {
    	cweight /= pdf_num->getVal(var);
      }
    }
    histNew->Fill(val,cweight);
  }
  histNew->GetXaxis()->SetTitle(varname);
  return histNew;
}

/*
TH1F *generateTemplate(TH1F *base, RooFormulaVar *pdf_num, RooRealVar *var, int dir=1){
  // Correction based on Bin-centers
  int nbins = base->GetNbinsX();
  TH1F *histNew = (TH1F*) base->Clone("hnew");
  histNew->SetName(base->GetName());
  for (int bi=0;bi<nbins;bi++){
    var->setVal(base->GetBinCenter(bi+1));
    double bc=0;
    if (dir>=0){
      //std::cout << var->getVal() << ", " << pdf_num->getVal() <<std::endl;
      bc = pdf_num->getVal()*base->GetBinContent(bi+1);
    } else {
      bc = base->GetBinContent(bi+1)/pdf_num->getVal();
    }

    histNew->SetBinContent(bi+1,bc);
  }
  
  return histNew;
}
*/
void getArgSetParameters(RooArgList & params,std::vector<double> &val){

  val.clear();
  TIterator* iter(params.createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
      //std::cout << "RooContainer:: -- Getting values from parameter -- " << std::endl;
      //rrv->Print();
      val.push_back(rrv->getVal());
  }
}

void setArgSetParameters(RooArgList & params,std::vector<double> &val){

  int i = 0;
  TIterator* iter(params.createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
      // each var must be shifted
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
      //std::cout << "RooContainer:: -- Setting values from parameter -- " << std::endl;
      //rrv->Print();
      std::cout << " Set parameter " << rrv->GetName() << " to " << val[i] << std::endl;
      rrv->setVal(val[i]);
      i++;
  }
}

TH1F *fillHistogramFromPdf(TH1F *base, RooAbsPdf *pdf, RooRealVar *x){
  TH1F *newHist = (TH1F*)base->Clone();
  // clear
  for (int b=1;b<=newHist->GetNbinsX();b++){
    newHist->SetBinContent(b,0);
  }
  const RooArgList args(*x);
  //pdf->fillHistogram(newHist,args,1.,0,true);
  for (int b=1;b<=base->GetNbinsX();b++){
     x->setRange(Form("bin_%d",b),base->GetBinLowEdge(b),base->GetBinLowEdge(b+1));
     RooAbsReal *inte = pdf->createIntegral(RooArgSet(*x),RooArgSet(*x),Form("bin_%d",b));
     newHist->SetBinContent(b,inte->getVal());
     newHist->SetBinError(b,0);
  }

  return newHist;
}

void generateVariations(TH1F *base_hist,RooFitResult *res_ptr, RooAbsPdf *pdf_ptr, RooRealVar *x,std::vector<TH1F> & v_th1f_,RooWorkspace *wspace){
// Generate variations from a simple getHistogram from a Pdf
  TMatrixD cov = res_ptr->covarianceMatrix();
  TVectorD eval;
  TMatrixD evec   = cov.EigenVectors(eval);
  evec.Print();
  std::cout << "Eigenvalues (squares of errors)"<< std::endl;
  eval.Print();
  // ---------------------------------------------------------------------------------------------------------------------

  // now we know that the eigenvalues of the covariance matrix must scale each parameter as given by the new paraemeters
  int n_par = eval.GetNoElements();
  std::cout << "Number of Parameters from the Fit -- " << n_par  << std::endl;
  RooArgList rooFloatParameters = res_ptr->floatParsFinal();  // Why not return a set ????!!!
  // Sigh 
  RooArgList rooParameters;
  TIterator *parsit = rooFloatParameters.createIterator();
  while (RooAbsArg *arg = (RooAbsArg*)parsit->Next()) {
     RooRealVar *vit = wspace->var(arg->GetName());
     rooParameters.add(*vit);
  }

  //RooArgSet rooParameters;
  //TIterator *parsit = rooParameterList->createIterator
  pdf_ptr->getParameters(*x)->Print();

  // fill vector with original parameters
  std::vector <double> original_values;
  getArgSetParameters(rooParameters,original_values);

  TH1F * temp_hist;

  // now loop over columns of the eigenvector Matrix
  std::vector<double> new_values;
  for (int par=0;par<n_par;par++){
    
    // this row in evec is the scales for the parameters

    double err = TMath::Sqrt(eval[par]);

    new_values.clear(); // make sure its empty before we fill it
    for (int i=0;i<n_par;i++){
      new_values.push_back(original_values[i] + evec[i][par]*err);	
    }

    std::cout << "Systematic from parameter "<< par << " +1 sigma" << std::endl;
    setArgSetParameters(rooParameters,new_values);

    temp_hist = (TH1F*) fillHistogramFromPdf(base_hist,pdf_ptr,x);
    temp_hist->Scale(base_hist->Integral()/temp_hist->Integral());
    temp_hist->SetName(Form("th1f_%s_%s_param%dUp",pdf_ptr->GetName(),res_ptr->GetName(),par));
    v_th1f_.push_back(*temp_hist);

    // now -1 sigma
    new_values.clear();  // need to clear again
    for (int i=0;i<n_par;i++){
        new_values.push_back(original_values[i] - evec[i][par]*err);	
    }

    std::cout << "Systematic from parameter "<< par << " -1 sigma" << std::endl;
    setArgSetParameters(rooParameters,new_values);

    temp_hist = (TH1F*) fillHistogramFromPdf(base_hist,pdf_ptr,x);
    temp_hist->Scale(base_hist->Integral()/temp_hist->Integral());
    temp_hist->SetName(Form("th1f_%s_%s_param%dDown",pdf_ptr->GetName(),res_ptr->GetName(),par));
    v_th1f_.push_back(*temp_hist);
 
    // after each parameter we reset the originals back
    std::cout << "Reset to Original Values " << std::endl;
    setArgSetParameters(rooParameters,original_values);
  } 
  std::cout << "Generated systematic varitions from " << res_ptr->GetName()<<std::endl;    
}

void generateVariations(TH1F *base_hist,RooFitResult *res_ptr, RooFormulaVar *pdf_ptr, RooRealVar *x,std::vector<TH1F> & v_th1f_,RooWorkspace *wspace, std::string dataSetName){
// Generate variations for a CORRECTION to a dataset  
  RooDataSet *data = (RooDataSet*)wspace->data(dataSetName.c_str());
  assert(data);
  // The "Pdf" will be a correction function and we will correct the Z->nunu/W->lnu histograms with it
  //
  // The ideas
  // 1. "Uncorrect" the histogram by inverting the correciton. 
  // 2. "Make +/- 1sigma variations from each of the parameters in the correction (r) then "recorrect" to produce new varied histos
  // 3. addtionally what if the correction is not just the same in control/signal -> can we make a systematic on this?
 
  TMatrixD cov = res_ptr->covarianceMatrix();
  TVectorD eval;
  TMatrixD evec   = cov.EigenVectors(eval);
  evec.Print();
  std::cout << "Eigenvalues (squares of errors)"<< std::endl;
  eval.Print();
  // ---------------------------------------------------------------------------------------------------------------------

  // now we know that the eigenvalues of the covariance matrix must scale each parameter as given by the new paraemeters
  int n_par = eval.GetNoElements();
  std::cout << "Number of Parameters from the Fit -- " << n_par  << std::endl;
  RooArgList rooFloatParameters = res_ptr->floatParsFinal();  // Why not return a set ????!!!
  // Sigh 
  RooArgList rooParameters;
  TIterator *parsit = rooFloatParameters.createIterator();
  while (RooAbsArg *arg = (RooAbsArg*)parsit->Next()) {
     RooRealVar *vit = wspace->var(arg->GetName());
     rooParameters.add(*vit);
  }

  //RooArgSet rooParameters;
  //TIterator *parsit = rooParameterList->createIterator
  pdf_ptr->getParameters(*x)->Print();

  // fill vector with original parameters
  std::vector <double> original_values;
  getArgSetParameters(rooParameters,original_values);

  TH1F * temp_hist;

  // now loop over columns of the eigenvector Matrix
  std::vector<double> new_values;
  for (int par=0;par<n_par;par++){
    
    // this row in evec is the scales for the parameters

    double err = TMath::Sqrt(eval[par]);

    new_values.clear(); // make sure its empty before we fill it
    for (int i=0;i<n_par;i++){
      new_values.push_back(original_values[i] + evec[i][par]*err);	
    }

    std::cout << "Systematic from parameter "<< par << " +1 sigma" << std::endl;
    setArgSetParameters(rooParameters,new_values);

    temp_hist = (TH1F*) generateTemplate(base_hist,pdf_ptr,*x,data);
    temp_hist->SetName(Form("th1f_%s_%s_param%dUp",dataSetName.c_str(),res_ptr->GetName(),par));
    v_th1f_.push_back(*temp_hist);

    // now -1 sigma
    new_values.clear();  // need to clear again
    for (int i=0;i<n_par;i++){
        new_values.push_back(original_values[i] - evec[i][par]*err);	
    }

    std::cout << "Systematic from parameter "<< par << " -1 sigma" << std::endl;
    setArgSetParameters(rooParameters,new_values);


    temp_hist = (TH1F*) generateTemplate(base_hist,pdf_ptr,*x,data);
    temp_hist->SetName(Form("th1f_%s_%s_param%dDown",dataSetName.c_str(),res_ptr->GetName(),par));
    v_th1f_.push_back(*temp_hist);
 
    // after each parameter we reset the originals back
    std::cout << "Reset to Original Values " << std::endl;
    setArgSetParameters(rooParameters,original_values);
  } 
  std::cout << "Generated systematic varitions from " << res_ptr->GetName()<<std::endl;    
}

#endif
