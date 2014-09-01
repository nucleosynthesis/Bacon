// Some statistics calculations
#include "TH1F.h"
#include "TMath.h"
#include "TGraph.h"
#include <vector>
#include <iostream>

double calculateExpectedSignificance(std::vector<double> signal, std::vector<double> background){
  
  double sterm=0;
  double logterms=0;
  unsigned int nchannel  = signal.size();

  if (background.size()!=nchannel) { 
  	std::cout << " background and signal vectors should be the same size!" << std::endl;
	return 0;
  }

  for (unsigned int i=0;i<nchannel;i++){
    if (!signal[i]>0) continue;
    logterms+=(signal[i]+background[i])*TMath::Log((signal[i]+background[i])/background[i]);
    sterm+=signal[i];
  }
  double sig =  1.4142*TMath::Sqrt(logterms - sterm);
  //std::cout << "Significance @ " << nchannel << " bins = " << sig << std::endl; 
  return sig;
}

double calculateExpectedCLs(double mu, std::vector<double> signal, std::vector<double> background){
  double logterms=0;
  unsigned int nchannel  = signal.size();
  if (background.size()!=nchannel) { 
  	std::cout << " background and signal vectors should be the same size!" << std::endl;
	return 0;
  }
  for (unsigned int i=0;i<nchannel;i++){
    if (!background[i]>0) continue;
    double bi = background[i];
    double si = mu*signal[i];
    logterms+= bi*(TMath::Log(si+bi) - TMath::Log(bi)) - si;  
  }
  double qmu = -2*logterms;
  //std::cout << qmu <<" ";
  double CLs = TMath::Prob(qmu,1);
  return CLs;
}

double calculateExpectedLimit(double rlow, double rhigh, std::vector<double> signal, std::vector<double> background, int np=20){
  TGraph *gr = new TGraph();
  double step = (rhigh-rlow)/np;
  double mu = rlow;
  for (int i=0;i<np;i++){
    double CLs = calculateExpectedCLs(mu,signal,background);
    gr->SetPoint(i,CLs,mu); 
    //std::cout << "At mu = "<<mu<< ", CLs = " << CLs << std::endl;
    mu+=step;
  }
  return gr->Eval(0.05);
}

double calculateExpectedLimit(double rlow,double rhigh, TH1F *signal, TH1F *background, int np=20){
  int nchannel = signal->GetNbinsX();
   
  std::vector<double> signal_v;
  std::vector<double> background_v;
  if (background->GetNbinsX()!=nchannel) { 
  	std::cout << " background and signal vectors should be the same size!" << std::endl;
	return 0;
  }
  for (int b=1;b<=nchannel;b++){
    signal_v.push_back(signal->GetBinContent(b));
    background_v.push_back(background->GetBinContent(b));
  }
  double lim = calculateExpectedLimit(rlow,rhigh,signal_v,background_v,np);
  return lim;
}


double phi(double x, double y){
  if (x==0 && y==0) return 0;
  return x*TMath::Log(y) - y;
}
double calculateSignificanceBkgUncert(std::vector<double> signal, std::vector<double> background, std::vector<double> backgroundE){
  // assuming a global tau, ideally should have one per bin
  // tau == nMC/nBackground
  unsigned int nchannel  = signal.size();
  if (background.size()!=nchannel) { 
  	std::cout << " background and signal vectors should be the same size!" << std::endl;
	return 0;
  }
  double logterms=0;
  for (unsigned int i=0;i<nchannel;i++){
    double s = signal[i];
    if (!s>0) continue;
    double b = background[i];
    double sigma = backgroundE[i];
    double tau = b/(sigma*sigma);
 //   std::cout << "\ni, tau, b, err, sqrt(b) = " << i << " " << tau << " " << b<<  " " <<sigma<< " " <<TMath::Sqrt(b) <<  std::endl;
    double m = tau*b;
    double n = s + b;
    double bh = (n + m)/(1+tau);
    double ln0=phi(m,tau*bh) + phi(n,bh) - phi(m,tau*b) - phi(n,s+b);
    //std::cout << " bin i, ln0 = " << ln0 << std::endl;
    
    logterms+=ln0;
  }
  return TMath::Sqrt(-2*logterms);
}

double calculateExpectedSignificance(TH1F *signal, TH1F *background, bool incbkgsys=false){

  int nchannel = signal->GetNbinsX();
   
  if (background->GetNbinsX()!=nchannel) { 
  	std::cout << " background and signal vectors should be the same size!" << std::endl;
	return 0;
  }

  std::vector<double> signal_v;
  std::vector<double> background_v;
  std::vector<double> background_vE;

  for (int b=1;b<=nchannel;b++){
    signal_v.push_back(signal->GetBinContent(b));
    background_v.push_back(background->GetBinContent(b));
    background_vE.push_back(background->GetBinError(b));
  }

  
  double sig =0; 
  if (incbkgsys){
	//double tau = background->GetEntries()/background->Integral();  
	sig = calculateSignificanceBkgUncert(signal_v,background_v,background_vE);
  }
  else sig = calculateExpectedSignificance(signal_v,background_v);
  return sig;
   
}

