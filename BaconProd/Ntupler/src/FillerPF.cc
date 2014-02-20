#include "BaconProd/Ntupler/interface/FillerPF.hh"
#include "BaconAna/DataFormats/interface/TPFPart.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <TClonesArray.h>
#include <TMath.h>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerPF::FillerPF():
  fPFName("particleFlow"),
  fPVName("offlinePrimaryVertices")
{
}

//--------------------------------------------------------------------------------------------------
FillerPF::~FillerPF(){}

//--------------------------------------------------------------------------------------------------
void FillerPF::fill(TClonesArray *array,TClonesArray *iVtxCol,
		    const edm::Event &iEvent) 
		   
{
  assert(array);
  // Get PF collection
  edm::Handle<reco::PFCandidateCollection> hPFProduct;
  iEvent.getByLabel(fPFName,hPFProduct);
  assert(hPFProduct.isValid());
  const reco::PFCandidateCollection *PFCol = hPFProduct.product();

  // Get vertex collection
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByLabel(fPVName,hVertexProduct);
  assert(hVertexProduct.isValid());
  const reco::VertexCollection *pvCol = hVertexProduct.product();

  TClonesArray &rArray = *array;
  int pId = 0; 
  for(reco::PFCandidateCollection::const_iterator itPF = PFCol->begin(); itPF!=PFCol->end(); itPF++) {
    pId++;
    // construct object and place in array
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TPFPart();
    baconhep::TPFPart *pPF = (baconhep::TPFPart*)rArray[index];

    //
    // Kinematics
    //==============================    
    pPF->pt  = itPF->pt();
    pPF->eta = itPF->eta();
    pPF->phi = itPF->phi();
    pPF->m   = itPF->mass();
    pPF->e   = itPF->energy();
    pPF->q   = itPF->charge();
    pPF->pfType = itPF->particleId();
    //pPF->deltaP = itPF->deltaP();
    pPF->ecalE  = itPF->ecalEnergy();
    pPF->hcalE  = itPF->hcalEnergy();
    //
    // TrackInfo
    //==============================
    const reco::TrackRef& pfTrack = itPF->trackRef();
    if(!pfTrack.isNonnull()) continue;
    int    ndof     = pfTrack->ndof();
    double chi2     = pfTrack->chi2();
    pPF->trkChi2 = TMath::Prob(chi2,ndof);
      const reco::Vertex *closestVtx = 0;
    for(reco::VertexCollection::const_iterator iV = pvCol->begin(); iV!=pvCol->end(); ++iV) {
      if(iV->trackWeight(itPF->trackRef())>0) {
	closestVtx  = &(*iV);
	break;
      }
    }
    if(closestVtx == 0) continue;
    pPF->vtxChi2 = closestVtx->trackWeight(itPF->trackRef());
    int lId = -1;
    for(int i0 = 0; i0 < iVtxCol->GetEntries(); i0++) { 
      baconhep::TVertex* pVertex = (TVertex*)(*iVtxCol)[i0];
      if(fabs(closestVtx->x() - pVertex->x) + 
	 fabs(closestVtx->y() - pVertex->y) + 
	 fabs(closestVtx->z() - pVertex->z) > 0.0001) continue;
      lId = i0;
      break;
    }
    pPF->vtxId = lId;
  } 
}
