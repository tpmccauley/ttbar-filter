// -*- C++ -*-
//
// Package:    TTBarFilter
// Class:      TTBarFilter
// 
/**\class TTBarFilter TTBarFilter.cc TTBarFilter/TTBarFilter/src/TTBarFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  thomas.mccauley@cern.ch
//         Created:  Wed Apr 13 13:35:17 CEST 2016
// $Id$
//
//

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include <iostream>
#include <string>
#include <fstream>

bool 
dR(std::vector<std::pair<double,double> >& etaphi, double& meta, double& mphi) 
{
  double dR;
  
  for ( std::vector<std::pair<double,double> >::iterator iep = etaphi.begin(), epend = etaphi.end();
        iep != epend; ++iep ) 
  {
    dR = (iep->first - meta)*(iep->first - meta);
    dR += (iep->second - mphi)*(iep->second - mphi);
    dR = sqrt(dR);
    
    if ( dR < 0.3 )
      return false;
  }
 
  return true;
}

class TTBarFilter : public edm::EDFilter 
{
public:
  explicit TTBarFilter(const edm::ParameterSet&);
  ~TTBarFilter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual bool beginRun(edm::Run&, edm::EventSetup const&);
  virtual bool endRun(edm::Run&, edm::EventSetup const&);
  virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  edm::InputTag muonInputTag_;
  edm::InputTag pfJetInputTag_;

  std::ofstream csvOut_;
  std::string csvFileName_;
  double minMuonPt_;
  double maxMuonEta_;

  double minJetPt_;
  unsigned int minNJets_;
  double maxJetEta_;
};

TTBarFilter::TTBarFilter(const edm::ParameterSet& iConfig)
  : muonInputTag_(iConfig.getParameter<edm::InputTag>("muonInputTag")),
    pfJetInputTag_(iConfig.getParameter<edm::InputTag>("pfJetInputTag")),
    csvFileName_(iConfig.getParameter<std::string>("csvFileName")),
    minMuonPt_(iConfig.getParameter<double>("minMuonPt")),
    maxMuonEta_(iConfig.getParameter<double>("maxMuonEta")),
    minJetPt_(30.0), minNJets_(4), maxJetEta_(2.4) 
{    
  csvOut_.open(csvFileName_.c_str());
}

TTBarFilter::~TTBarFilter()
{}

bool
TTBarFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel(muonInputTag_, muons);

  if ( ! muons.isValid() ) 
  {    
    std::cerr<<"SingleMuonFilter: invalid muon collection"<<std::endl;
    return false;
  }

  edm::Handle<reco::PFJetCollection> pfJets;
  iEvent.getByLabel(pfJetInputTag_, pfJets);

  if ( ! pfJets.isValid() ) 
  {  
    std::cerr<<"SingleMuonFilter: invalid PFJet collection"<<std::endl;
    return false;
  }
  
  if ( pfJets->size() < minNJets_ ) 
    return false;

  double jetPt, jetEta, jetPhi;
  unsigned int nJets = 0;
  std::vector<std::pair<double,double> > jetEtaPhi;

  for ( reco::PFJetCollection::const_iterator ij = pfJets->begin(), end = pfJets->end();
        ij != end; ++ ij ) 
  {    
    jetPt = ij->et();
    jetEta = ij->eta();
    jetPhi = ij->phi();
    
    if ( jetPt > minJetPt_ && fabs(jetEta) < maxJetEta_ ) 
    {   
      nJets++;
      jetEtaPhi.push_back(std::make_pair(jetEta,jetPhi));
    }
  }
 
  if ( nJets < minNJets_ ) 
    return false;

  double pt, eta, phi;
  double relIso;
  
  unsigned int nLoose = 0;
  unsigned int nTight = 0;

  for ( reco::MuonCollection::const_iterator it = muons->begin(), end = muons->end();
        it != end; ++it ) 
  {
    if ( (*it).isGlobalMuon() ) 
    {   
      pt = (*it).globalTrack()->pt();      
      eta = (*it).globalTrack()->eta();
      phi = (*it).globalTrack()->phi();
      relIso = ((*it).isolationR03().sumPt + (*it).isolationR03().emEt + (*it).isolationR03().hadEt) / pt;
   
      if ( (*it).combinedMuon()->normalizedChi2() < 10 && 
           (*it).combinedMuon()->hitPattern().numberOfValidMuonHits() > 0 &&
           pt > minMuonPt_ &&  
           fabs(eta) < maxMuonEta_ && 
           relIso < 0.05 && 
           dR(jetEtaPhi, eta, phi) ) 
      {
        nTight++;
        continue;
      }
      
      if ( pt > 10 && relIso < 0.2 ) 
      {  
        nLoose++;
      }
    }
  }
  
  if ( nTight == 1 && nLoose == 0 ) 
  {  
    csvOut_<<"'"<< iEvent.id().run() <<":"<< iEvent.id().event() <<"'"<<std::endl;
    return true;
  }

  return false;
}

void 
TTBarFilter::beginJob()
{}

void 
TTBarFilter::endJob() {
}

bool 
TTBarFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

bool 
TTBarFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

bool 
TTBarFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

bool 
TTBarFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

void
TTBarFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(TTBarFilter);
