/*
 * PFCompleteFiller.h
 *
 *  Created on: Sep 25, 2017
 *      Author: hqu
 */

#ifndef NTUPLER_INTERFACE_PFCOMPLETEFILLER_H_
#define NTUPLER_INTERFACE_PFCOMPLETEFILLER_H_

#include "DeepNTuples/BTagHelpers/interface/TrackInfoBuilder.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"

namespace deepntuples {

class PFCompleteFiller: public NtupleBase {
public:
  PFCompleteFiller() : PFCompleteFiller("", 0.8) {}
  PFCompleteFiller(std::string branchName, double jetR=0.8) : NtupleBase(branchName, jetR) {}
  virtual ~PFCompleteFiller() {}

  // get input parameters from the cfg file
  virtual void readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector && cc) override;

  // read event content or event setup for each event
  virtual void readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

protected:
  // declare the data branches (name, type, default values)
  virtual void book() override;
  // fill the branches
  virtual bool fill(const pat::Jet &jet, size_t jetidx, const JetHelper &jet_helper) override;

private:
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::Handle<reco::VertexCollection> vertices;

  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
  edm::Handle<reco::VertexCompositePtrCandidateCollection> SVs;

  edm::ESHandle<TransientTrackBuilder> builder_;

  // muon
  edm::EDGetTokenT<edm::View<pat::Muon>>           muonToken_;
  edm::Handle<edm::View<pat::Muon>>                muons;

  // electron
  edm::EDGetTokenT<edm::View<pat::Electron>>       electronToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >           vetoIdToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >           looseIdToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >           mediumIdToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >           tightIdToken_;
  edm::Handle<edm::View<pat::Electron>>            electrons;
  edm::Handle<edm::ValueMap<bool> >                veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> >                loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> >                medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> >                tight_id_decisions;

  bool fillElectronVars_ = false;
  bool fillMuonVars_ = false;
};

} /* namespace deepntuples */

#endif /* NTUPLER_INTERFACE_PFCOMPLETEFILLER_H_ */
