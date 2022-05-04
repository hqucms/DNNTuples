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
  PFCompleteFiller() : PFCompleteFiller("", 0.4, 0.4) {}
  PFCompleteFiller(std::string branchName, double jetR=0.4, double pfcandR=0.4) : NtupleBase(branchName, jetR, pfcandR) {}
  virtual ~PFCompleteFiller() {}

  // get input parameters from the cfg file
  virtual void readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector && cc) override;

  // read event content or event setup for each event
  virtual void readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

protected:
  // declare the data branches (name, type, default values)
  virtual void book() override;
  // fill the branches
  virtual bool fill(const reco::VertexCompositePtrCandidate &sv, size_t svidx, const edm::Handle<edm::View<reco::Candidate>> &candHandle) override;

private:
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::Handle<reco::VertexCollection> vertices;

  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
  edm::Handle<reco::VertexCompositePtrCandidateCollection> SVs;

  edm::ESHandle<TransientTrackBuilder> builder_;
};

} /* namespace deepntuples */

#endif /* NTUPLER_INTERFACE_PFCOMPLETEFILLER_H_ */
