/*
 * TrackFiller.h
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#ifndef NTUPLEAK8_INTERFACE_TRACKFILLER_H_
#define NTUPLEAK8_INTERFACE_TRACKFILLER_H_

#include "DeepNTuples/BTagHelpers/interface/TrackInfoBuilder.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"

namespace deepntuples {

class TrackFiller: public NtupleBase {
public:
  TrackFiller() : TrackFiller("", 0.8) {}
  TrackFiller(std::string branchName, double jetR=0.8) : NtupleBase(branchName, jetR) {}
  virtual ~TrackFiller() {}

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

};

} /* namespace deepntuples */

#endif /* NTUPLEAK8_INTERFACE_TRACKFILLER_H_ */
