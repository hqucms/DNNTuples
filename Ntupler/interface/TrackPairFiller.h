#ifndef NTUPLER_INTERFACE_TRACKPAIRFILLER_H_
#define NTUPLER_INTERFACE_TRACKPAIRFILLER_H_

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoBTag/FeatureTools/interface/TrackPairInfoBuilder.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"

namespace deepntuples {

  class TrackPairFiller : public NtupleBase {
  public:
    TrackPairFiller() : TrackPairFiller("", 0.8) {}
    TrackPairFiller(std::string branchName, double jetR = 0.8) : NtupleBase(branchName, jetR) {}
    virtual ~TrackPairFiller() {}

    // get input parameters from the cfg file
    virtual void readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) override;

    // read event content or event setup for each event
    virtual void readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  protected:
    // declare the data branches (name, type, default values)
    virtual void book() override;
    // fill the branches
    virtual bool fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) override;

  private:
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilderToken_;
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::Handle<reco::VertexCollection> vertices;

    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
    edm::Handle<reco::VertexCompositePtrCandidateCollection> SVs;

    edm::ESHandle<TransientTrackBuilder> builder_;
  };

} /* namespace deepntuples */

#endif /* NTUPLER_INTERFACE_TRACKPAIRFILLER_H_ */