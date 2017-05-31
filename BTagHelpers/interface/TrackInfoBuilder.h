#ifndef BTAGHELPERS_INTERFACE_TRACKINFOBUILDER_H_
#define BTAGHELPERS_INTERFACE_TRACKINFOBUILDER_H_

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

namespace deepntuples {

class TrackInfoBuilder {
public:
  TrackInfoBuilder() {}
  virtual ~TrackInfoBuilder() {}

public:
  void buildTrackInfo(const edm::ESHandle<TransientTrackBuilder> builder,
      const pat::PackedCandidate &pfcand, const pat::Jet &jet, const reco::Vertex &pv);

  const float& getTrackMomentum() const {return trackMomentum_;}
  const float& getTrackEta() const {return trackEta_;}
  const float& getTrackEtaRel() const {return trackEtaRel_;}
  const float& getTrackPtRel() const {return trackPtRel_;}
  const float& getTrackPPar() const {return trackPPar_;}
  const float& getTrackDeltaR() const {return trackDeltaR_;}
  const float& getTrackPtRatio() const {return trackPtRatio_;}
  const float& getTrackPParRatio() const {return trackPParRatio_;}
  const float& getTrackSip2dVal() const {return trackSip2dVal_;}
  const float& getTrackSip2dSig() const {return trackSip2dSig_;}
  const float& getTrackSip3dVal() const {return trackSip3dVal_;}
  const float& getTrackSip3dSig() const {return trackSip3dSig_;}
  const float& getTrackJetDistVal() const {return trackJetDistVal_;}
  const float& getTrackJetDistSig() const {return trackJetDistSig_;}


private:

  float trackMomentum_ = 0;
  float trackEta_ = 0;
  float trackEtaRel_ = 0;
  float trackPtRel_ = 0;
  float trackPPar_ = 0;
  float trackDeltaR_ = 0;
  float trackPtRatio_ = 0;
  float trackPParRatio_ = 0;
  float trackSip2dVal_ = 0;
  float trackSip2dSig_ = 0;
  float trackSip3dVal_ = 0;
  float trackSip3dSig_ = 0;
  float trackJetDistVal_ = 0;
  float trackJetDistSig_ = 0;

};

} /* namespace deepntuples */

#endif /* BTAGHELPERS_INTERFACE_TRACKINFOBUILDER_H_ */
