#include <cmath>
#include "TVector3.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "DeepNTuples/NtupleCommons/interface/InfinityCatcher.h"
#include "DeepNTuples/BTagHelpers/interface/TrackInfoBuilder.h"

namespace deepntuples {

//void TrackInfoBuilder::buildTrackInfo(const edm::ESHandle<TransientTrackBuilder> builder,
//    const pat::PackedCandidate &pfcand, const pat::Jet &jet, const reco::Vertex& pv) {
// NEW:  Use SVs
void TrackInfoBuilder::buildTrackInfo(const edm::ESHandle<TransientTrackBuilder> builder,
    const pat::PackedCandidate &pfcand, const reco::VertexCompositePtrCandidate &sv, const reco::Vertex& pv) {

  const auto* bestTrk = pfcand.bestTrack();
  if (!bestTrk){
    trackMomentum_ = 0;
    trackEta_ = 0;
    trackEtaRel_ = 0;
    trackPtRel_ = 0;
    trackPPar_ = 0;
    trackDeltaR_ = 0;
    trackPtRatio_ = 0;
    trackPParRatio_ = 0;
    trackSip2dVal_ = 0;
    trackSip2dSig_ = 0;
    trackSip3dVal_ = 0;
    trackSip3dSig_ = 0;
    trackJetDistVal_ = 0;
    trackJetDistSig_ = 0;
    trackDecayLengthVal_ = 0;
    trackDecayLengthSig_ = 0;
    return;
  }

  // DataFormats/BTauReco/interface/IPTagInfo.h

  math::XYZVector svDir = sv.momentum().Unit();
  GlobalVector svXYZVector(sv.px(),sv.py(),sv.pz());

  const auto &trk = pfcand.pseudoTrack();
  reco::TransientTrack transientTrack(builder->build(trk));
  Measurement1D meas_ip2d = IPTools::signedTransverseImpactParameter(transientTrack, svXYZVector, pv).second;
  Measurement1D meas_ip3d = IPTools::signedImpactParameter3D(transientTrack, svXYZVector, pv).second;
  Measurement1D svdist = IPTools::jetTrackDistance(transientTrack, svXYZVector, pv).second;
  Measurement1D decay_length = IPTools::signedDecayLength3D(transientTrack, svXYZVector, pv).second;
  math::XYZVector trackMom = trk.momentum();
  double trackMag = std::sqrt(trackMom.Mag2());

  TVector3 trackMom3(trackMom.x(),trackMom.y(),trackMom.z());
  TVector3 svDir3(svDir.x(),svDir.y(),svDir.z());

  trackMomentum_ = catchInfs(trackMag);
  trackEta_ = catchInfs(trackMom.Eta());
  trackEtaRel_ = catchInfs(reco::btau::etaRel(svDir, trackMom));
  trackPtRel_ = catchInfs(trackMom3.Perp(svDir3));
  trackPPar_ = catchInfs(svDir.Dot(trackMom));
  trackDeltaR_ = catchInfs(reco::deltaR(trackMom, svDir));
  trackPtRatio_ = catchInfs(trackMom3.Perp(svDir3) / trackMag);
  trackPParRatio_ = catchInfs(svDir.Dot(trackMom) / trackMag);
  trackSip2dVal_ = catchInfs(meas_ip2d.value());
  trackSip2dSig_ = catchInfs(meas_ip2d.significance());
  trackSip3dVal_ = catchInfs(meas_ip3d.value());
  trackSip3dSig_ = catchInfs(meas_ip3d.significance());
  trackJetDistVal_ = catchInfs(svdist.value());
  trackJetDistSig_ = catchInfs(svdist.significance());
  trackDecayLengthVal_ = catchInfs(decay_length.value());
  trackDecayLengthSig_ = catchInfs(decay_length.significance());
}

} /* namespace deepntuples */

