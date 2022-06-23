/*
 * PFCompleteFiller.cc
 *
 *  Created on: Sep 25, 2017
 *      Author: hqu
 */

#include <unordered_map>
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DeepNTuples/Ntupler/interface/PFCompleteFiller.h"

namespace deepntuples {

void PFCompleteFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  transientTrackBuilderToken_ = cc.esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));
  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
}

void PFCompleteFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);
  builder_ = iSetup.getHandle(transientTrackBuilderToken_);
}

void PFCompleteFiller::book() {

  data.add<int>("n_pfcands", 0);
  data.add<float>("npfcands", 0);

  // no puppi scaled
  data.addMulti<float>("pfcand_px");
  data.addMulti<float>("pfcand_py");
  data.addMulti<float>("pfcand_pz");
  data.addMulti<float>("pfcand_energy");

  data.addMulti<float>("pfcand_pt_nopuppi");
  data.addMulti<float>("pfcand_pt_log_nopuppi");
  data.addMulti<float>("pfcand_e_log_nopuppi");

  data.addMulti<float>("pfcand_phirel");
  data.addMulti<float>("pfcand_etarel");
  data.addMulti<float>("pfcand_puppiw");
  data.addMulti<float>("pfcand_abseta");

  data.addMulti<float>("pfcand_charge");
  data.addMulti<float>("pfcand_isMu");
  data.addMulti<float>("pfcand_isEl");
  data.addMulti<float>("pfcand_isChargedHad");
  data.addMulti<float>("pfcand_isGamma");
  data.addMulti<float>("pfcand_isNeutralHad");

  // for neutral
  data.addMulti<float>("pfcand_hcalFrac");
  data.addMulti<float>("pfcand_hcalFracCalib");

  // for charged
  data.addMulti<float>("pfcand_VTX_ass");
  data.addMulti<float>("pfcand_fromPV");
  data.addMulti<float>("pfcand_lostInnerHits");
  data.addMulti<float>("pfcand_trackHighPurity");

  // impact parameters
  data.addMulti<float>("pfcand_dz");
  data.addMulti<float>("pfcand_dzsig");
  data.addMulti<float>("pfcand_dxy");
  data.addMulti<float>("pfcand_dxysig");

  // track quality
  data.addMulti<float>("pfcand_normchi2");
  data.addMulti<float>("pfcand_quality");

  data.addMulti<float>("pfcand_qoverp");
  data.addMulti<float>("pfcand_qoverpError");

  data.addMulti<float>("pfcand_nValidPixelHits");
  data.addMulti<float>("pfcand_nValidPixelBarrelHits");
  data.addMulti<float>("pfcand_nValidPixelEndcapHits");
  data.addMulti<float>("pfcand_nValidStripHits");
  data.addMulti<float>("pfcand_nValidStripTIBHits");
  data.addMulti<float>("pfcand_nValidStripTIDHits");
  data.addMulti<float>("pfcand_nValidStripTOBHits");
  data.addMulti<float>("pfcand_nValidStripTECHits");

  // data.addMulti<float>("pfcand_nLostPixelHits");
  // data.addMulti<float>("pfcand_nLostPixelBarrelHits");
  // data.addMulti<float>("pfcand_nLostPixelEndcapHits");
  // data.addMulti<float>("pfcand_nLostStripHits");
  // data.addMulti<float>("pfcand_nLostStripTIBHits");
  // data.addMulti<float>("pfcand_nLostStripTIDHits");
  // data.addMulti<float>("pfcand_nLostStripTOBHits");
  // data.addMulti<float>("pfcand_nLostStripTECHits");

  data.addMulti<float>("pfcand_pixelLayers");
  data.addMulti<float>("pfcand_pixelBarrelLayers");
  data.addMulti<float>("pfcand_pixelEndcapLayers");
  data.addMulti<float>("pfcand_stripLayers");
  data.addMulti<float>("pfcand_stripTIBLayers");
  data.addMulti<float>("pfcand_stripTIDLayers");
  data.addMulti<float>("pfcand_stripTOBLayers");
  data.addMulti<float>("pfcand_stripTECLayers");

  // track btag info
  // data.addMulti<float>("pfcand_btagMomentum");
  // data.addMulti<float>("pfcand_btagEta");
  data.addMulti<float>("pfcand_btagEtaRel");
  data.addMulti<float>("pfcand_btagPtRel");
  // data.addMulti<float>("pfcand_btagPPar");
  // data.addMulti<float>("pfcand_btagDeltaR");
  data.addMulti<float>("pfcand_btagPtRatio");
  data.addMulti<float>("pfcand_btagPParRatio");
  data.addMulti<float>("pfcand_btagSip3dVal");
  data.addMulti<float>("pfcand_btagSip3dSig");
  data.addMulti<float>("pfcand_btagJetDistVal");
//  data.addMulti<float>("pfcand_btagJetDistSig"); // always gives 0?

}

bool PFCompleteFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  const auto& pfCands = jet_helper.getJetConstituents();

  data.fill<int>("n_pfcands", pfCands.size());
  data.fill<float>("npfcands", pfCands.size());

  float etasign = jet.eta()>0 ? 1 : -1;

  for (const auto& cand : pfCands){

    const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));

    // basic kinematics, valid for both charged and neutral
    // not puppi weighted
    data.fillMulti<float>("pfcand_px", packed_cand->px());
    data.fillMulti<float>("pfcand_py", packed_cand->py());
    data.fillMulti<float>("pfcand_pz", packed_cand->pz());
    data.fillMulti<float>("pfcand_energy", packed_cand->energy());

    data.fillMulti<float>("pfcand_pt_nopuppi", packed_cand->pt());
    data.fillMulti<float>("pfcand_pt_log_nopuppi", catchInfs(std::log(packed_cand->pt()), -99));
    data.fillMulti<float>("pfcand_e_log_nopuppi", catchInfs(std::log(packed_cand->energy()), -99));

    data.fillMulti<float>("pfcand_phirel", reco::deltaPhi(*packed_cand, jet));
    data.fillMulti<float>("pfcand_etarel", etasign * (packed_cand->eta() - jet.eta()));
    data.fillMulti<float>("pfcand_abseta", std::abs(packed_cand->eta()));

    data.fillMulti<float>("pfcand_puppiw", jet_helper.getPuppiWeight(cand));

    data.fillMulti<float>("pfcand_charge", packed_cand->charge());
    data.fillMulti<float>("pfcand_isEl", std::abs(packed_cand->pdgId())==11);
    data.fillMulti<float>("pfcand_isMu", std::abs(packed_cand->pdgId())==13);
    data.fillMulti<float>("pfcand_isChargedHad", std::abs(packed_cand->pdgId())==211);
    data.fillMulti<float>("pfcand_isGamma", std::abs(packed_cand->pdgId())==22);
    data.fillMulti<float>("pfcand_isNeutralHad", std::abs(packed_cand->pdgId())==130);

    // for neutral
    float hcal_fraction = 0.;
    if (packed_cand->pdgId() == 1 || packed_cand->pdgId() == 130) {
      hcal_fraction = packed_cand->hcalFraction();
    } else if (packed_cand->isIsolatedChargedHadron()) {
      hcal_fraction = packed_cand->rawHcalFraction();
    }
    data.fillMulti<float>("pfcand_hcalFrac", hcal_fraction);
    data.fillMulti<float>("pfcand_hcalFracCalib", packed_cand->hcalFraction());

    // for charged
    data.fillMulti<float>("pfcand_VTX_ass", packed_cand->pvAssociationQuality());
    data.fillMulti<float>("pfcand_fromPV", packed_cand->fromPV());
    data.fillMulti<float>("pfcand_lostInnerHits", packed_cand->lostInnerHits());
    data.fillMulti<float>("pfcand_trackHighPurity", packed_cand->trackHighPurity());

    // impact parameters
    data.fillMulti<float>("pfcand_dz", catchInfs(packed_cand->dz()));
    data.fillMulti<float>("pfcand_dzsig", packed_cand->bestTrack() ? catchInfs(packed_cand->dz()/packed_cand->dzError()) : 0);
    data.fillMulti<float>("pfcand_dxy", catchInfs(packed_cand->dxy()));
    data.fillMulti<float>("pfcand_dxysig", packed_cand->bestTrack() ? catchInfs(packed_cand->dxy()/packed_cand->dxyError()) : 0);

    if (packed_cand->bestTrack()){
      const auto *trk = packed_cand->bestTrack();
      data.fillMulti<float>("pfcand_normchi2", catchInfs(trk->normalizedChi2()));
      data.fillMulti<float>("pfcand_quality", trk->qualityMask());

      data.fillMulti<float>("pfcand_qoverp", trk->qoverp());
      data.fillMulti<float>("pfcand_qoverpError", trk->qoverpError());

      data.fillMulti<float>("pfcand_nValidPixelHits", trk->hitPattern().numberOfValidPixelHits());
      data.fillMulti<float>("pfcand_nValidPixelBarrelHits", trk->hitPattern().numberOfValidPixelBarrelHits());
      data.fillMulti<float>("pfcand_nValidPixelEndcapHits", trk->hitPattern().numberOfValidPixelEndcapHits());
      data.fillMulti<float>("pfcand_nValidStripHits", trk->hitPattern().numberOfValidStripHits());
      data.fillMulti<float>("pfcand_nValidStripTIBHits", trk->hitPattern().numberOfValidStripTIBHits());
      data.fillMulti<float>("pfcand_nValidStripTIDHits", trk->hitPattern().numberOfValidStripTIDHits());
      data.fillMulti<float>("pfcand_nValidStripTOBHits", trk->hitPattern().numberOfValidStripTOBHits());
      data.fillMulti<float>("pfcand_nValidStripTECHits", trk->hitPattern().numberOfValidStripTECHits());

      // data.fillMulti<float>("pfcand_nLostPixelHits", trk->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS));
      // data.fillMulti<float>("pfcand_nLostPixelBarrelHits", trk->hitPattern().numberOfLostPixelBarrelHits(reco::HitPattern::MISSING_INNER_HITS));
      // data.fillMulti<float>("pfcand_nLostPixelEndcapHits", trk->hitPattern().numberOfLostPixelEndcapHits(reco::HitPattern::MISSING_OUTER_HITS));
      // data.fillMulti<float>("pfcand_nLostStripHits", trk->hitPattern().numberOfLostStripHits(reco::HitPattern::TRACK_HITS));
      // data.fillMulti<float>("pfcand_nLostStripTIBHits", trk->hitPattern().numberOfLostStripTIBHits(reco::HitPattern::TRACK_HITS));
      // data.fillMulti<float>("pfcand_nLostStripTIDHits", trk->hitPattern().numberOfLostStripTIDHits(reco::HitPattern::TRACK_HITS));
      // data.fillMulti<float>("pfcand_nLostStripTOBHits", trk->hitPattern().numberOfLostStripTOBHits(reco::HitPattern::TRACK_HITS));
      // data.fillMulti<float>("pfcand_nLostStripTECHits", trk->hitPattern().numberOfLostStripTECHits(reco::HitPattern::TRACK_HITS));

      data.fillMulti<float>("pfcand_pixelLayers", trk->hitPattern().pixelLayersWithMeasurement());
      data.fillMulti<float>("pfcand_pixelBarrelLayers", trk->hitPattern().pixelBarrelLayersWithMeasurement());
      data.fillMulti<float>("pfcand_pixelEndcapLayers", trk->hitPattern().pixelEndcapLayersWithMeasurement());
      data.fillMulti<float>("pfcand_stripLayers", trk->hitPattern().stripLayersWithMeasurement());
      data.fillMulti<float>("pfcand_stripTIBLayers", trk->hitPattern().stripTIBLayersWithMeasurement());
      data.fillMulti<float>("pfcand_stripTIDLayers", trk->hitPattern().stripTIDLayersWithMeasurement());
      data.fillMulti<float>("pfcand_stripTOBLayers", trk->hitPattern().stripTOBLayersWithMeasurement());
      data.fillMulti<float>("pfcand_stripTECLayers", trk->hitPattern().stripTECLayersWithMeasurement());

    } else {
      data.fillMulti<float>("pfcand_normchi2", 999);
      data.fillMulti<float>("pfcand_quality", 0);

      data.fillMulti<float>("pfcand_qoverp", 0);
      data.fillMulti<float>("pfcand_qoverpError", 0);

      data.fillMulti<float>("pfcand_nValidPixelHits", -1);
      data.fillMulti<float>("pfcand_nValidPixelBarrelHits", -1);
      data.fillMulti<float>("pfcand_nValidPixelEndcapHits", -1);
      data.fillMulti<float>("pfcand_nValidStripHits", -1);
      data.fillMulti<float>("pfcand_nValidStripTIBHits", -1);
      data.fillMulti<float>("pfcand_nValidStripTIDHits", -1);
      data.fillMulti<float>("pfcand_nValidStripTOBHits", -1);
      data.fillMulti<float>("pfcand_nValidStripTECHits", -1);

      // data.fillMulti<float>("pfcand_nLostPixelHits", -5);
      // data.fillMulti<float>("pfcand_nLostPixelBarrelHits", -5);
      // data.fillMulti<float>("pfcand_nLostPixelEndcapHits", -5);
      // data.fillMulti<float>("pfcand_nLostStripHits", -5);
      // data.fillMulti<float>("pfcand_nLostStripTIBHits", -5);
      // data.fillMulti<float>("pfcand_nLostStripTIDHits", -5);
      // data.fillMulti<float>("pfcand_nLostStripTOBHits", -5);
      // data.fillMulti<float>("pfcand_nLostStripTECHits", -5);

      data.fillMulti<float>("pfcand_pixelLayers", -1);
      data.fillMulti<float>("pfcand_pixelBarrelLayers", -1);
      data.fillMulti<float>("pfcand_pixelEndcapLayers", -1);
      data.fillMulti<float>("pfcand_stripLayers", -1);
      data.fillMulti<float>("pfcand_stripTIBLayers", -1);
      data.fillMulti<float>("pfcand_stripTIDLayers", -1);
      data.fillMulti<float>("pfcand_stripTOBLayers", -1);
      data.fillMulti<float>("pfcand_stripTECLayers", -1);
    }

    // build track info map
    TrackInfoBuilder trkinfo;
    trkinfo.buildTrackInfo(builder_, *packed_cand, jet, vertices->at(0));

    // data.fillMulti<float>("pfcand_btagMomentum", catchInfs(trkinfo.getTrackMomentum()));
    // data.fillMulti<float>("pfcand_btagEta", catchInfs(trkinfo.getTrackEta()));
    data.fillMulti<float>("pfcand_btagEtaRel", catchInfs(trkinfo.getTrackEtaRel()));
    data.fillMulti<float>("pfcand_btagPtRel", catchInfs(trkinfo.getTrackPtRel()));
    // data.fillMulti<float>("pfcand_btagPPar", catchInfs(trkinfo.getTrackPPar()));
    // data.fillMulti<float>("pfcand_btagDeltaR", catchInfs(trkinfo.getTrackDeltaR()));
    data.fillMulti<float>("pfcand_btagPtRatio", catchInfs(trkinfo.getTrackPtRatio()));
    data.fillMulti<float>("pfcand_btagPParRatio", catchInfs(trkinfo.getTrackPParRatio()));
    data.fillMulti<float>("pfcand_btagSip3dVal", catchInfs(trkinfo.getTrackSip3dVal()));
    data.fillMulti<float>("pfcand_btagSip3dSig", catchInfs(trkinfo.getTrackSip3dSig()));
    data.fillMulti<float>("pfcand_btagJetDistVal", catchInfs(trkinfo.getTrackJetDistVal()));
//    data.fillMulti<float>("pfcand_btagJetDistSig", catchInfs(trkinfo.getTrackJetDistSig()));


  }


  return true;
}

} /* namespace deepntuples */
