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
  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
}

void PFCompleteFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder_);
}

void PFCompleteFiller::book() {

  data.add<int>("n_pfcands", 0);
  data.add<float>("npfcands", 0);

  // no puppi scaled
  data.addMulti<float>("pfcand_pt_nopuppi");
  data.addMulti<float>("pfcand_pt_log_nopuppi");
  data.addMulti<float>("pfcand_e_log_nopuppi");

  data.addMulti<float>("pfcand_phirel");
  data.addMulti<float>("pfcand_etarel");
  data.addMulti<float>("pfcand_deltaR");
  data.addMulti<float>("pfcand_puppiw");
  data.addMulti<float>("pfcand_abseta");

  data.addMulti<float>("pfcand_drminsvin"); // restricted to within the jet cone

  data.addMulti<float>("pfcand_charge");
  data.addMulti<float>("pfcand_isMu");
  data.addMulti<float>("pfcand_isEl");
  data.addMulti<float>("pfcand_isChargedHad");
  data.addMulti<float>("pfcand_isGamma");
  data.addMulti<float>("pfcand_isNeutralHad");

  data.addMulti<float>("pfcand_px");
  data.addMulti<float>("pfcand_py");
  data.addMulti<float>("pfcand_pz");
  data.addMulti<float>("pfcand_energy");

  // for neutral
  data.addMulti<float>("pfcand_hcalFrac");
  data.addMulti<float>("pfcand_hcalFracCalib");

  // for charged
  data.addMulti<float>("pfcand_VTX_ass");
  data.addMulti<float>("pfcand_fromPV");
  data.addMulti<float>("pfcand_nValidHits");
  data.addMulti<float>("pfcand_nValidPixelHits");
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

  // track covariance
  data.addMulti<float>("pfcand_dptdpt");
  data.addMulti<float>("pfcand_detadeta");
  data.addMulti<float>("pfcand_dphidphi");
  data.addMulti<float>("pfcand_dxydxy");
  data.addMulti<float>("pfcand_dzdz");
  data.addMulti<float>("pfcand_dxydz");
  data.addMulti<float>("pfcand_dphidxy");
  data.addMulti<float>("pfcand_dlambdadz");

  // track btag info
  data.addMulti<float>("pfcand_btagMomentum"); // was commented
  data.addMulti<float>("pfcand_btagEta");  // wc
  data.addMulti<float>("pfcand_btagEtaRel");
  data.addMulti<float>("pfcand_btagPtRel");
  data.addMulti<float>("pfcand_btagPPar");  //wc
  data.addMulti<float>("pfcand_btagDeltaR");  //wc
  data.addMulti<float>("pfcand_btagPtRatio");
  data.addMulti<float>("pfcand_btagPParRatio");
  data.addMulti<float>("pfcand_btagSip2dVal");
  data.addMulti<float>("pfcand_btagSip2dSig");
  data.addMulti<float>("pfcand_btagSip3dVal");
  data.addMulti<float>("pfcand_btagSip3dSig");
  data.addMulti<float>("pfcand_btagJetDistVal");
  data.addMulti<float>("pfcand_btagDecayLengthVal");
  data.addMulti<float>("pfcand_btagDecayLengthSig");

}

//bool PFCompleteFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {
// CHANGING:  Now matching to svs instead of jets
bool PFCompleteFiller::fill(const reco::VertexCompositePtrCandidate &sv, size_t svidx, const edm::Handle<edm::View<reco::Candidate>> &candHandle) {
  std::vector<const pat::PackedCandidate *> pfcands;
  for (const auto &cand : *candHandle) {
    const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&cand);
    double dr = reco::deltaR(*packed_cand, sv);
    if (dr < pfcandR_) {
      pfcands.push_back(packed_cand);
    }
  }
  std::sort(pfcands.begin(), pfcands.end(), [](const pat::PackedCandidate *a, const pat::PackedCandidate *b) {
    return a->pt() > b->pt();
  });
  data.fill<int>("n_pfcands", pfcands.size());
  data.fill<float>("npfcands", pfcands.size());

  // loop over pfcands; find all matches w/in 0.4
  //for (const auto& cand : candHandle) {
  for (const auto *packed_cand : pfcands) {
    float etasign = sv.eta() > 0 ? 1 : -1;
    // save gen vars for each sv:
    // basic kinematics, valid for both charged and neutral
    // not puppi weighted
    data.fillMulti<float>("pfcand_pt_nopuppi", packed_cand->pt());
    data.fillMulti<float>("pfcand_pt_log_nopuppi", catchInfs(std::log(packed_cand->pt()), -99));
    data.fillMulti<float>("pfcand_e_log_nopuppi", catchInfs(std::log(packed_cand->energy()), -99));

    data.fillMulti<float>("pfcand_phirel", reco::deltaPhi(*packed_cand, sv));

    data.fillMulti<float>("pfcand_etarel", etasign * (packed_cand->eta() - sv.eta()));

    //std::cout << " TEMP:  dR is " << reco::deltaR(*packed_cand, sv) << ", jet dR=" << jetR_ << std::endl;

    data.fillMulti<float>("pfcand_deltaR", reco::deltaR(*packed_cand, sv));
    data.fillMulti<float>("pfcand_abseta", std::abs(packed_cand->eta()));

    data.fillMulti<float>("pfcand_puppiw", packed_cand->puppiWeight());

    data.fillMulti<float>("pfcand_charge", packed_cand->charge());
    data.fillMulti<float>("pfcand_isEl", std::abs(packed_cand->pdgId()) == 11);
    data.fillMulti<float>("pfcand_isMu", std::abs(packed_cand->pdgId()) == 13);
    data.fillMulti<float>("pfcand_isChargedHad", std::abs(packed_cand->pdgId()) == 211);
    data.fillMulti<float>("pfcand_isGamma", std::abs(packed_cand->pdgId()) == 22);
    data.fillMulti<float>("pfcand_isNeutralHad", std::abs(packed_cand->pdgId()) == 130);
    data.fillMulti<float>("pfcand_px", packed_cand->px());
    data.fillMulti<float>("pfcand_py", packed_cand->py());
    data.fillMulti<float>("pfcand_pz", packed_cand->pz());
    data.fillMulti<float>("pfcand_energy", packed_cand->energy());

    // *BELOW:*  Switchin packed_cand (ptr) -> SV (not ptr) where possible
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
    data.fillMulti<float>("pfcand_dzsig",
                          packed_cand->bestTrack() ? catchInfs(packed_cand->dz() / packed_cand->dzError()) : 0);
    data.fillMulti<float>("pfcand_dxy", catchInfs(packed_cand->dxy()));
    data.fillMulti<float>("pfcand_dxysig",
                          packed_cand->bestTrack() ? catchInfs(packed_cand->dxy() / packed_cand->dxyError()) : 0);

    if (packed_cand->bestTrack()) {
      const auto *trk = packed_cand->bestTrack();
      data.fillMulti<float>("pfcand_normchi2", catchInfs(trk->normalizedChi2()));
      data.fillMulti<float>("pfcand_quality", trk->qualityMask());
      data.fillMulti<float>("pfcand_nValidHits", trk->hitPattern().numberOfValidHits());
      data.fillMulti<float>("pfcand_nValidPixelHits", trk->hitPattern().numberOfValidPixelHits());

      // track covariance
      auto cov = [&](unsigned i, unsigned j) { return catchInfs(trk->covariance(i, j)); };
      data.fillMulti<float>("pfcand_dptdpt", cov(0, 0));
      data.fillMulti<float>("pfcand_detadeta", cov(1, 1));
      data.fillMulti<float>("pfcand_dphidphi", cov(2, 2));
      data.fillMulti<float>("pfcand_dxydxy", cov(3, 3));
      data.fillMulti<float>("pfcand_dzdz", cov(4, 4));
      data.fillMulti<float>("pfcand_dxydz", cov(3, 4));
      data.fillMulti<float>("pfcand_dphidxy", cov(2, 3));
      data.fillMulti<float>("pfcand_dlambdadz", cov(1, 4));
    } else {
      data.fillMulti<float>("pfcand_normchi2", 999);
      data.fillMulti<float>("pfcand_quality", 0);
      data.fillMulti<float>("pfcand_nValidHits", 0);
      data.fillMulti<float>("pfcand_nValidPixelHits", 0);

      data.fillMulti<float>("pfcand_dptdpt", 0);
      data.fillMulti<float>("pfcand_detadeta", 0);
      data.fillMulti<float>("pfcand_dphidphi", 0);
      data.fillMulti<float>("pfcand_dxydxy", 0);
      data.fillMulti<float>("pfcand_dzdz", 0);
      data.fillMulti<float>("pfcand_dxydz", 0);
      data.fillMulti<float>("pfcand_dphidxy", 0);
      data.fillMulti<float>("pfcand_dlambdadz", 0);
    }

    // build track info map
    TrackInfoBuilder trkinfo;
    // NEW:  Was (builder_, *packed_cand, jet, vertices->at(0));
    trkinfo.buildTrackInfo(builder_, *packed_cand, sv, vertices->at(0));  // jet -> sv

    data.fillMulti<float>("pfcand_btagMomentum", catchInfs(trkinfo.getTrackMomentum()));
    data.fillMulti<float>("pfcand_btagEta", catchInfs(trkinfo.getTrackEta()));
    data.fillMulti<float>("pfcand_btagEtaRel", catchInfs(trkinfo.getTrackEtaRel()));
    data.fillMulti<float>("pfcand_btagPtRel", catchInfs(trkinfo.getTrackPtRel()));
    data.fillMulti<float>("pfcand_btagPPar", catchInfs(trkinfo.getTrackPPar()));
    data.fillMulti<float>("pfcand_btagDeltaR", catchInfs(trkinfo.getTrackDeltaR()));
    data.fillMulti<float>("pfcand_btagPtRatio", catchInfs(trkinfo.getTrackPtRatio()));
    data.fillMulti<float>("pfcand_btagPParRatio", catchInfs(trkinfo.getTrackPParRatio()));
    data.fillMulti<float>("pfcand_btagSip2dVal", catchInfs(trkinfo.getTrackSip2dVal()));
    data.fillMulti<float>("pfcand_btagSip2dSig", catchInfs(trkinfo.getTrackSip2dSig()));
    data.fillMulti<float>("pfcand_btagSip3dVal", catchInfs(trkinfo.getTrackSip3dVal()));
    data.fillMulti<float>("pfcand_btagSip3dSig", catchInfs(trkinfo.getTrackSip3dSig()));
    data.fillMulti<float>("pfcand_btagJetDistVal", catchInfs(trkinfo.getTrackJetDistVal()));
    data.fillMulti<float>("pfcand_btagDecayLengthVal", catchInfs(trkinfo.getTrackDecayLengthVal()));
    data.fillMulti<float>("pfcand_btagDecayLengthSig", catchInfs(trkinfo.getTrackDecayLengthSig()));
  }

  return true;
}

} /* namespace deepntuples */
