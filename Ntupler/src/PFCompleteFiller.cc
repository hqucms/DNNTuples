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

  data.add<int>("n_parts", 0);
  data.add<float>("nparts", 0);

  // basic kinematics
  data.addMulti<float>("part_ptrel");
  data.addMulti<float>("part_erel");
  data.addMulti<float>("part_phirel");
  data.addMulti<float>("part_etarel");
  data.addMulti<float>("part_deltaR");
  data.addMulti<float>("part_puppiw");
  data.addMulti<float>("part_pt");
  data.addMulti<float>("part_abseta");
  data.addMulti<float>("part_mass");

  data.addMulti<float>("part_ptrel_log");
  data.addMulti<float>("part_erel_log");
  data.addMulti<float>("part_pt_log");

  data.addMulti<float>("part_drminsv");
  data.addMulti<float>("part_drsubjet1");
  data.addMulti<float>("part_drsubjet2");

  data.addMulti<float>("part_charge");
  data.addMulti<float>("part_isMu");
  data.addMulti<float>("part_isEl");
  data.addMulti<float>("part_isChargedHad");
  data.addMulti<float>("part_isGamma");
  data.addMulti<float>("part_isNeutralHad");

  // for neutral
  data.addMulti<float>("part_hcalFrac");

  // for charged
  data.addMulti<float>("part_VTX_ass");
  data.addMulti<float>("part_fromPV");
  data.addMulti<float>("part_lostInnerHits");
  data.addMulti<float>("part_trackHighPurity");

  // impact parameters
  data.addMulti<float>("part_dz");
  data.addMulti<float>("part_dzsig");
  data.addMulti<float>("part_dxy");
  data.addMulti<float>("part_dxysig");

  // track quality
  data.addMulti<float>("part_normchi2");
  data.addMulti<float>("part_quality");

  // track covariance
  data.addMulti<float>("part_dptdpt");
  data.addMulti<float>("part_detadeta");
  data.addMulti<float>("part_dphidphi");
  data.addMulti<float>("part_dxydxy");
  data.addMulti<float>("part_dzdz");
  data.addMulti<float>("part_dxydz");
  data.addMulti<float>("part_dphidxy");
  data.addMulti<float>("part_dlambdadz");

  // track btag info
  data.addMulti<float>("part_btagMomentum");
  data.addMulti<float>("part_btagEta");
  data.addMulti<float>("part_btagEtaRel");
  data.addMulti<float>("part_btagPtRel");
  data.addMulti<float>("part_btagPPar");
  data.addMulti<float>("part_btagDeltaR");
  data.addMulti<float>("part_btagPtRatio");
  data.addMulti<float>("part_btagPParRatio");
  data.addMulti<float>("part_btagSip2dVal");
  data.addMulti<float>("part_btagSip2dSig");
  data.addMulti<float>("part_btagSip3dVal");
  data.addMulti<float>("part_btagSip3dSig");
  data.addMulti<float>("part_btagJetDistVal");
  data.addMulti<float>("part_btagJetDistSig"); // always gives 0?

}

bool PFCompleteFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  std::vector<const pat::PackedCandidate*> pfCands;
  std::unordered_map<const pat::PackedCandidate*, TrackInfoBuilder> trackInfoMap;
  for (const auto * cand : jet_helper.getJetConstituents()){
    pfCands.push_back(cand);
    // build track info map
    trackInfoMap[cand];
    trackInfoMap[cand].buildTrackInfo(builder_, *cand, jet, vertices->at(0));
  }

  // sort -- default is by pt
  data.fill<int>("n_parts", pfCands.size());
  data.fill<float>("nparts", pfCands.size());

  float etasign = jet.eta()>0 ? 1 : -1;

  for (const auto *cand : pfCands){

    auto puppiP4 = cand->p4();

    // for jets stored in MiniAOD, need to rescale p4 by the puppi weight
    if (!jet_helper.hasPuppiWeightedDaughters()) puppiP4 *= cand->puppiWeight();

    // basic kinematics, valid for both charged and neutral
    data.fillMulti<float>("part_pt", puppiP4.pt());
    data.fillMulti<float>("part_ptrel", puppiP4.pt()/jet.pt());
    data.fillMulti<float>("part_erel", puppiP4.energy()/jet.energy());
    data.fillMulti<float>("part_phirel", reco::deltaPhi(puppiP4, jet));
    data.fillMulti<float>("part_etarel", etasign * (puppiP4.eta() - jet.eta()));
    data.fillMulti<float>("part_deltaR", reco::deltaR(puppiP4, jet));
    data.fillMulti<float>("part_abseta", std::abs(puppiP4.eta()));
    data.fillMulti<float>("part_mass", puppiP4.mass());

    data.fillMulti<float>("part_ptrel_log", catchInfs(std::log(puppiP4.pt()/jet.pt()), -99));
    data.fillMulti<float>("part_erel_log", catchInfs(std::log(puppiP4.energy()/jet.energy()), -99));
    data.fillMulti<float>("part_pt_log", catchInfs(std::log(puppiP4.pt()), -99));

    data.fillMulti<float>("part_puppiw", cand->puppiWeight());

    double minDR = 999;
    for (const auto &sv : *SVs){
      double dr = reco::deltaR(*cand, sv);
      if (dr < minDR) minDR = dr;
    }
    data.fillMulti<float>("part_drminsv", minDR==999 ? -1 : minDR);

    const auto& subjets = jet_helper.getSubJets();
    data.fillMulti<float>("part_drsubjet1", subjets.size()>0 ? reco::deltaR(*cand, *subjets.at(0)) : -1);
    data.fillMulti<float>("part_drsubjet2", subjets.size()>1 ? reco::deltaR(*cand, *subjets.at(1)) : -1);

    data.fillMulti<float>("part_charge", cand->charge());
    data.fillMulti<float>("part_isEl", std::abs(cand->pdgId())==11);
    data.fillMulti<float>("part_isMu", std::abs(cand->pdgId())==13);
    data.fillMulti<float>("part_isChargedHad", std::abs(cand->pdgId())==211);
    data.fillMulti<float>("part_isGamma", std::abs(cand->pdgId())==22);
    data.fillMulti<float>("part_isNeutralHad", std::abs(cand->pdgId())==130);

    // for neutral
    data.fillMulti<float>("part_hcalFrac", cand->hcalFraction());

    // for charged
    data.fillMulti<float>("part_VTX_ass", cand->pvAssociationQuality());
    data.fillMulti<float>("part_fromPV", cand->fromPV());
    data.fillMulti<float>("part_lostInnerHits", cand->lostInnerHits());
    data.fillMulti<float>("part_trackHighPurity", cand->trackHighPurity());

    // impact parameters
    data.fillMulti<float>("part_dz", catchInfs(cand->dz()));
    data.fillMulti<float>("part_dzsig", cand->bestTrack() ? catchInfs(cand->dz()/cand->dzError()) : 0);
    data.fillMulti<float>("part_dxy", catchInfs(cand->dxy()));
    data.fillMulti<float>("part_dxysig", cand->bestTrack() ? catchInfs(cand->dxy()/cand->dxyError()) : 0);

    if (cand->bestTrack()){
      const auto *trk = cand->bestTrack();
      data.fillMulti<float>("part_normchi2", catchInfs(trk->normalizedChi2()));
      data.fillMulti<float>("part_quality", trk->qualityMask());

      // track covariance
      auto cov = [&](unsigned i, unsigned j) {
        return catchInfs(trk->covariance(i, j));
      };
      data.fillMulti<float>("part_dptdpt", cov(0,0));
      data.fillMulti<float>("part_detadeta", cov(1,1));
      data.fillMulti<float>("part_dphidphi", cov(2,2));
      data.fillMulti<float>("part_dxydxy", cov(3,3));
      data.fillMulti<float>("part_dzdz", cov(4,4));
      data.fillMulti<float>("part_dxydz", cov(3,4));
      data.fillMulti<float>("part_dphidxy", cov(2,3));
      data.fillMulti<float>("part_dlambdadz", cov(1,4));
    }else{
      data.fillMulti<float>("part_normchi2", 999);
      data.fillMulti<float>("part_quality", 0);

      data.fillMulti<float>("part_dptdpt", 0);
      data.fillMulti<float>("part_detadeta", 0);
      data.fillMulti<float>("part_dphidphi", 0);
      data.fillMulti<float>("part_dxydxy", 0);
      data.fillMulti<float>("part_dzdz", 0);
      data.fillMulti<float>("part_dxydz", 0);
      data.fillMulti<float>("part_dphidxy", 0);
      data.fillMulti<float>("part_dlambdadz", 0);
    }

    const auto &trkinfo = trackInfoMap.at(cand);
    data.fillMulti<float>("part_btagMomentum", trkinfo.getTrackMomentum());
    data.fillMulti<float>("part_btagEta", trkinfo.getTrackEta());
    data.fillMulti<float>("part_btagEtaRel", trkinfo.getTrackEtaRel());
    data.fillMulti<float>("part_btagPtRel", trkinfo.getTrackPtRel());
    data.fillMulti<float>("part_btagPPar", trkinfo.getTrackPPar());
    data.fillMulti<float>("part_btagDeltaR", trkinfo.getTrackDeltaR());
    data.fillMulti<float>("part_btagPtRatio", trkinfo.getTrackPtRatio());
    data.fillMulti<float>("part_btagPParRatio", trkinfo.getTrackPParRatio());
    data.fillMulti<float>("part_btagSip2dVal", trkinfo.getTrackSip2dVal());
    data.fillMulti<float>("part_btagSip2dSig", trkinfo.getTrackSip2dSig());
    data.fillMulti<float>("part_btagSip3dVal", trkinfo.getTrackSip3dVal());
    data.fillMulti<float>("part_btagSip3dSig", trkinfo.getTrackSip3dSig());
    data.fillMulti<float>("part_btagJetDistVal", trkinfo.getTrackJetDistVal());
    data.fillMulti<float>("part_btagJetDistSig", trkinfo.getTrackJetDistSig());


  }


  return true;
}

} /* namespace deepntuples */
