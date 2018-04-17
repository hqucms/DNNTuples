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
  fillElectronVars_ = iConfig.getUntrackedParameter<bool>("fillElectronVars", false);
  fillMuonVars_ = iConfig.getUntrackedParameter<bool>("fillMuonVars", false);
  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
  if (fillMuonVars_) muonToken_ = cc.consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));
  if (fillElectronVars_) {
    electronToken_ = cc.consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"));
    vetoIdToken_ = cc.consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("eleVetoIds"));
    looseIdToken_ = cc.consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("eleLooseIds"));
    mediumIdToken_ = cc.consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("eleMediumIds"));
    tightIdToken_ = cc.consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("eleTightIds"));
  }
}

void PFCompleteFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);
  if (fillMuonVars_) iEvent.getByToken(muonToken_, muons);
  if (fillElectronVars_) {
    iEvent.getByToken(electronToken_, electrons);
    iEvent.getByToken(vetoIdToken_, veto_id_decisions);
    iEvent.getByToken(looseIdToken_, loose_id_decisions);
    iEvent.getByToken(mediumIdToken_, medium_id_decisions);
    iEvent.getByToken(tightIdToken_, tight_id_decisions);
  }

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
//  data.addMulti<float>("part_btagJetDistSig"); // always gives 0

  // muon info
  data.addMulti<float>("part_muonIsLoose");
  data.addMulti<float>("part_muonIsMedium");
  data.addMulti<float>("part_muonIsTight");
  data.addMulti<float>("part_muonIsHighPt");
  data.addMulti<float>("part_muonSegmentCompatibility");
  data.addMulti<float>("part_muonNumberOfMatchedStations");

  // electron info
  data.addMulti<float>("part_electronVetoId");
  data.addMulti<float>("part_electronLooseId");
  data.addMulti<float>("part_electronMediumId");
  data.addMulti<float>("part_electronTightId");
  data.addMulti<float>("part_electronR9");
  data.addMulti<float>("part_electronSigmaIetaIeta");
  data.addMulti<float>("part_electronHadronicOverEm");
  data.addMulti<float>("part_electronPassConversionVeto");

}

bool PFCompleteFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  std::vector<const pat::PackedCandidate*> pfCands;
  std::unordered_map<const pat::PackedCandidate*, TrackInfoBuilder> trackInfoMap;
  std::unordered_map<const pat::PackedCandidate*, edm::Ptr<pat::Muon>> muonMap;
  std::unordered_map<const pat::PackedCandidate*, edm::Ptr<pat::Electron>> electronMap;
  for (const auto * cand : jet_helper.getJetConstituents()){
    pfCands.push_back(cand);
    // build track info map
    trackInfoMap[cand];
    trackInfoMap[cand].buildTrackInfo(builder_, *cand, jet, vertices->at(0));
    // map to pat muons
    muonMap[cand];
    if (fillMuonVars_){
      if (std::abs(cand->pdgId())==13){
        for (unsigned i=0; i<muons->size(); ++i){
          if (reco::deltaR2(*cand, muons->at(i)) < 1.e-6){
            muonMap[cand] = muons->ptrAt(i); break;
          }
        }
      }
    }
    // map to pat electrons
    electronMap[cand];
    if (fillElectronVars_){
      if (std::abs(cand->pdgId())==11){
        for (unsigned i=0; i<electrons->size(); ++i){
          if (reco::deltaR2(*cand, electrons->at(i)) < 1.e-6){
            electronMap[cand] = electrons->ptrAt(i); break;
          }
        }
      }
    }

  }

  // sort -- default is by pt
//  std::sort(pfCands.begin(), pfCands.end(), [&](const pat::PackedCandidate *p1, const pat::PackedCandidate *p2){
//    return trackInfoMap.at(p1).getTrackSip2dSig() > trackInfoMap.at(p2).getTrackSip2dSig();
//  });

  data.fill<int>("n_parts", pfCands.size());
  data.fill<float>("nparts", pfCands.size());

  float etasign = jet.eta()>0 ? 1 : -1;

  for (const auto *cand : pfCands){

    // basic kinematics, valid for both charged and neutral
    data.fillMulti<float>("part_ptrel", cand->pt()/jet.pt());
    data.fillMulti<float>("part_erel", cand->energy()/jet.energy());
    data.fillMulti<float>("part_phirel", reco::deltaPhi(*cand, jet));
    data.fillMulti<float>("part_etarel", etasign * (cand->eta() - jet.eta()));
    data.fillMulti<float>("part_deltaR", reco::deltaR(*cand, jet));
    data.fillMulti<float>("part_puppiw", cand->puppiWeight());
    data.fillMulti<float>("part_pt", cand->pt());
    data.fillMulti<float>("part_abseta", std::abs(cand->eta()));
    data.fillMulti<float>("part_mass", cand->mass());

    data.fillMulti<float>("part_ptrel_log", catchInfs(std::log(cand->pt()/jet.pt()), -99));
    data.fillMulti<float>("part_erel_log", catchInfs(std::log(cand->energy()/jet.energy()), -99));
    data.fillMulti<float>("part_pt_log", catchInfs(std::log(cand->pt()), -99));

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
//    data.fillMulti<float>("part_btagJetDistSig", trkinfo.getTrackJetDistSig());

    // muon info
    const auto mu = muonMap.at(cand);
    if (mu.isNonnull()){
      data.fillMulti<float>("part_muonIsLoose", mu->isLooseMuon());
      data.fillMulti<float>("part_muonIsMedium", mu->isMediumMuon());
      data.fillMulti<float>("part_muonIsTight", mu->isTightMuon(vertices->at(0)));
      data.fillMulti<float>("part_muonIsHighPt", mu->isHighPtMuon(vertices->at(0)));
      data.fillMulti<float>("part_muonSegmentCompatibility", mu->segmentCompatibility());
      data.fillMulti<float>("part_muonNumberOfMatchedStations", mu->numberOfMatchedStations());
    }else{
      data.fillMulti<float>("part_muonIsLoose", 0);
      data.fillMulti<float>("part_muonIsMedium", 0);
      data.fillMulti<float>("part_muonIsTight", 0);
      data.fillMulti<float>("part_muonIsHighPt", 0);
      data.fillMulti<float>("part_muonSegmentCompatibility", 0);
      data.fillMulti<float>("part_muonNumberOfMatchedStations", 0);
    }

    // electron info
    const auto ele = electronMap.at(cand);
    if (ele.isNonnull()){
      data.fillMulti<float>("part_electronVetoId", (*veto_id_decisions)[ele]);
      data.fillMulti<float>("part_electronLooseId", (*loose_id_decisions)[ele]);
      data.fillMulti<float>("part_electronMediumId", (*medium_id_decisions)[ele]);
      data.fillMulti<float>("part_electronTightId", (*tight_id_decisions)[ele]);
      data.fillMulti<float>("part_electronR9", ele->full5x5_r9());
      data.fillMulti<float>("part_electronSigmaIetaIeta", ele->full5x5_sigmaIetaIeta());
      data.fillMulti<float>("part_electronHadronicOverEm", ele->hadronicOverEm());
      data.fillMulti<float>("part_electronPassConversionVeto", ele->passConversionVeto());
    }else{
      data.fillMulti<float>("part_electronVetoId", 0);
      data.fillMulti<float>("part_electronLooseId", 0);
      data.fillMulti<float>("part_electronMediumId", 0);
      data.fillMulti<float>("part_electronTightId", 0);
      data.fillMulti<float>("part_electronR9", 0);
      data.fillMulti<float>("part_electronSigmaIetaIeta", 0);
      data.fillMulti<float>("part_electronHadronicOverEm", 0);
      data.fillMulti<float>("part_electronPassConversionVeto", 0);
    }

  }


  return true;
}

} /* namespace deepntuples */
