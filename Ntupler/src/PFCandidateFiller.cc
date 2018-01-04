/*
 * PFCandidateFiller.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "DeepNTuples/Ntupler/interface/PFCandidateFiller.h"

namespace deepntuples {

void PFCandidateFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
}

void PFCandidateFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(svToken_, SVs);
}

void PFCandidateFiller::book() {

  data.add<int>("n_pfcands", 0);
  data.add<float>("npfcands", 0);

  // basic kinematics, valid for both charged and neutral
  data.addMulti<float>("pfcand_ptrel");
  data.addMulti<float>("pfcand_erel");
  data.addMulti<float>("pfcand_phirel");
  data.addMulti<float>("pfcand_etarel");
  data.addMulti<float>("pfcand_deltaR");
  data.addMulti<float>("pfcand_puppiw");
  data.addMulti<float>("pfcand_mass");

  data.addMulti<float>("pfcand_drminsv");
  data.addMulti<float>("pfcand_drsubjet1");
  data.addMulti<float>("pfcand_drsubjet2");

  data.addMulti<float>("pfcand_charge");
  data.addMulti<float>("pfcand_isMu");
  data.addMulti<float>("pfcand_isEl");
  data.addMulti<float>("pfcand_isGamma");
  data.addMulti<float>("pfcand_isChargedHad");
  data.addMulti<float>("pfcand_isNeutralHad");

  // for neutral
  data.addMulti<float>("pfcand_hcalFrac");

}

bool PFCandidateFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  const auto& jetConstituents = jet_helper.getJetConstituents();

  data.fill<int>("n_pfcands", jetConstituents.size());
  data.fill<float>("npfcands", jetConstituents.size());

  float etasign = jet.eta()>0 ? 1 : -1;

  for (const auto *pfcand : jetConstituents){

    // basic kinematics, valid for both charged and neutral
    data.fillMulti<float>("pfcand_ptrel", pfcand->pt()/jet.pt());
    data.fillMulti<float>("pfcand_erel", pfcand->energy()/jet.energy());
    data.fillMulti<float>("pfcand_phirel", reco::deltaPhi(*pfcand, jet));
    data.fillMulti<float>("pfcand_etarel", etasign * (pfcand->eta() - jet.eta()));
    data.fillMulti<float>("pfcand_deltaR", reco::deltaR(*pfcand, jet));
    data.fillMulti<float>("pfcand_puppiw", pfcand->puppiWeight());
    data.fillMulti<float>("pfcand_mass", pfcand->mass());

    double minDR = 999;
    for (const auto &sv : *SVs){
      double dr = reco::deltaR(*pfcand, sv);
      if (dr < minDR) minDR = dr;
    }
    data.fillMulti<float>("pfcand_drminsv", minDR==999 ? -1 : minDR);

    const auto& subjets = jet_helper.getSubJets();
    data.fillMulti<float>("pfcand_drsubjet1", subjets.size()>0 ? reco::deltaR(*pfcand, *subjets.at(0)) : -1);
    data.fillMulti<float>("pfcand_drsubjet2", subjets.size()>1 ? reco::deltaR(*pfcand, *subjets.at(1)) : -1);

    data.fillMulti<float>("pfcand_charge", pfcand->charge());
    data.fillMulti<float>("pfcand_isEl", std::abs(pfcand->pdgId())==11);
    data.fillMulti<float>("pfcand_isMu", std::abs(pfcand->pdgId())==13);
    data.fillMulti<float>("pfcand_isGamma", std::abs(pfcand->pdgId())==22);
    data.fillMulti<float>("pfcand_isChargedHad", std::abs(pfcand->pdgId())==211);
    data.fillMulti<float>("pfcand_isNeutralHad", std::abs(pfcand->pdgId())==130);

    // for neutral
    data.fillMulti<float>("pfcand_hcalFrac", pfcand->hcalFraction());

  }


  return true;
}

} /* namespace deepntuples */
