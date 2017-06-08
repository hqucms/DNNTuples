/*
 * FatJetInfoFiller.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "DeepNTuples/NtupleAK8/interface/FatJetInfoFiller.h"

namespace deepntuples {

void FatJetInfoFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  genParticlesToken_ = cc.consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  fjTagInfoName = iConfig.getParameter<std::string>("fjTagInfoName");
  for (const auto &flv : iConfig.getUntrackedParameter<std::vector<unsigned>>("fjKeepFlavors", {})){
    keepFlavors_.push_back(static_cast<FatJetMatching::FatJetFlavor>(flv));
  }
}

void FatJetInfoFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(genParticlesToken_, genParticlesHandle);
}

void FatJetInfoFiller::book() {
  // truth labels
  data.add<int>("fj_label", 0);
  data.add<int>("fj_isLight", 0);
  data.add<int>("fj_isTop", 0);
  data.add<int>("fj_isW", 0);
  data.add<int>("fj_isZ", 0);
  data.add<int>("fj_isH", 0);

  // JMAR label
  data.add<int>("fj_labelJMAR", 0);

  // gen-matched particle (top/W/etc.)
  data.add<float>("fj_gen_pt", 0);
  data.add<float>("fj_gen_eta", 0);

  // fatjet kinematics
  data.add<float>("fj_pt", 0);
  data.add<float>("fj_eta", 0);
  data.add<float>("fj_phi", 0);
  data.add<float>("fj_mass", 0);

  // substructure
  data.add<float>("fj_tau1", 0);
  data.add<float>("fj_tau2", 0);
  data.add<float>("fj_tau3", 0);
  data.add<float>("fj_tau21", 0);
  data.add<float>("fj_tau32", 0);

  // soft drop
  data.add<float>("fj_sdmass", 0);

  // subjets: soft drop gives up to 2 subjets
  data.add<float>("fj_n_sdsubjets", 0);

  data.add<float>("fj_sdsj1_pt", 0);
  data.add<float>("fj_sdsj1_eta", 0);
  data.add<float>("fj_sdsj1_phi", 0);
  data.add<float>("fj_sdsj1_mass", 0);
  data.add<float>("fj_sdsj1_csv", 0);
  data.add<float>("fj_sdsj1_ptD", 0);
  data.add<float>("fj_sdsj1_axis1", 0);
  data.add<float>("fj_sdsj1_axis2", 0);
  data.add<float>("fj_sdsj1_mult", 0);

  data.add<float>("fj_sdsj2_pt", 0);
  data.add<float>("fj_sdsj2_eta", 0);
  data.add<float>("fj_sdsj2_phi", 0);
  data.add<float>("fj_sdsj2_mass", 0);
  data.add<float>("fj_sdsj2_csv", 0);
  data.add<float>("fj_sdsj2_ptD", 0);
  data.add<float>("fj_sdsj2_axis1", 0);
  data.add<float>("fj_sdsj2_axis2", 0);
  data.add<float>("fj_sdsj2_mult", 0);

  // some variables used in a baseline tagger
  data.add<float>("fj_ptDR", 0);
  data.add<float>("fj_relptdiff", 0);
  data.add<float>("fj_sdn2", 0);


  //double-b
  data.add<float>("fj_doubleb", 0);

  //flavor info
  data.add<int>("fj_isBB", 0);
  data.add<int>("fj_isNonBB", 0);
  data.add<int>("fj_nbHadrons", 0);
  data.add<int>("fj_ncHadrons", 0);

  //double-b inputs
  data.add<float>("fj_z_ratio", 0);
  data.add<float>("fj_trackSipdSig_3", 0);
  data.add<float>("fj_trackSipdSig_2", 0);
  data.add<float>("fj_trackSipdSig_1", 0);
  data.add<float>("fj_trackSipdSig_0", 0);
  data.add<float>("fj_trackSipdSig_1_0", 0);
  data.add<float>("fj_trackSipdSig_0_0", 0);
  data.add<float>("fj_trackSipdSig_1_1", 0);
  data.add<float>("fj_trackSipdSig_0_1", 0);
  data.add<float>("fj_trackSip2dSigAboveCharm_0", 0);
  data.add<float>("fj_trackSip2dSigAboveBottom_0", 0);
  data.add<float>("fj_trackSip2dSigAboveBottom_1", 0);
  data.add<float>("fj_tau1_trackEtaRel_0", 0);
  data.add<float>("fj_tau1_trackEtaRel_1", 0);
  data.add<float>("fj_tau1_trackEtaRel_2", 0);
  data.add<float>("fj_tau0_trackEtaRel_0", 0);
  data.add<float>("fj_tau0_trackEtaRel_1", 0);
  data.add<float>("fj_tau0_trackEtaRel_2", 0);
  data.add<float>("fj_tau_vertexMass_0", 0);
  data.add<float>("fj_tau_vertexEnergyRatio_0", 0);
  data.add<float>("fj_tau_vertexDeltaR_0", 0);
  data.add<float>("fj_tau_flightDistance2dSig_0", 0);
  data.add<float>("fj_tau_vertexMass_1", 0);
  data.add<float>("fj_tau_vertexEnergyRatio_1", 0);
  data.add<float>("fj_tau_flightDistance2dSig_1", 0);
  data.add<float>("fj_jetNTracks", 0);
  data.add<float>("fj_nSV", 0);


}

bool FatJetInfoFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  auto genmatch = fjmatch_.flavor(&jet, *genParticlesHandle);

  // get only specific types of jet from a given sample
  if (keepFlavors_.size() && std::find(keepFlavors_.begin(), keepFlavors_.end(), genmatch.first)==keepFlavors_.end()){
    return false;
  }

  // truth labels
  data.fill<int>("fj_label", genmatch.first);
  data.fill<int>("fj_isLight", genmatch.first == FatJetMatching::Default);
  data.fill<int>("fj_isTop", genmatch.first == FatJetMatching::Top);
  data.fill<int>("fj_isW", genmatch.first == FatJetMatching::W);
  data.fill<int>("fj_isZ", genmatch.first == FatJetMatching::Z);
  data.fill<int>("fj_isH", genmatch.first == FatJetMatching::H);

  // JMAR label
  data.fill<int>("fj_labelJMAR", fjmatch_.flavorJMAR(&jet, *genParticlesHandle, 0.6).first);

  // gen-matched particle (top/W/etc.)
  data.fill<float>("fj_gen_pt", genmatch.second ? genmatch.second->pt() : -999);
  data.fill<float>("fj_gen_eta", genmatch.second ? genmatch.second->eta() : -999);

  // fatjet kinematics
  data.fill<float>("fj_pt", jet.pt());
  data.fill<float>("fj_eta", jet.eta());
  data.fill<float>("fj_phi", jet.phi());
  data.fill<float>("fj_mass", jet.mass());

  // substructure
  float tau1 = jet.userFloat("NjettinessAK8:tau1");
  float tau2 = jet.userFloat("NjettinessAK8:tau2");
  float tau3 = jet.userFloat("NjettinessAK8:tau3");
  data.fill<float>("fj_tau1", tau1);
  data.fill<float>("fj_tau2", tau2);
  data.fill<float>("fj_tau3", tau3);
  data.fill<float>("fj_tau21", tau1 > 0 ? tau2/tau1 : 1.01);
  data.fill<float>("fj_tau32", tau2 > 0 ? tau3/tau2 : 1.01);

  // soft drop
  data.fill<float>("fj_sdmass", jet.userFloat("ak8PFJetsCHSSoftDropMass"));

  // subjets: soft drop gives up to 2 subjets
  const auto& subjets = jet_helper.getSubJets();
  data.fill<float>("fj_n_sdsubjets", subjets.size());

  if (subjets.size() > 0){
    const auto &sj1 = subjets.at(0);
    JetHelper jh1(sj1);
    data.fill<float>("fj_sdsj1_pt", sj1->pt());
    data.fill<float>("fj_sdsj1_eta", sj1->eta());
    data.fill<float>("fj_sdsj1_phi", sj1->phi());
    data.fill<float>("fj_sdsj1_mass", sj1->mass());
    data.fill<float>("fj_sdsj1_csv", sj1->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    data.fill<float>("fj_sdsj1_ptD", jh1.ptD());
    data.fill<float>("fj_sdsj1_axis1", jh1.axis1());
    data.fill<float>("fj_sdsj1_axis2", jh1.axis2());
    data.fill<float>("fj_sdsj1_mult", jh1.mult());

    if (subjets.size() > 1){
      const auto &sj2 = subjets.at(1);
      JetHelper jh2(sj2);
      data.fill<float>("fj_sdsj2_pt", sj2->pt());
      data.fill<float>("fj_sdsj2_eta", sj2->eta());
      data.fill<float>("fj_sdsj2_phi", sj2->phi());
      data.fill<float>("fj_sdsj2_mass", sj2->mass());
      data.fill<float>("fj_sdsj2_csv", sj2->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      data.fill<float>("fj_sdsj2_ptD", jh2.ptD());
      data.fill<float>("fj_sdsj2_axis1", jh2.axis1());
      data.fill<float>("fj_sdsj2_axis2", jh2.axis2());
      data.fill<float>("fj_sdsj2_mult", jh2.mult());

      // some variables used in a baseline tagger
      float deltaR = reco::deltaR(*sj1, *sj2);
      float var_sd_0 = sj2->pt()/(sj1->pt()+sj2->pt());
      data.fill<float>("fj_ptDR", jet.pt() * deltaR);
      data.fill<float>("fj_relptdiff", std::abs(sj1->pt()-sj2->pt()) / jet.pt());
      data.fill<float>("fj_sdn2", var_sd_0/std::pow(deltaR,-2));
    }
  }

  // --------
  // double-b

  const auto *bdsvTagInfo = jet.tagInfoBoostedDoubleSV(fjTagInfoName);
  assert(bdsvTagInfo);
  const auto &vars = bdsvTagInfo->taggingVariables();

  data.fill<float>("fj_doubleb", jet.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));

  //flavor info
  data.fill<int>("fj_isBB", jet.jetFlavourInfo().getbHadrons().size() >= 2);
  data.fill<int>("fj_isNonBB", jet.jetFlavourInfo().getbHadrons().size() < 2);
  data.fill<int>("fj_nbHadrons", jet.jetFlavourInfo().getbHadrons().size());
  data.fill<int>("fj_ncHadrons", jet.jetFlavourInfo().getcHadrons().size());

  //double-b inputs
  data.fill<float>("fj_z_ratio", vars.get(reco::btau::z_ratio));
  data.fill<float>("fj_trackSipdSig_3", vars.get(reco::btau::trackSip3dSig_3));
  data.fill<float>("fj_trackSipdSig_2", vars.get(reco::btau::trackSip3dSig_2));
  data.fill<float>("fj_trackSipdSig_1", vars.get(reco::btau::trackSip3dSig_1));
  data.fill<float>("fj_trackSipdSig_0", vars.get(reco::btau::trackSip3dSig_0));
  data.fill<float>("fj_trackSipdSig_1_0", vars.get(reco::btau::tau2_trackSip3dSig_0));
  data.fill<float>("fj_trackSipdSig_0_0", vars.get(reco::btau::tau1_trackSip3dSig_0));
  data.fill<float>("fj_trackSipdSig_1_1", vars.get(reco::btau::tau2_trackSip3dSig_1));
  data.fill<float>("fj_trackSipdSig_0_1", vars.get(reco::btau::tau1_trackSip3dSig_1));
  data.fill<float>("fj_trackSip2dSigAboveCharm_0", vars.get(reco::btau::trackSip2dSigAboveCharm));
  data.fill<float>("fj_trackSip2dSigAboveBottom_0", vars.get(reco::btau::trackSip2dSigAboveBottom_0));
  data.fill<float>("fj_trackSip2dSigAboveBottom_1", vars.get(reco::btau::trackSip2dSigAboveBottom_1));
  data.fill<float>("fj_tau1_trackEtaRel_0", vars.get(reco::btau::tau2_trackEtaRel_0));
  data.fill<float>("fj_tau1_trackEtaRel_1", vars.get(reco::btau::tau2_trackEtaRel_1));
  data.fill<float>("fj_tau1_trackEtaRel_2", vars.get(reco::btau::tau2_trackEtaRel_2));
  data.fill<float>("fj_tau0_trackEtaRel_0", vars.get(reco::btau::tau1_trackEtaRel_0));
  data.fill<float>("fj_tau0_trackEtaRel_1", vars.get(reco::btau::tau1_trackEtaRel_1));
  data.fill<float>("fj_tau0_trackEtaRel_2", vars.get(reco::btau::tau1_trackEtaRel_2));
  data.fill<float>("fj_tau_vertexMass_0", vars.get(reco::btau::tau1_vertexMass));
  data.fill<float>("fj_tau_vertexEnergyRatio_0", vars.get(reco::btau::tau1_vertexEnergyRatio));
  data.fill<float>("fj_tau_vertexDeltaR_0", vars.get(reco::btau::tau1_vertexDeltaR));
  data.fill<float>("fj_tau_flightDistance2dSig_0", vars.get(reco::btau::tau1_flightDistance2dSig));
  data.fill<float>("fj_tau_vertexMass_1", vars.get(reco::btau::tau2_vertexMass));
  data.fill<float>("fj_tau_vertexEnergyRatio_1", vars.get(reco::btau::tau2_vertexEnergyRatio));
  data.fill<float>("fj_tau_flightDistance2dSig_1", vars.get(reco::btau::tau2_flightDistance2dSig));
  data.fill<float>("fj_jetNTracks", vars.get(reco::btau::jetNTracks));
  data.fill<float>("fj_nSV", vars.get(reco::btau::jetNSecondaryVertices));

  return true;
}

} /* namespace deepntuples */

