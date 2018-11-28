/*
 * FatJetInfoFiller.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "DeepNTuples/Ntupler/interface/FatJetInfoFiller.h"
#include <string>
#include <algorithm>

namespace deepntuples {

void FatJetInfoFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  genParticlesToken_ = cc.consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  fjTagInfoName = iConfig.getParameter<std::string>("fjTagInfoName");
  useReclusteredJets_ = iConfig.getParameter<bool>("useReclusteredJets");
  isQCDSample_ = iConfig.getUntrackedParameter<bool>("isQCDSample", false);
  isTrainSample_ = iConfig.getUntrackedParameter<bool>("isTrainSample", false);
  fjName = iConfig.getParameter<std::string>("jetType") + std::to_string(int(10*jetR_));
}

void FatJetInfoFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(genParticlesToken_, genParticlesHandle);
}

void FatJetInfoFiller::book() {
  // truth labels
  data.add<int>("fj_label", 0);
  data.add<int>("fj_isTop", 0);
  data.add<int>("fj_isW", 0);
  data.add<int>("fj_isZ", 0);
  data.add<int>("fj_isH", 0);
  data.add<int>("fj_isQCD", 0);

  data.add<int>("label_Top_bcq", 0);
  data.add<int>("label_Top_bqq", 0);
  data.add<int>("label_Top_bc",  0);
  data.add<int>("label_Top_bq",  0);

  data.add<int>("label_W_cq",    0);
  data.add<int>("label_W_qq",    0);

  data.add<int>("label_Z_bb",    0);
  data.add<int>("label_Z_cc",    0);
  data.add<int>("label_Z_qq",    0);

  data.add<int>("label_H_bb",    0);
  data.add<int>("label_H_cc",    0);
  data.add<int>("label_H_qqqq",  0);
  data.add<int>("label_H_tautau",0);

  data.add<int>("label_QCD_bb",  0);
  data.add<int>("label_QCD_cc",  0);
  data.add<int>("label_QCD_b",   0);
  data.add<int>("label_QCD_c",   0);
  data.add<int>("label_QCD_others", 0);

  data.add<int>("sample_isQCD", 0);

  data.add<int>("sample_useReclusteredJets", useReclusteredJets_);
  data.add<int>("sample_hasPuppiWeightedDaughters", 0);

//  // legacy labels
//  data.add<int>("fj_labelLegacy", 0);

  // JMAR label
  data.add<int>("fj_labelJMAR", 0);
  data.add<float>("fjJMAR_gen_pt", 0);
  data.add<float>("fjJMAR_gen_eta", 0);
  data.add<float>("fjJMAR_gen_phi", 0);
  data.add<int>("fjJMAR_gen_pdgid", 0);

  // gen-matched particle (top/W/etc.)
  data.add<float>("fj_gen_pt", 0);
  data.add<float>("fj_gen_eta", 0);

  // --- jet energy/mass regression ---
  data.add<float>("fj_genjet_pt", 0);
  data.add<float>("fj_genOverReco_pt", 1); // default to 1 if not gen-matched
  data.add<float>("fj_genOverReco_pt_null", 0); // default to 0 if not gen-matched
  data.add<float>("fj_genjet_mass", 0);
  data.add<float>("fj_genOverReco_mass", 1); // default to 1 if not gen-matched
  data.add<float>("fj_genOverReco_mass_null", 0); // default to 0 if not gen-matched
  data.add<float>("fj_genjet_sdmass", 0);
  data.add<float>("fj_genjet_sdmass_sqrt", 0);
  data.add<float>("fj_genOverReco_sdmass", 1); // default to 1 if not gen-matched
  data.add<float>("fj_genOverReco_sdmass_null", 0); // default to 0 if not gen-matched
  // ----------------------------------

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
  data.add<float>("fj_sdmass_fromsubjets", 0);
  data.add<float>("fj_corrsdmass", 0);

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

  data.fill<int>("sample_hasPuppiWeightedDaughters", jet_helper.hasPuppiWeightedDaughters());

  // JMAR label
  {
    auto jmar = fjmatch_.flavorJMAR(&jet, *genParticlesHandle, 0.6);
    data.fill<int>("fj_labelJMAR", jmar.first);
    data.fill<float>("fjJMAR_gen_pt", jmar.second ? jmar.second->pt() : -999);
    data.fill<float>("fjJMAR_gen_eta", jmar.second ? jmar.second->eta() : -999);
    data.fill<float>("fjJMAR_gen_phi", jmar.second ? jmar.second->phi() : -999);
    data.fill<int>("fjJMAR_gen_pdgid", jmar.second ? jmar.second->pdgId() : -999);
  }

  // ----------------------------------------------------------------
//  auto fjlabel = fjmatch_.flavorLabel(&jet, *genParticlesHandle, 0.6);
  auto fjlabel = fjmatch_.flavorLabel(&jet, *genParticlesHandle, jetR_);

  data.fill<int>("fj_label", fjlabel.first);

  data.fill<int>("fj_isTop", fjlabel.first >= FatJetMatching::Top_all && fjlabel.first < FatJetMatching::W_all);
  data.fill<int>("fj_isW",   fjlabel.first >= FatJetMatching::W_all && fjlabel.first < FatJetMatching::Z_all);
  data.fill<int>("fj_isZ",   fjlabel.first >= FatJetMatching::Z_all && fjlabel.first < FatJetMatching::H_all);
  data.fill<int>("fj_isH",   fjlabel.first >= FatJetMatching::H_all && fjlabel.first < FatJetMatching::QCD_all);
  data.fill<int>("fj_isQCD", fjlabel.first >= FatJetMatching::QCD_all);

  // veto unmatched jets in signal samples for training
  if (isTrainSample_ && !isQCDSample_ && fjlabel.first >= FatJetMatching::QCD_all)
    return false;

  data.fill<int>("label_Top_bcq", fjlabel.first == FatJetMatching::Top_bcq);
  data.fill<int>("label_Top_bqq", fjlabel.first == FatJetMatching::Top_bqq);
  data.fill<int>("label_Top_bc",  fjlabel.first == FatJetMatching::Top_bc);
  data.fill<int>("label_Top_bq",  fjlabel.first == FatJetMatching::Top_bq);

  data.fill<int>("label_W_cq",    fjlabel.first == FatJetMatching::W_cq);
  data.fill<int>("label_W_qq",    fjlabel.first == FatJetMatching::W_qq);

  data.fill<int>("label_Z_bb",    fjlabel.first == FatJetMatching::Z_bb);
  data.fill<int>("label_Z_cc",    fjlabel.first == FatJetMatching::Z_cc);
  data.fill<int>("label_Z_qq",    fjlabel.first == FatJetMatching::Z_qq);

  data.fill<int>("label_H_bb",    fjlabel.first == FatJetMatching::H_bb);
  data.fill<int>("label_H_cc",    fjlabel.first == FatJetMatching::H_cc);
  data.fill<int>("label_H_qqqq",  fjlabel.first == FatJetMatching::H_qqqq);
  data.fill<int>("label_H_tautau",fjlabel.first == FatJetMatching::H_tautau);

  data.fill<int>("label_QCD_bb",  fjlabel.first == FatJetMatching::QCD_bb);
  data.fill<int>("label_QCD_cc",  fjlabel.first == FatJetMatching::QCD_cc);
  data.fill<int>("label_QCD_b",   fjlabel.first == FatJetMatching::QCD_b);
  data.fill<int>("label_QCD_c",   fjlabel.first == FatJetMatching::QCD_c);
  data.fill<int>("label_QCD_others", fjlabel.first == FatJetMatching::QCD_others);

  data.fill<int>("sample_isQCD",  isQCDSample_);


  // gen-matched particle (top/W/etc.)
  data.fill<float>("fj_gen_pt", fjlabel.second ? fjlabel.second->pt() : -999);
  data.fill<float>("fj_gen_eta", fjlabel.second ? fjlabel.second->eta() : -999);

  // ----------------------------------

  // fatjet kinematics
  data.fill<float>("fj_pt", jet.pt());
  data.fill<float>("fj_eta", jet.eta());
  data.fill<float>("fj_phi", jet.phi());
  data.fill<float>("fj_mass", jet.mass());

  // substructure
  float tau1 = jet.userFloat("Njettiness" + fjName +"Puppi:tau1");
  float tau2 = jet.userFloat("Njettiness" + fjName +"Puppi:tau2");
  float tau3 = jet.userFloat("Njettiness" + fjName +"Puppi:tau3");
  data.fill<float>("fj_tau1", tau1);
  data.fill<float>("fj_tau2", tau2);
  data.fill<float>("fj_tau3", tau3);
  data.fill<float>("fj_tau21", tau1 > 0 ? tau2/tau1 : 1.01);
  data.fill<float>("fj_tau32", tau2 > 0 ? tau3/tau2 : 1.01);

  // soft drop
  std::string name(fjName);
  std::transform(name.begin(), name.end(), name.begin(), ::tolower);
  auto msd_uncorr = jet.userFloat(name +"PFJetsPuppiSoftDropMass");
  data.fill<float>("fj_sdmass", msd_uncorr);
  data.fill<float>("fj_sdmass_fromsubjets", jet.groomedMass());

  // subjets: soft drop gives up to 2 subjets
  const auto& subjets = jet_helper.getSubJets();

  data.fill<float>("fj_n_sdsubjets", subjets.size());
  data.fill<float>("fj_corrsdmass", jet_helper.getCorrectedPuppiSoftDropMass(subjets).second);

  if (subjets.size() > 0){
    const auto &sj1 = subjets.at(0);
    data.fill<float>("fj_sdsj1_pt", sj1->pt());
    data.fill<float>("fj_sdsj1_eta", sj1->eta());
    data.fill<float>("fj_sdsj1_phi", sj1->phi());
    data.fill<float>("fj_sdsj1_mass", sj1->mass());
    data.fill<float>("fj_sdsj1_csv", sj1->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

    if (subjets.size() > 1){
      const auto &sj2 = subjets.at(1);
      data.fill<float>("fj_sdsj2_pt", sj2->pt());
      data.fill<float>("fj_sdsj2_eta", sj2->eta());
      data.fill<float>("fj_sdsj2_phi", sj2->phi());
      data.fill<float>("fj_sdsj2_mass", sj2->mass());
      data.fill<float>("fj_sdsj2_csv", sj2->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

      // some variables used in a baseline tagger
      float deltaR = reco::deltaR(*sj1, *sj2);
      float var_sd_0 = sj2->pt()/(sj1->pt()+sj2->pt());
      data.fill<float>("fj_ptDR", jet.pt() * deltaR);
      data.fill<float>("fj_relptdiff", std::abs(sj1->pt()-sj2->pt()) / jet.pt());
      data.fill<float>("fj_sdn2", var_sd_0/std::pow(deltaR,-2));
    }
  }

  // ----------------------------------------------------------------

  // --- jet energy/mass regression ---
  const auto *genjet = jet_helper.genjetWithNu();
  if (genjet){
    // jet here points to the uncorrected jet
    data.fill<float>("fj_genjet_pt", genjet->pt());
    data.fill<float>("fj_genOverReco_pt", catchInfs(genjet->pt() / jet.pt(), 1));
    data.fill<float>("fj_genOverReco_pt_null", catchInfs(genjet->pt() / jet.pt(), 0));
    data.fill<float>("fj_genjet_mass", genjet->mass());
    data.fill<float>("fj_genOverReco_mass", catchInfs(genjet->mass() / jet.mass(), 1));
    data.fill<float>("fj_genOverReco_mass_null", catchInfs(genjet->mass() / jet.mass(), 0));
  }
  const auto *sdgenjet = jet_helper.genjetWithNuSoftDrop();
  if (sdgenjet){
    // jet here points to the uncorrected jet
    auto pos = [](double x){ return x<0 ? 0 : x; };
    data.fill<float>("fj_genjet_sdmass", pos(sdgenjet->mass()));
    data.fill<float>("fj_genjet_sdmass_sqrt", std::sqrt(pos(sdgenjet->mass())));
    data.fill<float>("fj_genOverReco_sdmass", catchInfs(pos(sdgenjet->mass()) / pos(msd_uncorr), 1));
    data.fill<float>("fj_genOverReco_sdmass_null", catchInfs(pos(sdgenjet->mass()) / pos(msd_uncorr), 0));
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

