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
  isPuppiJets_ = iConfig.getParameter<bool>("isPuppiJets");
  isQCDSample_ = iConfig.getUntrackedParameter<bool>("isQCDSample", false);
  sample_use_pythia_ = iConfig.getParameter<bool>("isPythia");
  sample_use_herwig_ = iConfig.getParameter<bool>("isHerwig");
  sample_use_madgraph_ = iConfig.getParameter<bool>("isMadGraph");
  isTrainSample_ = iConfig.getUntrackedParameter<bool>("isTrainSample", false);
}

void FatJetInfoFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(genParticlesToken_, genParticlesHandle);
}

void FatJetInfoFiller::book() {
  data.add<int>("sample_isPuppiJets", isPuppiJets_);
  data.add<int>("sample_isQCD", isQCDSample_);
  data.add<int>("sample_use_pythia", sample_use_pythia_);
  data.add<int>("sample_use_herwig", sample_use_herwig_);
  data.add<int>("sample_use_madgraph", sample_use_madgraph_); // MG can be interfaced w/ either pythia or herwig

  // truth labels
  data.add<int>("jet_partonFlavour", 0);
  data.add<int>("jet_hadronFlavour", 0);
  data.add<int>("jet_nBHadrons", 0);
  data.add<int>("jet_nCHadrons", 0);

  data.add<float>("jet_genjet_pt", 0);
  data.add<float>("jet_genjet_pt_withNu", 0);
}

bool FatJetInfoFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  // truth labels
  data.fill<int>("jet_partonFlavour", jet.partonFlavour());
  data.fill<int>("jet_hadronFlavour", jet.hadronFlavour());
  data.fill<int>("jet_nBHadrons", jet.jetFlavourInfo().getbHadrons().size());
  data.fill<int>("jet_nCHadrons", jet.jetFlavourInfo().getcHadrons().size());

  data.fill<float>("jet_genjet_pt", jet.genJet() ? jet.genJet()->pt() : -1);
  data.fill<float>("jet_genjet_pt_withNu", jet_helper.genjetWithNu() ? jet_helper.genjetWithNu()->pt() : -1);

  return true;
}

} /* namespace deepntuples */

