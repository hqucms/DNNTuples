/*
 * JetInfoAK8.h
 *
 *  Created on: May 23, 2017
 *      Author: hqu
 */

#ifndef NTUPLER_INTERFACE_JETINFOFILLER_H_
#define NTUPLER_INTERFACE_JETINFOFILLER_H_

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"
#include "DeepNTuples/BTagHelpers/interface/FlavorDefinition.h"

namespace deepntuples {

class JetInfoFiller: public NtupleBase {
public:
  JetInfoFiller() : JetInfoFiller("") {}
  JetInfoFiller(std::string branchName, double jetR=0.8) : NtupleBase(branchName, jetR), flavorDef(jetR),
      jetIdTight(PFJetIDSelectionFunctor::SUMMER18PUPPI, PFJetIDSelectionFunctor::TIGHT),
      jetIdTightLepVeto(PFJetIDSelectionFunctor::SUMMER18PUPPI, PFJetIDSelectionFunctor::TIGHTLEPVETO){}
  virtual ~JetInfoFiller() {}

  // get input parameters from the cfg file
  virtual void readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector && cc) override;

  // read event content or event setup for each event
  virtual void readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

protected:
  // declare the data branches (name, type, default values)
  virtual void book() override;
  // fill the branches
  virtual bool fill(const pat::Jet &jet, size_t jetidx, const JetHelper &jet_helper) override;

private:
  FlavorDefinition flavorDef;

  double minPt_ = 0;
  double maxPt_ = 0;
  double maxAbsEta_ = 0;
  std::vector<std::string> btag_discriminators_;

  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
  edm::EDGetTokenT<double> rhoToken_;

  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
  edm::Handle<double> rhoInfo;

  unsigned event_ = 0;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;

  PFJetIDSelectionFunctor jetIdTight;
  PFJetIDSelectionFunctor jetIdTightLepVeto;

};

} /* namespace deepntuples */

#endif /* NTUPLER_INTERFACE_JetInfoFiller_H_ */
