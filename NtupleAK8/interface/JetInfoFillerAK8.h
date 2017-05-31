/*
 * JetInfoAK8.h
 *
 *  Created on: May 23, 2017
 *      Author: hqu
 */

#ifndef NTUPLEAK8_INTERFACE_JETINFOFILLERAK8_H_
#define NTUPLEAK8_INTERFACE_JETINFOFILLERAK8_H_

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"
#include "DeepNTuples/BTagHelpers/interface/FlavorDefinition.h"

namespace deepntuples {

class JetInfoFillerAK8: public NtupleBase {
public:
  JetInfoFillerAK8() : JetInfoFillerAK8("") {}
  JetInfoFillerAK8(std::string branchName, double jetR=0.8) : NtupleBase(branchName, jetR), flavorDef(jetR) {}
  virtual ~JetInfoFillerAK8() {}

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

//  edm::EDGetTokenT<edm::Association<reco::GenJetCollection>> genJetMatchReclusterToken_;
//  edm::EDGetTokenT<edm::Association<reco::GenJetCollection>> genJetMatchWithNuToken_;
//  edm::Handle<edm::Association<reco::GenJetCollection>> genJetMatchRecluster;
//  edm::Handle<edm::Association<reco::GenJetCollection>> genJetMatchWithNu;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;


//  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
//  edm::Handle<pat::MuonCollection> muonsHandle;
//
//  edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;
//  edm::Handle<pat::ElectronCollection> electronsHandle;

};

} /* namespace deepntuples */

#endif /* NTUPLEAK8_INTERFACE_JETINFOFILLERAK8_H_ */
