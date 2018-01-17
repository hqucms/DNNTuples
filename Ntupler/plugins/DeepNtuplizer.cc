/*
 * DeepNtuplizerAK8.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DeepNTuples/NtupleCommons/interface/TreeWriter.h"

#include "DeepNTuples/Ntupler/interface/JetInfoFiller.h"
#include "DeepNTuples/Ntupler/interface/FatJetInfoFiller.h"
#include "DeepNTuples/Ntupler/interface/PFCandidateFiller.h"
#include "DeepNTuples/Ntupler/interface/TrackFiller.h"
#include "DeepNTuples/Ntupler/interface/SVFiller.h"
#include "DeepNTuples/Ntupler/interface/PFCompleteFiller.h"


using namespace deepntuples;

class DeepNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DeepNtuplizer(const edm::ParameterSet&);
  ~DeepNtuplizer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;

  edm::Service<TFileService> fs;
  TreeWriter *treeWriter = nullptr;

  NtupleBase* addModule(NtupleBase *m){
    modules_.push_back(m);
    return m;
  }
  std::vector<NtupleBase*> modules_;
};

DeepNtuplizer::DeepNtuplizer(const edm::ParameterSet& iConfig):
    jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets")))
{

  // read configuration parameters
  const double jetR = iConfig.getParameter<double>("jetR");

  // register modules
  JetInfoFiller *jetinfo = new JetInfoFiller("", jetR);
  addModule(jetinfo);

  FatJetInfoFiller *fjinfo = new FatJetInfoFiller("", jetR);
  addModule(fjinfo);

  PFCandidateFiller *pfcands = new PFCandidateFiller("", jetR);
  addModule(pfcands);

  TrackFiller *tracks = new TrackFiller("", jetR);
  addModule(tracks);

  SVFiller *sv = new SVFiller("", jetR);
  addModule(sv);

  PFCompleteFiller *parts = new PFCompleteFiller("", jetR);
  addModule(parts);

  // read config and init modules
  for(auto& m: modules_)
    m->readConfig(iConfig, consumesCollector());

}

DeepNtuplizer::~DeepNtuplizer()
{
  for(auto *m : modules_)
    delete m;
}


// ------------ method called for each event  ------------
void DeepNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  for(auto *m : modules_){
    m->readEvent(iEvent, iSetup);
  }

  edm::Handle<edm::View<pat::Jet>> jets;
  iEvent.getByToken(jetToken_, jets);

  for (unsigned idx=0; idx<jets->size(); ++idx){
    bool write_ = true;

    const auto& jet = jets->at(idx); // need to keep the JEC for puppi sdmass corr
    JetHelper jet_helper(&jet);

    for (auto *m : modules_){
      if (!m->fillBranches(jet.correctedJet("Uncorrected"), idx, jet_helper)){
        write_ = false;
        break;
      }
    }

    if (write_) {
      treeWriter->fill();
    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void DeepNtuplizer::beginJob() {
  if( !fs ){
    throw edm::Exception( edm::errors::Configuration,
        "TFile Service is not registered in cfg file" );
  }
  treeWriter = new TreeWriter(fs->make<TTree>("tree" ,"tree"));

  for(auto *m : modules_)
    m->initBranches(treeWriter);

}

// ------------ method called once each job just after ending the event loop  ------------
void DeepNtuplizer::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DeepNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeepNtuplizer);
