/*
 * DeepNtuplizer.cc
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

//#include "DeepNTuples/Ntupler/interface/JetInfoFiller.h"
//#include "DeepNTuples/Ntupler/interface/FatJetInfoFiller.h"
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

  double jetR = -1;
  double pfcandR = -1;
  bool isPuppi = true;

  edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> candToken_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection>> genJetWithNuMatchToken_;
  // NEW:
  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;

  edm::Service<TFileService> fs;
  TreeWriter *treeWriter = nullptr;

  NtupleBase* addModule(NtupleBase *m){
    modules_.push_back(m);
    return m;
  }
  std::vector<NtupleBase*> modules_;
};

DeepNtuplizer::DeepNtuplizer(const edm::ParameterSet& iConfig):
    jetR(iConfig.getParameter<double>("jetR")),
    pfcandR(iConfig.getParameter<double>("pfcandR")),
    isPuppi(iConfig.getParameter<bool>("isPuppiJets")),
    jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
    candToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pfcands"))),
    genJetWithNuMatchToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("genJetsMatch"))),
    svToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs")))
{

  // register modules
  SVFiller *sv = new SVFiller("", jetR, pfcandR);  // note:  jetR currently does nothing when passed in; should probably stay that way
  // (jetR performs different roles in sv/pf fillers)
  addModule(sv);

  PFCompleteFiller *parts = new PFCompleteFiller("", jetR, pfcandR);
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

  edm::Handle<edm::View<reco::Candidate>> candHandle;
  iEvent.getByToken(candToken_, candHandle);

  edm::Handle<edm::Association<reco::GenJetCollection>> genJetWithNuMatchHandle;
  iEvent.getByToken(genJetWithNuMatchToken_, genJetWithNuMatchHandle);

  // NEW
  edm::Handle<reco::VertexCompositePtrCandidateCollection> SVs;
  iEvent.getByToken(svToken_, SVs);

  // NEW:  Loop over SVs instead of jets
  //for (unsigned idx=0; idx<jets->size(); ++idx){
  for (unsigned idx=0; idx<SVs->size(); ++idx) {
    bool write_ = true;

    const auto& sv = SVs->at(idx);
    //const auto& jet = jets->at(idx); // need to keep the JEC for puppi sdmass corr
    //JetHelper jet_helper(&jet, candHandle, isPuppi);
    //jet_helper.setGenjetWithNu((*genJetWithNuMatchHandle)[jets->refAt(idx)]);

    for (auto *m : modules_){
      //if (!m->fillBranches(jet.correctedJet("Uncorrected"), idx, jet_helper)){
      if (!m->fillBranches(sv, idx, candHandle)) {
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
  fs->file().SetCompressionAlgorithm(ROOT::kLZ4);
  fs->file().SetCompressionLevel(4);
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
