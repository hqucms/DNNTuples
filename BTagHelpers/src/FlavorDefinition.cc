#include "DeepNTuples/BTagHelpers/interface/FlavorDefinition.h"

namespace deepntuples {

void FlavorDefinition::setGenParticles(const reco::GenParticleCollection& genParticles) {
  neutrinosLepB.clear();
  neutrinosLepB_C.clear();

  for (const auto &gen : genParticles) {
    if(abs(gen.pdgId())==12||abs(gen.pdgId())==14||abs(gen.pdgId())==16) {
      const reco::GenParticle* mother = dynamic_cast<const reco::GenParticle*>(gen.mother());
      if(mother) {
        if( (abs(mother->pdgId())>500 && abs(mother->pdgId())<600)
            || (abs(mother->pdgId())>5000 && abs(mother->pdgId())<6000) ) {
          neutrinosLepB.push_back(&gen);
        }
        if( (abs(mother->pdgId())>400 && abs(mother->pdgId())<500)
            || (abs(mother->pdgId())>4000 && abs(mother->pdgId())<5000) ) {
          neutrinosLepB_C.push_back(&gen);
        }
      }
      else {
        std::cout << "No mother" << std::endl;
      }
    }
  }

}

JetFlavor FlavorDefinition::jet_flavour(const pat::Jet& jet, bool usePhysForLightAndUndefined) const {
  int hflav = abs(jet.hadronFlavour());
  int pflav = abs(jet.partonFlavour());
  int physflav = 0;
  if(jet.genParton()) physflav=abs(jet.genParton()->pdgId());
  std::size_t nbs = jet.jetFlavourInfo().getbHadrons().size();
  std::size_t ncs = jet.jetFlavourInfo().getcHadrons().size();

  if(hflav == 5) { //B jet
    if(nbs > 1) return JetFlavor::BB;
    else if(nbs == 1) {
      for (const auto *gen : neutrinosLepB){
        if(reco::deltaR(*gen, jet) < jetR_) {
          return JetFlavor::LeptonicB;
        }
      }
      for (const auto *gen : neutrinosLepB_C){
        if(reco::deltaR(*gen, jet) < jetR_) {
          return JetFlavor::LeptonicB_C;
        }
      }
      return JetFlavor::B;
    }
    else {
      if(usePhysForLightAndUndefined){
        if(physflav == 21) return JetFlavor::G;
        else if(physflav == 3) return JetFlavor::S;
        else if(physflav == 2 || physflav ==1) return JetFlavor::UD;
        else return JetFlavor::UNDEFINED;
      }
      else return JetFlavor::UNDEFINED;
    }
  }
  else if(hflav == 4) { //C jet
    if(ncs > 1) return JetFlavor::CC;
    else return JetFlavor::C;
  }
  else { //not a heavy jet
    if(std::abs(pflav) == 4 || std::abs(pflav) == 5 || nbs || ncs) {
      if(usePhysForLightAndUndefined){
        if(physflav == 21) return JetFlavor::G;
        else if(physflav == 3) return JetFlavor::S;
        else if(physflav == 2 || physflav ==1) return JetFlavor::UD;
        else return JetFlavor::UNDEFINED;
      }
      else return JetFlavor::UNDEFINED;
    }
    else if(usePhysForLightAndUndefined){
      if(physflav == 21) return JetFlavor::G;
      else if(physflav == 3) return JetFlavor::S;
      else if(physflav == 2 || physflav ==1) return JetFlavor::UD;
      else return JetFlavor::UNDEFINED;
    }
    else {
      if(pflav == 21) return JetFlavor::G;
      else if(pflav == 3) return JetFlavor::S;
      else if(pflav == 2 || pflav ==1) return JetFlavor::UD;
      else return JetFlavor::UNDEFINED;
    }
  }

}

std::vector<std::size_t> FlavorDefinition::jet_muonsIds(const pat::Jet& jet, const std::vector<pat::Muon>& event_muons) const {
  std::vector <std::size_t> muonsIds;
  for (std::size_t i = 0; i < event_muons.size(); i++) {
    const auto & muon = event_muons.at(i);
    if(reco::deltaR(muon, jet) < jetR_) muonsIds.emplace_back(i);
  }
  return muonsIds;

}

std::vector<std::size_t> FlavorDefinition::jet_electronsIds(const pat::Jet& jet, const std::vector<pat::Electron>& event_electrons) const {
  std::vector <std::size_t> electronsIds;
  for (std::size_t i = 0; i < event_electrons.size(); i++) {
    const auto & electron = event_electrons.at(i);
    if(reco::deltaR(electron, jet) < jetR_) electronsIds.emplace_back(i);
  }
  return electronsIds;
}

} /* namespace deepntuples */

