/*
 * FlavorDefinition.h
 *
 *      Author: mverzett
 *      Modified: hqu
 */

#ifndef BTAGHELPERS_INTERFACE_FLAVORDEFINITION_H_
#define BTAGHELPERS_INTERFACE_FLAVORDEFINITION_H_

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

namespace deepntuples {

enum JetFlavor {UNDEFINED, G, UD, S, C, B, BB, LeptonicB, LeptonicB_C};

class FlavorDefinition {
public:
  FlavorDefinition() {}
  FlavorDefinition(double jetR) : jetR_(jetR) {}
  virtual ~FlavorDefinition() {}

  void setGenParticles(const reco::GenParticleCollection &genParticles);

  JetFlavor jet_flavour(const pat::Jet& jet, bool usePhysForLightAndUndefined=false) const;

  std::vector<std::size_t> jet_muonsIds(const pat::Jet& jet, const std::vector<pat::Muon>& event_muons) const;
  std::vector<std::size_t> jet_electronsIds(const pat::Jet& jet, const std::vector<pat::Electron>& event_electrons) const;


private:
  double jetR_ = 0.4;
  std::vector<const reco::GenParticle*> neutrinosLepB;
  std::vector<const reco::GenParticle*> neutrinosLepB_C;

};

} /* namespace deepntuples */

#endif /* BTAGHELPERS_INTERFACE_FLAVORDEFINITION_H_ */
