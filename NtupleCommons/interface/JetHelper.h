/*
 * JetHelper.hh
 *
 *  Created on: Jan 27, 2017
 *      Author: hqu
 */

#ifndef NTUPLECOMMONS_INTERFACE_JETHELPER_H_
#define NTUPLECOMMONS_INTERFACE_JETHELPER_H_

#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

namespace deepntuples {

class JetHelper {
public:
  JetHelper() {}
  JetHelper(const pat::Jet *jet, bool has_puppi_weighted_daughters);

  virtual ~JetHelper() {}

  // set genjet (clustered w/ neutrino)
  void setGenjetWithNu(const reco::GenJetRef &genjetRef) { genjetWithNu_ = (genjetRef.isNull() ? nullptr : &(*genjetRef)); }
  void setGenjetWithNuSoftDrop(const reco::GenJetRef &genjetRef) { genjetWithNuSoftDrop_ = (genjetRef.isNull() ? nullptr : &(*genjetRef)); }
  // ------

  // return jet constituents (PF candidates)
  const std::vector<const pat::PackedCandidate*>& getJetConstituents() const { return daughters_; }
  unsigned int numberOfDaughters() const { return daughters_.size(); }
  bool hasPuppiWeightedDaughters() const { return has_puppi_weighted_daughters_; }

  const pat::Jet& jet() const { return *jet_; }
  const std::vector<const pat::Jet*>& getSubJets() const { return subjets_; }

  const reco::GenJet* genjetWithNu() const { return genjetWithNu_; }
  const reco::GenJet* genjetWithNuSoftDrop() const { return genjetWithNuSoftDrop_; }

  std::pair<double, double> getCorrectedPuppiSoftDropMass(const std::vector<const pat::Jet*> &puppisubjets) const; // tmp


private:
  void initializeConstituents();


private:
  // data members
  const pat::Jet *jet_ = nullptr;
  bool has_puppi_weighted_daughters_ = false;
  const reco::GenJet *genjetWithNu_ = nullptr;
  const reco::GenJet *genjetWithNuSoftDrop_ = nullptr;
  std::vector<const pat::Jet*> subjets_;
  std::vector<const pat::PackedCandidate*> daughters_;

};

} /* namespace deepntuples */

#endif /* NTUPLECOMMONS_INTERFACE_JETHELPER_H_ */
