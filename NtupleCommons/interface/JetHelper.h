/*
 * JetHelper.hh
 *
 *  Created on: Jan 27, 2017
 *      Author: hqu
 */

#ifndef NTUPLECOMMONS_INTERFACE_JETHELPER_H_
#define NTUPLECOMMONS_INTERFACE_JETHELPER_H_

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

namespace deepntuples {

class JetHelper {
public:
  JetHelper() {}
  JetHelper(const pat::Jet *jet);

  virtual ~JetHelper() {}

  // set subjets
  void setSubjets(const std::vector<pat::Jet>& sdjets, double R);

  // set genjet (clustered w/ neutrino)
  void setGenjetWithNu(const reco::GenJetRef &genjetRef) { genjetWithNu_ = (genjetRef.isNull() ? nullptr : &(*genjetRef)); }
  void setGenjetWithNuSoftDrop(const reco::GenJetRef &genjetRef) { genjetWithNuSoftDrop_ = (genjetRef.isNull() ? nullptr : &(*genjetRef)); }
  // ------

  // return jet constituents (PF candidates)
  const std::vector<const pat::PackedCandidate*>& getJetConstituents() const { return daughters_; }
  unsigned int numberOfDaughters() const { return daughters_.size(); }

  const std::vector<const pat::PackedCandidate*>& getGroomedJetConstituents() const { return daughtersGroomed_; }

  const std::vector<const pat::Jet*>& getSubJets() const { return subjets_; }

  const reco::GenJet* genjetWithNu() const { return genjetWithNu_; }
  const reco::GenJet* genjetWithNuSoftDrop() const { return genjetWithNuSoftDrop_; }

  std::pair<double, double> getCorrectedPuppiSoftDropMass(const std::vector<const pat::Jet*> &puppisubjets) const; // tmp

  // quark/gluon discrimination variables
  double ptD()   const { return ptD_;          }
  double axis1() const { return axis1_;        }
  double axis2() const { return axis2_;        }
  int    mult()  const { return multiplicity_; }


private:
  void initializeConstituents();
  void computeQG(bool useQualityCut = false);


private:
  // data members
  const pat::Jet *jet_ = nullptr;
  const reco::GenJet *genjetWithNu_ = nullptr;
  const reco::GenJet *genjetWithNuSoftDrop_ = nullptr;
  std::vector<const pat::Jet*> subjets_;
  std::vector<const pat::PackedCandidate*> daughters_;
  std::vector<const pat::PackedCandidate*> daughtersGroomed_;

  double ptD_ = -1;
  double axis1_ = -1;
  double axis2_ = -1;
  int multiplicity_ = -1;


};

} /* namespace deepntuples */

#endif /* NTUPLECOMMONS_INTERFACE_JETHELPER_H_ */
