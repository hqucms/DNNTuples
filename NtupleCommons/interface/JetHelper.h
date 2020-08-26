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

#include <map>

namespace deepntuples {

class JetHelper {
public:
  JetHelper() {}
  JetHelper(const pat::Jet *jet, const edm::Handle<reco::CandidateView> &pfcands, bool isPuppi=true);

  virtual ~JetHelper() {}

  // set genjet (clustered w/ neutrino)
  void setGenjetWithNu(const reco::GenJetRef &genjetRef) { genjetWithNu_ = (genjetRef.isNull() ? nullptr : &(*genjetRef)); }
  // ------

  // return jet constituents (PF candidates)
  const std::vector<reco::CandidatePtr>& getJetConstituents() const { return daughters_; }
  unsigned int numberOfDaughters() const { return daughters_.size(); }
  float getPuppiWeight(const reco::CandidatePtr &cand) const {
    auto iter = puppi_wgt_cache_.find(cand.key());
    if (iter == puppi_wgt_cache_.end()){
      throw cms::Exception("[JetHelper::getPuppiWeight] Cannot get puppi wgt!");
    }
    return iter->second;
  }

  const pat::Jet& jet() const { return *jet_; }
  const reco::GenJet* genjetWithNu() const { return genjetWithNu_; }


private:
  void initializeConstituents(const edm::Handle<reco::CandidateView> &pfcands);


private:
  // data members
  const pat::Jet *jet_ = nullptr;
  const bool isPuppi_ = true;
  const reco::GenJet *genjetWithNu_ = nullptr;
  std::vector<reco::CandidatePtr> daughters_;
  std::map<reco::CandidatePtr::key_type, float> puppi_wgt_cache_;

};

} /* namespace deepntuples */

#endif /* NTUPLECOMMONS_INTERFACE_JETHELPER_H_ */
