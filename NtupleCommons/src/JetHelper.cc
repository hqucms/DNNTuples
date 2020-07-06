/*
 * JetHelper.cc
 *
 *  Created on: Jan 27, 2017
 *      Author: hqu
 */

#include "DeepNTuples/NtupleCommons/interface/JetHelper.h"

namespace deepntuples {

JetHelper::JetHelper(const pat::Jet* jet, const edm::Handle<reco::CandidateView> &pfcands) : jet_(jet) {
  if (!jet) throw cms::Exception("[JetHelper::JetHelper] Null pointer for input jet!");
  initializeConstituents(pfcands);
}

void JetHelper::initializeConstituents(const edm::Handle<reco::CandidateView> &pfcands) {
  daughters_.clear();
  puppi_wgt_cache_.clear();

  // get all consitituents
  for (unsigned idau=0; idau<jet_->numberOfDaughters(); ++idau){
    auto dauPtr = jet_->daughterPtr(idau);
    if (dauPtr->numberOfDaughters()>0){
      // is a subjet
      const auto *sj = dynamic_cast<const pat::Jet*>(&(*dauPtr));
      if (!sj) throw cms::Exception("[JetHelper::initializeConstituents] Cannot convert to subjet!");
      // add all daughters
      for (unsigned k=0; k<sj->numberOfDaughters(); ++k){
        const auto& candPtr = sj->daughterPtr(k);
        const auto *cand = dynamic_cast<const pat::PackedCandidate*>(&(*candPtr));
        if (cand->puppiWeight() < 0.01) continue; // [94X] ignore particles w/ extremely low puppi weights
        // Here we get the original PackedCandidate as stored in MiniAOD (i.e., not puppi weighted)
        // https://github.com/cms-sw/cmssw/pull/28035
        daughters_.push_back(pfcands->ptrAt(dauPtr.key()));
        // For the Puppi weight, we get it from the new candidate in case it is recomputed
        puppi_wgt_cache_[dauPtr.key()] = cand->puppiWeight();
      }
    }else{
      const auto& candPtr = dauPtr;
      const auto *cand = dynamic_cast<const pat::PackedCandidate*>(&(*candPtr));
      if (cand->puppiWeight() < 0.01) continue; // [94X] ignore particles w/ extremely low puppi weights
      // Here we get the original PackedCandidate as stored in MiniAOD (i.e., not puppi weighted)
      // https://github.com/cms-sw/cmssw/pull/28035
      daughters_.push_back(pfcands->ptrAt(dauPtr.key()));
      // For the Puppi weight, we get it from the new candidate in case it is recomputed
      puppi_wgt_cache_[dauPtr.key()] = cand->puppiWeight();
    }
  }
  // sort by original pt
  std::sort(daughters_.begin(), daughters_.end(), [](const reco::CandidatePtr &a, const reco::CandidatePtr &b){ return a->pt() > b->pt(); });

}

} /* namespace deepntuples */
