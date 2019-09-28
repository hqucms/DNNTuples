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
  if (jet->nSubjetCollections() == 0) throw cms::Exception("[JetHelper::JetHelper] No subjet collection for input jet!");
  initializeConstituents(pfcands);
}

void JetHelper::initializeConstituents(const edm::Handle<reco::CandidateView> &pfcands) {
  subjets_.clear();
  uncorr_subjets_.clear();
  daughters_.clear();
  puppi_wgt_cache_.clear();

  // get subjets
  auto subjets = jet_->subjets();
  for (const auto &sj : subjets){
    subjets_.push_back(&(*sj));
    uncorr_subjets_.push_back(&(*sj));
  }
  // sort subjets by pt
  std::sort(subjets_.begin(), subjets_.end(),
      [](const pat::Jet* p1, const pat::Jet* p2){return p1->pt()>p2->pt();});

  std::sort(uncorr_subjets_.begin(), uncorr_subjets_.end(),
      [](const pat::Jet* p1, const pat::Jet* p2){return p1->correctedP4("Uncorrected").pt()>p2->correctedP4("Uncorrected").pt();});

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

std::pair<double, double> JetHelper::getCorrectedPuppiSoftDropMass(const std::vector<const pat::Jet*> &puppisubjets) const {
  double sdpuppimass = 0;
  if (puppisubjets.size()==1){
    sdpuppimass = puppisubjets[0]->correctedP4(0).mass();
  }else if (puppisubjets.size()>=2){
    sdpuppimass = (puppisubjets[0]->correctedP4(0) + puppisubjets[1]->correctedP4(0)).mass();
  }
  double pt = jet_->pt();
  double eta = jet_->eta();
  double gencorr = 1.006261 + ((-1.061605) * pow(pt*0.079990,-1.204538));
  double recocorr = 1;
  if (std::abs(eta) <= 1.3){
    recocorr = 1.093020+(-0.000150068)*pt+(3.44866e-07)*pow(pt,2)+(-2.68100e-10)*pow(pt,3)+(8.67440e-14)*pow(pt,4)+(-1.00114e-17)*pow(pt,5);
  }else{
    recocorr = 1.272115+(-0.000571640)*pt+(8.37289e-07)*pow(pt,2)+(-5.20433e-10)*pow(pt,3)+(1.45375e-13)*pow(pt,4)+(-1.50389e-17)*pow(pt,5);
  }
  return std::make_pair(sdpuppimass, sdpuppimass*gencorr*recocorr);
}

} /* namespace deepntuples */
