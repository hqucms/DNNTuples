/*
 * JetHelper.cc
 *
 *  Created on: Jan 27, 2017
 *      Author: hqu
 */

#include "DeepNTuples/NtupleCommons/interface/JetHelper.h"

namespace deepntuples {

JetHelper::JetHelper(const pat::Jet* jet, bool has_puppi_weighted_daughters) : jet_(jet), has_puppi_weighted_daughters_(has_puppi_weighted_daughters) {
  if (!jet) throw cms::Exception("[JetHelper::JetHelper] Null pointer for input jet!");
  if (jet->nSubjetCollections() == 0) throw cms::Exception("[JetHelper::JetHelper] No subjet collection for input jet!");
  initializeConstituents();
}

void JetHelper::initializeConstituents() {
  // get subjets
  auto subjets = jet_->subjets();
  for (const auto &sj : subjets){
    subjets_.push_back(&(*sj));
  }
  // sort subjets by pt
  std::sort(subjets_.begin(), subjets_.end(),
      [](const pat::Jet* p1, const pat::Jet* p2){return p1->pt()>p2->pt();});

  // get all consitituents
  for (unsigned idau=0; idau<jet_->numberOfDaughters(); ++idau){
    const auto *dau = jet_->daughter(idau);
    if (dau->numberOfDaughters()>0){
      // is a subjet; add all daughters
      for (unsigned k=0; k<dau->numberOfDaughters(); ++k){
        const auto *cand = dynamic_cast<const pat::PackedCandidate*>(dau->daughter(k));
        if (cand->puppiWeight() < 0.01) continue; // [94X] ignore particles w/ extremely low puppi weights
        daughters_.push_back(cand);
      }
    }else{
      const auto *cand = dynamic_cast<const pat::PackedCandidate*>(dau);
      if (cand->puppiWeight() < 0.01) continue; // [94X] ignore particles w/ extremely low puppi weights
      daughters_.push_back(cand);
    }
  }
  // sort by puppi weighted pt
  if (has_puppi_weighted_daughters_){
    std::sort(daughters_.begin(), daughters_.end(),
        [](const pat::PackedCandidate* p1, const pat::PackedCandidate* p2){return p1->pt() > p2->pt();});
  } else {
    std::sort(daughters_.begin(), daughters_.end(),
        [](const pat::PackedCandidate* p1, const pat::PackedCandidate* p2){return p1->puppiWeight()*p1->pt() > p2->puppiWeight()*p2->pt();});
  }

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
