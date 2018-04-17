/*
 * JetHelper.cc
 *
 *  Created on: Jan 27, 2017
 *      Author: hqu
 */

#include "DeepNTuples/NtupleCommons/interface/JetHelper.h"

namespace deepntuples {

JetHelper::JetHelper(const pat::Jet* jet) : jet_(jet) {
  if (!jet) throw cms::Exception("[JetHelper::JetHelper] Null pointer for input jet!");
  initializeConstituents();
  computeQG();
}

void JetHelper::computeQG(bool useQualityCut) {
  // RecoJets/JetProducers/plugins/QGTagger.cc
  // Modified to use all constituents for ak8 jets

  float sum_weight = 0., sum_deta = 0., sum_dphi = 0., sum_deta2 = 0., sum_dphi2 = 0., sum_detadphi = 0., sum_pt = 0.;
  int mult = 0;

  //Loop over the jet constituents
  for(const auto *daughter : daughters_){
    auto part = static_cast<const pat::PackedCandidate*>(daughter);

    if(part->charge()){
      if(!(part->fromPV() > 1 && part->trackHighPurity())) continue;
      if(useQualityCut){
        if((part->dz()*part->dz())/(part->dzError()*part->dzError()) > 25.) continue;
        if((part->dxy()*part->dxy())/(part->dxyError()*part->dxyError()) < 25.) ++mult;
      } else ++mult;
    } else {
      if(part->pt() < 1.0) continue;
      ++mult;
    }

    float deta = daughter->eta() - jet_->eta();
    float dphi = reco::deltaPhi(daughter->phi(), jet_->phi());
    float partPt = daughter->pt();
    float weight = partPt*partPt;

    sum_weight += weight;
    sum_pt += partPt;
    sum_deta += deta*weight;
    sum_dphi += dphi*weight;
    sum_deta2 += deta*deta*weight;
    sum_detadphi += deta*dphi*weight;
    sum_dphi2 += dphi*dphi*weight;
  }

  //Calculate axis2 and ptD
  float a = 0., b = 0., c = 0.;
  float ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
  if(sum_weight > 0){
    ave_deta = sum_deta/sum_weight;
    ave_dphi = sum_dphi/sum_weight;
    ave_deta2 = sum_deta2/sum_weight;
    ave_dphi2 = sum_dphi2/sum_weight;
    a = ave_deta2 - ave_deta*ave_deta;
    b = ave_dphi2 - ave_dphi*ave_dphi;
    c = -(sum_detadphi/sum_weight - ave_deta*ave_dphi);
  }
  float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
  axis1_ = (a+b+delta > 0 ?  sqrt(0.5*(a+b+delta)) : 0);
  axis2_ = (a+b-delta > 0 ?  sqrt(0.5*(a+b-delta)) : 0);
  ptD_   = (sum_weight > 0 ? sqrt(sum_weight)/sum_pt : 0);
  multiplicity_ = mult;
}

void JetHelper::initializeConstituents() {
  try {
    // ak8: use default subjets collection
    auto subjets = jet_->subjets();
    for (const auto &sj : subjets){
	  subjets_.push_back(&(*sj));
      for (unsigned idau=0; idau<sj->numberOfDaughters(); ++idau){
        daughtersGroomed_.push_back(dynamic_cast<const pat::PackedCandidate*>(sj->daughter(idau)));
      }
    }
    std::sort(daughtersGroomed_.begin(), daughtersGroomed_.end(),
    		[](const pat::PackedCandidate* p1, const pat::PackedCandidate* p2){return p1->pt()>p2->pt();});
    std::sort(subjets_.begin(), subjets_.end(),
    		[](const pat::Jet* p1, const pat::Jet* p2){return p1->pt()>p2->pt();});

    // Then get all constituents
    for (unsigned idau=0; idau<jet_->numberOfDaughters(); ++idau){
      const auto *dau = jet_->daughter(idau);
      if (dau->numberOfDaughters()>0){
        // is a subjet; add all daughters
        for (unsigned k=0; k<dau->numberOfDaughters(); ++k){
          daughters_.push_back(dynamic_cast<const pat::PackedCandidate*>(dau->daughter(k)));
        }
      }else{
        daughters_.push_back(dynamic_cast<const pat::PackedCandidate*>(dau));
      }
    }
    std::sort(daughters_.begin(), daughters_.end(),
    		[](const pat::PackedCandidate* p1, const pat::PackedCandidate* p2){return p1->pt()>p2->pt();});

  }catch(const cms::Exception &e){
    // ak4
    for (unsigned idau=0; idau<jet_->numberOfDaughters(); ++idau){
      daughters_.push_back(dynamic_cast<const pat::PackedCandidate*>(jet_->daughter(idau)));
    }
  }
}

void JetHelper::setSubjets(const std::vector<pat::Jet>& sdjets, double R) {
  subjets_.clear();
  for (const pat::Jet &sj : sdjets) {
    // sdjets stores the soft-drop AK8 jets, with the actual subjets stored as daughters
    // PhysicsTools/PatAlgos/python/slimming/applySubstructure_cff.py
    // PhysicsTools/PatUtils/plugins/JetSubstructurePacker.cc
    if (reco::deltaR(sj, *jet_) < R) {
      for ( size_t ida = 0; ida < sj.numberOfDaughters(); ++ida ) {
        auto candPtr =  sj.daughterPtr(ida);
        subjets_.push_back(&(*edm::Ptr<pat::Jet>(candPtr)));
      }
      break;
    }
  }
}

} /* namespace deepntuples */

std::pair<double, double> deepntuples::JetHelper::getCorrectedPuppiSoftDropMass(const std::vector<const pat::Jet*> &puppisubjets) const {
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

