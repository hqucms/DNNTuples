/*
 * SVFiller.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DeepNTuples/Ntupler/interface/SVFiller.h"

namespace deepntuples {



// TEMPORARY, debugging only
void printGenInfoHeader() {
  using namespace std;
  cout    << right << setw(6) << "#" << " " << setw(10) << "pdgId"
      << "  " << "Chg" << "  " << setw(10) << "Mass" << "  " << setw(48) << " Momentum"
      << left << "  " << setw(10) << "Mothers" << " " << setw(30) << "Daughters" << endl;
}

void printGenParticleInfo(const reco::GenParticle* genParticle, const int idx) {
  using namespace std;
  cout  << right << setw(3) << genParticle->status();
  cout  << right << setw(3) << idx << " " << setw(10) << genParticle->pdgId() << "  ";
  cout  << right << "  " << setw(3) << genParticle->charge() << "  " << TString::Format("%10.3g", genParticle->mass() < 1e-5 ? 0 : genParticle->mass());
  cout  << left << setw(50) << TString::Format("  (E=%6.4g pT=%6.4g eta=%7.3g phi=%7.3g)", genParticle->energy(), genParticle->pt(), genParticle->eta(), genParticle->phi());

  TString                     mothers;
  for (unsigned int iMom = 0; iMom < genParticle->numberOfMothers(); ++iMom) {
    if (mothers.Length())     mothers        += ",";
    mothers   += genParticle->motherRef(iMom).key();
  }
  cout << "  " << setw(10) << mothers;
  TString                     daughters;
  for (unsigned int iDau = 0; iDau < genParticle->numberOfDaughters(); ++iDau) {
    if (daughters.Length())   daughters      += ",";
    daughters += genParticle->daughterRef(iDau).key();
  }
  cout << " " << setw(30) << daughters << endl;
}




int hadronFlavorID(int id) {
  if ((std::abs(id) > 400 && std::abs(id) < 500) || (std::abs(id) > 4000 && std::abs(id) < 5000)) {
    return 4;
  } else if ((std::abs(id) > 500 && std::abs(id) < 600) || (std::abs(id) > 5000 && std::abs(id) < 6000)) {
    return 5;
  } else {
    return -1;
  }
}

int hadronFlavor(const reco::GenParticle gp) { // was GenParticle* gp
  int id = hadronFlavorID(gp.pdgId());
  if (id == 4) { // c -- may be b->c
    reco::GenParticle const* gp_ = &gp;
    //std::cout << "in hadrFlav, had id=4" << std::endl;
    while (gp_->numberOfMothers() != 0) {
      /*if (gp_->numberOfMothers() > 1) {
        std::cout << "    IN HADRONFLAVOR:  Found >1 mothers: " << gp_->numberOfMothers() << std::endl;
        for (unsigned i=0; i<gp_->numberOfMothers(); i++) {
          std::cout << "       mother " << i << " pdgID: " << hadronFlavorID(gp_->mother(i)->pdgId()) << ", nmothers: " << gp_->mother(i)->numberOfMothers() << ", mother ID=" << hadronFlavorID(gp_->mother(i)->mother(0)->pdgId()) << std::endl;
        }
      }*/
      if (hadronFlavorID(gp_->pdgId()) == 5) {
        //std::cout << "  found mother, hadID 5" << std::endl;
        return 10;
      }
      //std::cout << "  mother hadID incorrect, checking next mother" << std::endl;
      gp_ = (reco::GenParticle*)(gp_->mother(0));
    }
  }
  return id;
}


void SVFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
  // NEW - add additional collections for matching
  //jetToken_ = cc.consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  jetToken_ = cc.consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"));
  genParticlesToken_ = cc.consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
}

void SVFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //std::cout << "\nBEGINING readEvent" <<std::endl;
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);
  iEvent.getByToken(jetToken_, jets);
  iEvent.getByToken(genParticlesToken_, particles);

  /*std::cout << "TEMP:  PRINT GEN PARTICLE INFO" << std::endl;
  printGenInfoHeader();
  for (unsigned ipart = 0; ipart<particles->size(); ++ipart){
    printGenParticleInfo(&(particles->at(ipart)), ipart);
  }
  std::cout << "DONE, now tracking" << std::endl;
  */

  // NEW:  Perform matching here!
  // (Can't do in DeepNtuplizer bc can't pass results; can't do in fill() bc need all SVs)
  // Note:  sv results can be stored as an array; SV order should not change
  matchedIDs = new int[SVs->size()];  // NOTE:  idxsv should give you the nth entry of matchedIDs

  //check num of good jets
  int n_good_jets = 0;
  for (unsigned j=0; j<jets->size(); j++) {
    const auto& jet = jets->at(j);
    if ((jet.pt() > 40) 
       && (std::abs(jet.eta()) < 2.5) ) {
       //&& (jet.hadronFlavour() & 4) ) {
      n_good_jets++;
    }
  }
  //std::cout << "Found " << n_good_jets << " good jets" << std::endl;
  //std::cout << "Init SVs: " << SVs->size() << std::endl;

  // if sv close to a jet, assign ID of -999 and ignore it
  int jet_free_svs = 0;
  for (unsigned i=0; i<SVs->size(); i++) {
    const auto& sv = SVs->at(i);
    matchedIDs[i] = -1;  // initialize
    for (unsigned j=0; j<jets->size(); j++) {
      const auto& jet = jets->at(j);
      if (  (reco::deltaR(jet, sv) < 0.4)
         && (jet.pt() > 40)
         && (std::abs(jet.eta()) < 2.5)) {
         //&& (jet.hadronFlavour() & 4) ) {  // not 100% sure whether this is correct
         //&& (jet.jetID() & 4) ) {  // don't know how to add this...
        matchedIDs[i] = -999;  // -999 = ignore
        //std::cout << "     SV " << i << " close to jet " << j << std::endl;
      } else if (  (reco::deltaR(jet, sv) < 0.50)
         && (jet.pt() > 35)
         && (std::abs(jet.eta()) < 2.7)) {
        //std::cout << "   close match" << std::endl;
      }
    }
    if (matchedIDs[i] == -1) {
      jet_free_svs += 1;
      //std::cout << "  SV outside jet" << std::endl;
    }/* else {
      std::cout << "  SV inside jet" << std::endl;
    }*/
  }

  // - create full list of (sv, hadr pairs); dr-sorted
  //std::cout << "Jet-free SVs: " << jet_free_svs << std::endl;
  std::vector<std::tuple<float, unsigned, unsigned>> pairList; // sv num, had num
  for (unsigned i=0; i<SVs->size(); i++) {
    const auto& sv = SVs->at(i);
    if (matchedIDs[i] == -999) continue; // if sv is near a jet, skip it
    //std::cout << "found sv " << i << std::endl;
    for (unsigned j=0; j<particles->size(); j++) {
      const auto& gp = particles->at(j);
      //if (hadronFlavor(gp) > 0) std::cout << "   hadr flav: " << hadronFlavor(gp) << ", dr=" << reco::deltaR(gp, sv) << std::endl;
      if (!gp.statusFlags().isLastCopy()) continue;
      if ((hadronFlavor(gp) == 4 || hadronFlavor(gp) == 5 || hadronFlavor(gp) == 10)
          && reco::deltaR(gp, sv) < 0.4) {
        //std::cout << "   MATCHED SV " << i << " to particle " << j <<", flav=" << hadronFlavor(gp) << ", dR=" << reco::deltaR(gp, sv) << std::endl;
        pairList.push_back(std::tuple<float, int, int>(deltaR(gp, sv), i, j));
      }
    }
  }
  std::sort(pairList.begin(), pairList.end(),
    [](const std::tuple<float, unsigned, unsigned> & a, const std::tuple<float, unsigned, unsigned> & b) -> bool 
    {return std::get<0>(a) < std::get<0>(b);}  // NOTE:  Flipped direction of <!  May have been wrong before...
    );

  // - take lowest-sep pair and match; remove sv AND hadr from consideration
  //std::cout << "Assigning lowest-sep pairs" << std::endl;
  while (pairList.size() > 0) {
    // pair witih lowest dR:
    std::tuple<float, unsigned, unsigned> sel_tuple = pairList[0];
    const auto& gp = particles->at(std::get<2>(sel_tuple));
    matchedIDs[std::get<1>(sel_tuple)] = hadronFlavor(gp); // matchedIDs[sv#] = gen hadron ID
    //std::cout << "  ASSIGNED SV " << std::get<1>(sel_tuple) << " to part " << std::get<2>(sel_tuple) << ", hadr type=" << hadronFlavor(gp) << std::endl;
    // remove all occurences of this sv, hadron
    //for (unsigned i = 0; i < pairList.size(); i++) { // loop through each element of pairList, check for removal
    unsigned ind = 0;
    while (ind < pairList.size()) {
      if (std::get<1>(pairList[ind]) == std::get<1>(sel_tuple) || std::get<2>(pairList[ind]) == std::get<2>(sel_tuple)) {
        // if an element gets erased, do NOT increment ind; the next element will get slotted into ind
        pairList.erase(pairList.begin()+ind);
      } else {
        ind++;
      }
    }
  }
  //std::cout << "Finished removing" << std::endl;

  // - find lights in unmatched SVs (if dR>0.8 for all hadr)
  // note: currently finding too many lights
  int unmatched = 0;
  for (unsigned i=0; i<SVs->size(); i++) {
    if (matchedIDs[i] == -1) unmatched += 1;
  }
  //std::cout << "*Unmatched SVs: " << unmatched << std::endl;

  std::vector<std::tuple<float, unsigned, unsigned>> pairList_; // sv num, had num
  for (unsigned i=0; i<SVs->size(); i++) {
    const auto& sv = SVs->at(i);
    bool isLight = true;
    if (matchedIDs[i] == -999) continue; // if sv is near a jet, skip it
    for (unsigned j=0; j<particles->size(); j++) {
      const auto& gp = particles->at(j);
      if ((hadronFlavor(gp) == 4 || hadronFlavor(gp) == 5 || hadronFlavor(gp) == 10)
          && reco::deltaR(gp, sv) < 0.8) {
        //pairList_.push_back(std::tuple<float, unsigned, unsigned>(deltaR(gp, sv), i, j));
        isLight = false;
      }
    }
    if (isLight) {
      //if (matchedIDs[i] != -1) std::cout << "ERROR ERROR ERROR!" << std::endl;
      matchedIDs[i] = 0;
      //std::cout << "  SV " << i << " is light" << std::endl;
    }
    else if (matchedIDs[i] == -1) {
      //std::cout << "  SV " << i << " is unassigned" << std::endl;
    }
  }
  /*
  //std::cout << "Finalizing matching" << std::endl;
  std::cout << "pairList_ size is " << pairList_.size() << std::endl;
  for (unsigned i=0; i<SVs->size(); i++) {
    // if SV not found in pairList, call it a light
    bool isLight = true;
    for (unsigned j=0; j<pairList_.size(); j++) {
      std::tuple<float, unsigned, unsigned> sel_tuple = pairList_[j];
      if (std::get<1>(sel_tuple) == i) isLight = false;
    }
    if (isLight) matchedIDs[i] = 0;
  }
  */


  //finally, move -999s to -1s (now okay to ignore these)
  //for (unsigned i = 0; i<SVs->size(); i++) {
  //  if (matchedIDs[i] == -999) matchedIDs[i] = -1;
  //}

  /*std::cout << "Matching results: [" << std::endl;
  for (unsigned i=0; i<SVs->size(); i++) {
    std::cout << matchedIDs[i] << ", ";
  }
  std::cout << "]\n\n";*/
  //std::cout << "Leaving readevent" << std::endl;
}

void SVFiller::book() {

  data.add<int>("n_sv", 0);
  data.add<float>("nsv", 0);

  // basic kinematics
  //data.addMulti<float>("sv_ptrel"); // old jet vars
  //data.addMulti<float>("sv_erel");
  //data.addMulti<float>("sv_phirel");
  //data.addMulti<float>("sv_etarel");
  //data.addMulti<float>("sv_deltaR");
  data.addMulti<float>("sv_pt");
  data.addMulti<float>("sv_abseta");
  data.addMulti<float>("sv_mass");
  data.addMulti<float>("sv_energy");

  data.addMulti<float>("sv_px");
  data.addMulti<float>("sv_py");
  data.addMulti<float>("sv_pz");

  //data.addMulti<float>("sv_ptrel_log");
  //data.addMulti<float>("sv_erel_log");
  data.addMulti<float>("sv_pt_log");
  data.addMulti<float>("sv_e_log");

  // sv properties
  data.addMulti<float>("sv_ntracks");
  data.addMulti<float>("sv_chi2");
  data.addMulti<float>("sv_ndf");
  data.addMulti<float>("sv_normchi2");
  data.addMulti<float>("sv_dxy");
  data.addMulti<float>("sv_dxyerr");
  data.addMulti<float>("sv_dxysig");
  data.addMulti<float>("sv_d3d");
  data.addMulti<float>("sv_d3derr");
  data.addMulti<float>("sv_d3dsig");
  data.addMulti<float>("sv_costhetasvpv"); //pAngle in nanoAOD
  data.addMulti<float>("sv_phi");
  data.addMulti<float>("sv_neardr");

  data.addMulti<float>("sv_gen_flavor");

}

//bool SVFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {
bool SVFiller::fill(const reco::VertexCompositePtrCandidate &sv, size_t svidx, const edm::Handle<edm::View<reco::Candidate>> candHandle){
  //match svs to candHandles


  const auto &pv = vertices->at(0);

  // New:  add dR info for jets near the SV
  float dr_min = 100;
  for (unsigned j=0; j<jets->size(); j++) {
    const auto& jet = jets->at(j);
    if ((jet.pt() > 40) &&
        (std::abs(jet.eta()) < 2.5)) {
      if (reco::deltaR(jet, sv) < dr_min) {
        dr_min = reco::deltaR(jet, sv);
      }
    }
  }
  if (dr_min < 0.1) return false;
  std::cout << dr_min << std::endl;

  data.fillMulti<float>("sv_gen_flavor", matchedIDs[svidx]);

  data.fillMulti<float>("sv_pt", sv.pt());
  data.fillMulti<float>("sv_abseta", std::abs(sv.eta()));
  data.fillMulti<float>("sv_mass", sv.mass());
  data.fillMulti<float>("sv_energy", sv.energy());
  //NEW
  data.fillMulti<float>("sv_px", sv.momentum().X());
  data.fillMulti<float>("sv_py", sv.momentum().Y());
  data.fillMulti<float>("sv_pz", sv.momentum().Z());

  //data.fillMulti<float>("sv_ptrel_log", catchInfs(std::log(sv->pt()/jet.pt()), -99));
  //data.fillMulti<float>("sv_erel_log", catchInfs(std::log(sv->energy()/jet.energy()), -99));
  data.fillMulti<float>("sv_pt_log", catchInfs(std::log(sv.pt()), -99));
  data.fillMulti<float>("sv_e_log", catchInfs(std::log(sv.energy()), -99));
    
  // sv properties
  data.fillMulti<float>("sv_ntracks", sv.numberOfDaughters());
  data.fillMulti<float>("sv_chi2", sv.vertexChi2());
  data.fillMulti<float>("sv_ndf", sv.vertexNdof());
  data.fillMulti<float>("sv_normchi2", catchInfs(sv.vertexNormalizedChi2()));
    
  const auto &dxy = vertexDxy(sv, pv);
  data.fillMulti<float>("sv_dxy", dxy.value());
  data.fillMulti<float>("sv_dxyerr", dxy.error());
  data.fillMulti<float>("sv_dxysig", dxy.significance());
    
  const auto &d3d = vertexD3d(sv, pv);
  data.fillMulti<float>("sv_d3d", d3d.value());
  data.fillMulti<float>("sv_d3derr", d3d.error());
  data.fillMulti<float>("sv_d3dsig", d3d.significance());
  data.fillMulti<float>("sv_costhetasvpv", vertexDdotP(sv, pv));
  data.fillMulti<float>("sv_phi", sv.phi());

  data.fillMulti<float>("sv_neardr", dr_min);


  return true;
}


Measurement1D SVFiller::vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

Measurement1D SVFiller::vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

float SVFiller::vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv)  {
    reco::Candidate::Vector p = sv.momentum();
    reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
    return p.Unit().Dot(d.Unit());
}


} /* namespace deepntuples */



