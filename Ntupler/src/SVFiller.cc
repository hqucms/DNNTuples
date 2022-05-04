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

int hadronFlavor(const reco::GenParticle& gp) { // was GenParticle* gp
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
  bHadronsToken_ = cc.consumes<reco::GenParticleRefVector>(iConfig.getParameter<edm::InputTag>("bHadrons"));
  cHadronsToken_ = cc.consumes<reco::GenParticleRefVector>(iConfig.getParameter<edm::InputTag>("cHadrons"));
}

void SVFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //std::cout << "\nBEGINING readEvent" <<std::endl;
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);
  iEvent.getByToken(jetToken_, jets);
  iEvent.getByToken(genParticlesToken_, particles);
  iEvent.getByToken(bHadronsToken_, bhadrons);
  iEvent.getByToken(cHadronsToken_, chadrons);

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
  lightdr    = new float[SVs->size()];
  hadrdr     = new float[SVs->size()];
  light_dq   = new int[SVs->size()];

  // if sv close to a jet, assign ID of -2 and ignore it
  
  for (unsigned i=0; i<SVs->size(); i++) {
    const auto& sv = SVs->at(i);
    matchedIDs[i] = -1;  // initialize
    lightdr[i] = 10;
    hadrdr[i] = 10;
    light_dq[i] = -1;
    for (unsigned j=0; j<jets->size(); j++) {
      const auto& jet = jets->at(j);
      if (  (reco::deltaR(jet, sv) < 0.4)
         && (jet.pt() > 40)
         && (std::abs(jet.eta()) < 2.5)) {
         //&& (jet.hadronFlavour() & 4) ) {  // not 100% sure whether this is correct
         //&& (jet.jetID() & 4) ) {  // don't know how to add this...
        matchedIDs[i] = -2;  // -2 = ignore
        //std::cout << "     SV " << i << " close to jet " << j << std::endl;
      }
    }
    //if (matchedIDs[i] == -2) {
    //  std::cout << "  SV " << i << " is unmatched" << std::endl;
    //}
  }

  /*std::cout << "CURRENT SV LIST:" << std::endl;
  for (unsigned i=0; i<SVs->size(); i++) {
    std::cout << "(" << i << ", " << matchedIDs[i] << "),  ";
  }
  std::cout << std::endl;
  */
  // - create full list of (sv, hadr pairs); dr-sorted
  //std::cout << "Jet-free SVs: " << jet_free_svs << std::endl;
  std::vector<std::tuple<float, unsigned, unsigned>> pairList; // sv num, had num
  for (unsigned i=0; i<SVs->size(); i++) {
    const auto& sv = SVs->at(i);
    if (matchedIDs[i] == -2) continue; // if sv is near a jet, skip it
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
    hadrdr[std::get<1>(sel_tuple)] = std::get<0>(sel_tuple);
    //if (hadronFlavor(gp) > 10) {
    //  std::cout << "  ASSIGNED SV " << std::get<1>(sel_tuple) << " to part " << std::get<2>(sel_tuple) << ", hadr type=" << hadronFlavor(gp) << std::endl;
    //}
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

  for (unsigned i=0; i<SVs->size(); i++) {
    const auto& sv = SVs->at(i);
    bool isLight = true;
    if (matchedIDs[i] == -2) continue; // if sv is near a jet, skip it
    float dr_dq = 1;  // dR of closest pfcand that DQed the SV as a light
    for (unsigned j=0; j<particles->size(); j++) {
      const auto& gp = particles->at(j);
      if ((hadronFlavor(gp) == 4 || hadronFlavor(gp) == 5 || hadronFlavor(gp) == 10)
          && reco::deltaR(gp, sv) < 0.8) {
        isLight = false;
        if (reco::deltaR(gp, sv) < dr_dq && reco::deltaR(gp, sv) > 0.4) {
          light_dq[i] = hadronFlavor(gp); // store flav of part that DQed the SV
          dr_dq = reco::deltaR(gp, sv);
        }
        if (reco::deltaR(gp, sv) < lightdr[i]) lightdr[i] = reco::deltaR(gp, sv);
      }
    }
    if (isLight) {
      //if (matchedIDs[i] != -1) std::cout << "ERROR ERROR ERROR!" << std::endl;
      matchedIDs[i] = 0;
      //std::cout << "  SV " << i << " is light" << std::endl;
    }
  }

  /*std::cout << "FINAL SV LIST:" << std::endl;
  for (unsigned i=0; i<SVs->size(); i++) {
    std::cout << matchedIDs[i] << ",  ";
  }
  std::cout << std::endl;
  */

  /*std::cout << "Matching results: [" << std::endl;
  for (unsigned i=0; i<SVs->size(); i++) {
    std::cout << matchedIDs[i] << ", ";
  }
  std::cout << "]\n\n";*/
  //std::cout << "Leaving readevent" << std::endl;
}

void SVFiller::book() {
  // config info, for bookkeeping
  data.add<float>("jetR", jetR_);
  data.add<float>("pfcandR", pfcandR_);

  data.add<int>("n_pv", 0);

  // basic kinematics
  data.add<float>("sv_pt", -100);
  data.add<float>("sv_eta", 4);
  data.add<float>("sv_phi", -5);
  data.add<float>("sv_mass", -5);
  data.add<float>("sv_abseta", 4);

  data.add<float>("sv_px", 999);
  data.add<float>("sv_py", 999);
  data.add<float>("sv_pz", 999);
  data.add<float>("sv_energy", -100);

  //data.add<float>("sv_ptrel_log");
  //data.add<float>("sv_erel_log");
  data.add<float>("sv_pt_log", -1);
  data.add<float>("sv_e_log", -1);

  // sv properties
  data.add<float>("sv_ntracks", -1);
  data.add<float>("sv_chi2", -5);
  data.add<float>("sv_ndf", -5);
  data.add<float>("sv_normchi2", -35000);
  data.add<float>("sv_dxy", -1);
  data.add<float>("sv_dxyerr", -1);
  data.add<float>("sv_dxysig", -100);
  data.add<float>("sv_d3d", -5);
  data.add<float>("sv_d3derr", -1);
  data.add<float>("sv_d3dsig", -100);
  data.add<float>("sv_costhetasvpv", -2); //pAngle in nanoAOD

  data.add<float>("sv_neardr", -1); // nearest jet with dR>0.1
  data.add<float>("sv_lightdr", -1); // lowest dR of nearby 
  data.add<float>("sv_hadrdr", -1); // if matched, dR w/ matched part; else 10.
  data.add<int>("sv_light_dq", -1); // flav of the particle that DQed the SV from being a light

  data.add<float>("sv_gen_flavor", -1);

  // NEW, for reweighting:
  data.add<int>("sv_sample_label", -1);
  // Number of c/bs in cone around SV (ignoring cs produced from a b)
  data.add<int>("sv_n_c", -1);
  data.add<int>("sv_n_b", -1);

}

bool SVFiller::fill(const reco::VertexCompositePtrCandidate &sv, size_t svidx, const edm::Handle<edm::View<reco::Candidate>> &candHandle){
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
  if (dr_min < jetR_) {
    return false;
  }

  int n_b = 0; // # bs in cone
  int n_c = 0; // # cs in cone (can be from b)

  for (unsigned j = 0; j < bhadrons->size(); j++) {
    if (reco::deltaR(*bhadrons->at(j), sv) < pfcandR_) n_b++;
  }
  for (unsigned j = 0; j < chadrons->size(); j++) {
    if (reco::deltaR(*chadrons->at(j), sv) < pfcandR_) n_c++;
  }

  //if (matchedIDs[svidx]!=10 && matchedIDs[svidx]!=5 && matchedIDs[svidx]!=4 && matchedIDs[svidx]!=0) {
  //  //std::cout << "ASSIGNED " << matchedIDs[svidx] << std::endl;
  //  matchedIDs[svidx] = -1;
  //}
  data.fill<float>("sv_gen_flavor", matchedIDs[svidx]);

  data.fill<int>("n_pv", vertices->size());

  data.fill<float>("sv_pt", sv.pt());
  data.fill<float>("sv_eta", sv.eta());
  data.fill<float>("sv_phi", sv.phi());
  data.fill<float>("sv_mass", sv.mass());
  data.fill<float>("sv_abseta", std::abs(sv.eta()));

  data.fill<float>("sv_px", sv.px());
  data.fill<float>("sv_py", sv.py());
  data.fill<float>("sv_pz", sv.pz());
  data.fill<float>("sv_energy", sv.energy());

  data.fill<float>("sv_pt_log", catchInfs(std::log(sv.pt()), -99));
  data.fill<float>("sv_e_log", catchInfs(std::log(sv.energy()), -99));
    
  // sv properties
  data.fill<float>("sv_ntracks", sv.numberOfDaughters());
  data.fill<float>("sv_chi2", sv.vertexChi2());
  data.fill<float>("sv_ndf", sv.vertexNdof());
  data.fill<float>("sv_normchi2", catchInfs(sv.vertexNormalizedChi2()));
    
  const auto &dxy = vertexDxy(sv, pv);
  data.fill<float>("sv_dxy", dxy.value());
  data.fill<float>("sv_dxyerr", dxy.error());
  data.fill<float>("sv_dxysig", dxy.significance());
    
  const auto &d3d = vertexD3d(sv, pv);
  data.fill<float>("sv_d3d", d3d.value());
  data.fill<float>("sv_d3derr", d3d.error());
  data.fill<float>("sv_d3dsig", d3d.significance());
  data.fill<float>("sv_costhetasvpv", vertexDdotP(sv, pv));

  data.fill<float>("sv_neardr", dr_min);
  data.fill<float>("sv_lightdr", lightdr[svidx]);
  data.fill<float>("sv_hadrdr", hadrdr[svidx]);
  data.fill<int>("sv_light_dq", light_dq[svidx]);

  data.fill<int>("sv_sample_label", 2);

  data.fill<int>("sv_n_c", n_c);
  data.fill<int>("sv_n_b", n_b);


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



