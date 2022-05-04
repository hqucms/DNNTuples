/*
 * SVFiller.h
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#ifndef NTUPLER_INTERFACE_SVFILLER_H_
#define NTUPLER_INTERFACE_SVFILLER_H_

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"

namespace deepntuples {

class SVFiller: public NtupleBase {
public:
  SVFiller() : SVFiller("", 0.4, 0.4) {}
  SVFiller(std::string branchName, double jetR=0.4, double pfcandR=0.4) : NtupleBase(branchName, jetR, pfcandR) {}
  virtual ~SVFiller() {}

  // get input parameters from the cfg file
  virtual void readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector && cc) override;

  // read event content or event setup for each event
  virtual void readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

protected:
  // declare the data branches (name, type, default values)
  virtual void book() override;
  // fill the branches
  //virtual bool fill(const pat::Jet &jet, size_t jetidx, const JetHelper &jet_helper) override;
  virtual bool fill(const reco::VertexCompositePtrCandidate &sv, size_t svidx, const edm::Handle<edm::View<reco::Candidate>> &candHandle) override;

private:
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
  edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleRefVector> bHadronsToken_;
  edm::EDGetTokenT<reco::GenParticleRefVector> cHadronsToken_;

  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<reco::VertexCompositePtrCandidateCollection> SVs;
  // NEW, matching
  edm::Handle<edm::View<pat::Jet>> jets;
  edm::Handle<reco::GenParticleCollection> particles;
  edm::Handle<reco::GenParticleRefVector> bhadrons;
  edm::Handle<reco::GenParticleRefVector> chadrons;

  int *matchedIDs; //array of all gen labels in the event
  float *lightdr;
  float *hadrdr; // array of all dR vals of matched hadrons
  int *light_dq; // if particle gets DQed from being a light, specify ID of particle that DQed it.

private:
  static Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv);
  static Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv);
  static float vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv);

};

} /* namespace deepntuples */

#endif /* NTUPLER_INTERFACE_SVFILLER_H_ */
