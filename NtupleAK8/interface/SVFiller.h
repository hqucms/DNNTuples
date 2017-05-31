/*
 * SVFiller.h
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#ifndef NTUPLEAK8_INTERFACE_SVFILLER_H_
#define NTUPLEAK8_INTERFACE_SVFILLER_H_

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"

namespace deepntuples {

class SVFiller: public NtupleBase {
public:
  SVFiller() : SVFiller("", 0.8) {}
  SVFiller(std::string branchName, double jetR=0.8) : NtupleBase(branchName, jetR) {}
  virtual ~SVFiller() {}

  // get input parameters from the cfg file
  virtual void readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector && cc) override;

  // read event content or event setup for each event
  virtual void readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

protected:
  // declare the data branches (name, type, default values)
  virtual void book() override;
  // fill the branches
  virtual bool fill(const pat::Jet &jet, size_t jetidx, const JetHelper &jet_helper) override;

private:
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;

  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<reco::VertexCompositePtrCandidateCollection> SVs;

private:
  static Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv);
  static Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv);
  static float vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv);

};

} /* namespace deepntuples */

#endif /* NTUPLEAK8_INTERFACE_SVFILLER_H_ */
