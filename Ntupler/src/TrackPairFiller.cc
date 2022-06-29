#include <unordered_map>
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DeepNTuples/Ntupler/interface/TrackPairFiller.h"

namespace deepntuples {

  void TrackPairFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
    transientTrackBuilderToken_ =
        cc.esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));
    vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
    svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
  }

  void TrackPairFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    iEvent.getByToken(vtxToken_, vertices);
    iEvent.getByToken(svToken_, SVs);
    builder_ = iSetup.getHandle(transientTrackBuilderToken_);
  }

  void TrackPairFiller::book() {
    data.addMulti<int>("idx_trk1");
    data.addMulti<int>("idx_trk2");

    // pca (1 to 2)
    data.addMulti<float>("trkpair_pca_distance");
    data.addMulti<float>("trkpair_pca_significance");

    // //   pcas (1 and 2) poistions
    // data.addMulti<float>("trkpair_pcaSeed_x1");
    // data.addMulti<float>("trkpair_pcaSeed_y1");
    // data.addMulti<float>("trkpair_pcaSeed_z1");

    // data.addMulti<float>("trkpair_pcaSeed_x2");
    // data.addMulti<float>("trkpair_pcaSeed_y2");
    // data.addMulti<float>("trkpair_pcaSeed_z2");

    // data.addMulti<float>("trkpair_pcaSeed_xerr1");
    // data.addMulti<float>("trkpair_pcaSeed_yerr1");
    // data.addMulti<float>("trkpair_pcaSeed_zerr1");

    // data.addMulti<float>("trkpair_pcaSeed_xerr2");
    // data.addMulti<float>("trkpair_pcaSeed_yerr2");
    // data.addMulti<float>("trkpair_pcaSeed_zerr2");

    // dot prod betweeen track and pca direction
    data.addMulti<float>("trkpair_dotprod1");
    data.addMulti<float>("trkpair_dotprod2");

    //pca distance form PV on both tracks
    data.addMulti<float>("trkpair_pca_dist1");
    data.addMulti<float>("trkpair_pca_dist2");

    //track track or dir dir dotprod
    data.addMulti<float>("trkpair_dotprod12_2D");
    data.addMulti<float>("trkpair_dotprod12_2DV");
    data.addMulti<float>("trkpair_dotprod12_3D");
    data.addMulti<float>("trkpair_dotprod12_3DV");

    //jet pca relative
    data.addMulti<float>("trkpair_pca_jetAxis_dist");
    data.addMulti<float>("trkpair_pca_jetAxis_dotprod");
    // data.addMulti<float>("trkpair_pca_jetAxis_dEta");
    // data.addMulti<float>("trkpair_pca_jetAxis_dPhi");
  }

  bool TrackPairFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {
    const auto& pv = vertices->at(0);
    const auto& pfCands = jet_helper.getJetConstituents();
    GlobalVector jetdirection(jet.px(), jet.py(), jet.pz());

    for (unsigned i = 0; i < pfCands.size(); ++i) {
      const auto& cand_i = pfCands[i];
      if (!cand_i->bestTrack())
        continue;
      const auto trk_i = builder_->build(cand_i);
      for (unsigned j = i + 1; j < pfCands.size(); ++j) {
        if (i == j)
          continue;
        const auto& cand_j = pfCands[j];
        if (!cand_j->bestTrack())
          continue;
        const auto trk_j = builder_->build(cand_j);

        data.fillMulti<int>("idx_trk1", i);
        data.fillMulti<int>("idx_trk2", j);

        btagbtvdeep::TrackPairInfoBuilder trkpairinfo;

        trkpairinfo.buildTrackPairInfo(&trk_i,
                                       &trk_j,
                                       vertices->at(0),
                                       -1,  // not used
                                       jetdirection,
                                       IPTools::absoluteImpactParameter3D(trk_i, pv),
                                       IPTools::absoluteTransverseImpactParameter(trk_i, pv));

        data.fillMulti<float>("trkpair_pca_distance", trkpairinfo.pca_distance());
        data.fillMulti<float>("trkpair_pca_significance", trkpairinfo.pca_significance());

        // //   pcas (1 and 2) poistions
        // data.fillMulti<float>("trkpair_pcaSeed_x1", trkpairinfo.pcaSeed_x());
        // data.fillMulti<float>("trkpair_pcaSeed_y1", trkpairinfo.pcaSeed_y());
        // data.fillMulti<float>("trkpair_pcaSeed_z1", trkpairinfo.pcaSeed_z());

        // data.fillMulti<float>("trkpair_pcaSeed_x2", trkpairinfo.pcaTrack_x());
        // data.fillMulti<float>("trkpair_pcaSeed_y2", trkpairinfo.pcaTrack_y());
        // data.fillMulti<float>("trkpair_pcaSeed_z2", trkpairinfo.pcaTrack_z());

        // data.fillMulti<float>("trkpair_pcaSeed_xerr1", trkpairinfo.pcaSeed_xerr());
        // data.fillMulti<float>("trkpair_pcaSeed_yerr1", trkpairinfo.pcaSeed_yerr());
        // data.fillMulti<float>("trkpair_pcaSeed_zerr1", trkpairinfo.pcaSeed_zerr());

        // data.fillMulti<float>("trkpair_pcaSeed_xerr2", trkpairinfo.pcaTrack_xerr());
        // data.fillMulti<float>("trkpair_pcaSeed_yerr2", trkpairinfo.pcaTrack_yerr());
        // data.fillMulti<float>("trkpair_pcaSeed_zerr2", trkpairinfo.pcaTrack_zerr());

        // dot prod betweeen track and pca direction
        data.fillMulti<float>("trkpair_dotprod1", trkpairinfo.dotprodTrack());
        data.fillMulti<float>("trkpair_dotprod2", trkpairinfo.dotprodSeed());

        //pca distance form PV on both tracks
        data.fillMulti<float>("trkpair_pca_dist1", trkpairinfo.pcaSeed_dist());
        data.fillMulti<float>("trkpair_pca_dist2", trkpairinfo.pcaTrack_dist());

        //track track or dir dir dotprod
        data.fillMulti<float>("trkpair_dotprod12_2D", trkpairinfo.dotprodTrackSeed2D());
        data.fillMulti<float>("trkpair_dotprod12_2DV", trkpairinfo.dotprodTrackSeed2DV());
        data.fillMulti<float>("trkpair_dotprod12_3D", trkpairinfo.dotprodTrackSeed3D());
        data.fillMulti<float>("trkpair_dotprod12_3DV", trkpairinfo.dotprodTrackSeed3DV());

        //jet pca relative
        data.fillMulti<float>("trkpair_pca_jetAxis_dist", trkpairinfo.pca_jetAxis_dist());
        data.fillMulti<float>("trkpair_pca_jetAxis_dotprod", trkpairinfo.pca_jetAxis_dotprod());
        // data.fillMulti<float>("trkpair_pca_jetAxis_dEta", trkpairinfo.pca_jetAxis_dEta());
        // data.fillMulti<float>("trkpair_pca_jetAxis_dPhi", trkpairinfo.pca_jetAxis_dPhi());
      }
    }

    return true;
  }

} /* namespace deepntuples */