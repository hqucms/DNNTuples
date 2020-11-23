#include <unordered_map>
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DeepNTuples/Ntupler/interface/TrackPairFiller.h"

namespace deepntuples {

void TrackPairFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
}

void TrackPairFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder_);
}

void TrackPairFiller::book() {

  //test
  data.add<int>("n_pfcandidates", 0);
  
  data.addMulti<int>("track1_index");
  data.addMulti<int>("track2_index");
  
  // pca (1 to 2)
  data.addMulti<float>("pca_distance");
  data.addMulti<float>("pca_significance");
  
  //   pcas (1 and 2) poistions
  data.addMulti<float>("pcaSeed_x1");
  data.addMulti<float>("pcaSeed_y1");
  data.addMulti<float>("pcaSeed_z1");
  
  data.addMulti<float>("pcaSeed_x2");
  data.addMulti<float>("pcaSeed_y2");
  data.addMulti<float>("pcaSeed_z2");
  
  data.addMulti<float>("pcaSeed_xerr1");
  data.addMulti<float>("pcaSeed_yerr1");
  data.addMulti<float>("pcaSeed_zerr1");
  
  data.addMulti<float>("pcaSeed_xerr2");
  data.addMulti<float>("pcaSeed_yerr2");
  data.addMulti<float>("pcaSeed_zerr2");
  
  // dot prod betweeen track and pca direction 
  data.addMulti<float>("dotprod1");
  data.addMulti<float>("dotprod2"); 
  
  //pca distance form PV on both tracks
  data.addMulti<float>("pca_dist1");
  data.addMulti<float>("pca_dist2"); 
  
  //track track or dir dir dotprod
  data.addMulti<float>("dotprod12_2D");
  data.addMulti<float>("dotprod12_2DV");
  data.addMulti<float>("dotprod12_3D");
  data.addMulti<float>("dotprod12_3DV");
 
  //jet pca relative
  data.addMulti<float>("pca_jetAxis_dist");
  data.addMulti<float>("pca_jetAxis_dotprod");
  data.addMulti<float>("pca_jetAxis_dEta");
  data.addMulti<float>("pca_jetAxis_dPhi_");
  
  

}

bool TrackPairFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  const auto& pfCands = jet_helper.getJetConstituents();

  data.fill<int>("n_pfcandidates", pfCands.size());
  
  std::vector<reco::TransientTrack> selectedTracks;
  std::vector<int> selectedTracks_pfidx;
  int counter = 0;
  
  
   for (const auto& cand : pfCands){

    //const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));

    if (cand->bestTrack()){
        selectedTracks.push_back( builder_->build(cand) );
        selectedTracks_pfidx.push_back(counter);
        
        
    }
    counter++;
       
    }


    for(std::vector<reco::TransientTrack>::const_iterator it = selectedTracks.begin(); it != selectedTracks.end(); it++){
        
        int index1 = it - selectedTracks.begin();
        
        for(std::vector<reco::TransientTrack>::const_iterator tt = selectedTracks.begin(); tt != selectedTracks.end(); tt++){
            
            int index2 = tt - selectedTracks.begin();
            
            std::cout << index1 << index2 << std::endl;
            
            if (index1!=index2){
            
            data.fillMulti<int>("track1_index", selectedTracks_pfidx[index1]);
            data.fillMulti<int>("track2_index", selectedTracks_pfidx[index2]);  
            
            TrackPairInfoBuilder trkpairinfo;
            
            trkpairinfo.buildTrackPairInfo(&(*it),&(*tt),vertices->at(0),jet);
            
            data.fillMulti<float>("pca_distance", trkpairinfo.pca_distance());
            data.fillMulti<float>("pca_significance", trkpairinfo.pca_significance());

            //   pcas (1 and 2) poistions
            data.fillMulti<float>("pcaSeed_x1", trkpairinfo.pcaSeed_x());
            data.fillMulti<float>("pcaSeed_y1", trkpairinfo.pcaSeed_y());
            data.fillMulti<float>("pcaSeed_z1", trkpairinfo.pcaSeed_z());

            data.fillMulti<float>("pcaSeed_x2", trkpairinfo.pcaTrack_x());
            data.fillMulti<float>("pcaSeed_y2", trkpairinfo.pcaTrack_y());
            data.fillMulti<float>("pcaSeed_z2", trkpairinfo.pcaTrack_z());

            data.fillMulti<float>("pcaSeed_xerr1", trkpairinfo.pcaSeed_xerr());
            data.fillMulti<float>("pcaSeed_yerr1", trkpairinfo.pcaSeed_yerr());
            data.fillMulti<float>("pcaSeed_zerr1", trkpairinfo.pcaSeed_zerr());

            data.fillMulti<float>("pcaSeed_xerr2", trkpairinfo.pcaTrack_xerr());
            data.fillMulti<float>("pcaSeed_yerr2", trkpairinfo.pcaTrack_yerr());
            data.fillMulti<float>("pcaSeed_zerr2", trkpairinfo.pcaTrack_zerr());

            // dot prod betweeen track and pca direction 
            data.fillMulti<float>("dotprod1", trkpairinfo.dotprodTrack());
            data.fillMulti<float>("dotprod2", trkpairinfo.dotprodSeed()); 

            //pca distance form PV on both tracks
            data.fillMulti<float>("pca_dist1", trkpairinfo.pcaSeed_dist());
            data.fillMulti<float>("pca_dist2", trkpairinfo.pcaTrack_dist()); 

            //track track or dir dir dotprod
            data.fillMulti<float>("dotprod12_2D", trkpairinfo.dotprodTrackSeed2D());
            data.fillMulti<float>("dotprod12_2DV", trkpairinfo.dotprodTrackSeed2DV());
            data.fillMulti<float>("dotprod12_3D", trkpairinfo.dotprodTrackSeed3D());
            data.fillMulti<float>("dotprod12_3DV", trkpairinfo.dotprodTrackSeed3DV());

            //jet pca relative
            data.fillMulti<float>("pca_jetAxis_dist", trkpairinfo.pca_jetAxis_dist());
            data.fillMulti<float>("pca_jetAxis_dotprod", trkpairinfo.pca_jetAxis_dotprod());
            data.fillMulti<float>("pca_jetAxis_dEta", trkpairinfo.pca_jetAxis_dEta());
            data.fillMulti<float>("pca_jetAxis_dPhi_", trkpairinfo.pca_jetAxis_dPhi());
            
            }
            
            
        }
        
        
    }


  return true;
}

} /* namespace deepntuples */
