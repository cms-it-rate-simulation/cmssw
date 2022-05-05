/*
   class Phase2PixelQCoreNtupleOriginal
*/
// DataFormats
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Phase2TrackerDigi/interface/QCore.h"
#include "DataFormats/Phase2TrackerDigi/interface/ROCBitStream.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Geometry
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"

// For ROOT
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>

// STD
#include <memory>
#include <string>
#include <iostream>

// CLHEP (for speed of light)
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"

//#define EDM_ML_DEBUG

using namespace std;
using namespace edm;
using namespace reco;

class Phase2PixelQCoreNtupleOriginal : public edm::one::EDAnalyzer<> {
public:
  explicit Phase2PixelQCoreNtupleOriginal(const edm::ParameterSet& conf);
  virtual ~Phase2PixelQCoreNtupleOriginal();
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& e, const edm::EventSetup& es);

protected:
  void fillEvt(const edm::Event&);

private:
  edm::EDGetTokenT<edm::DetSetVector<QCore> > qcore_token_;
  edm::EDGetTokenT<edm::DetSetVector<ROCBitStream> > bitstream_token_;
  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> pixelDigi_token_;

  //--- Structures for ntupling:
  struct evt {
    int run;
    int evtnum;

    void init();
  } evt_;


  edm::Service<TFileService> tFileService;
  TTree* pixeltree_;
  TTree* pixeltreeOnTrack_;
};


Phase2PixelQCoreNtupleOriginal::Phase2PixelQCoreNtupleOriginal(edm::ParameterSet const& conf)
    : trackerHitAssociatorConfig_(conf, consumesCollector()),
      pixelRecHits_token_(consumes<edmNew::DetSetVector<SiPixelRecHit>>(edm::InputTag("siPixelRecHits"))),

      qcore_token_(consumes<edm::DetSetVector<QCore> >(edm::InputTag("PixelQCore"))),
      bitstream_token_(consumes<edm::DetSetVector<ROCBitStream> >(edm::InputTag("PixelQCore"))),
      pixelDigi_token_(consumes<edm::DetSetVector<PixelDigi>>(conf.getParameter<edm::InputTag>("siPixelDigi"))),
      pixeltree_(0),
      pixeltreeOnTrack_(0) {}

Phase2PixelQCoreNtupleOriginal::~Phase2PixelQCoreNtupleOriginal() {}

void Phase2PixelQCoreNtupleOriginal::endJob() {}

void Phase2PixelQCoreNtupleOriginal::beginJob() {
  pixeltree_ = tFileService->make<TTree>("PixelNtupleOriginal", "Pixel hit analyzer ntuple");
  pixeltreeOnTrack_ = tFileService->make<TTree>("PixelNtupleOriginalOnTrack", "On-Track Pixel hit analyzer ntuple");

}

// Functions that gets called by framework every event
void Phase2PixelQCoreNtupleOriginal::analyze(const edm::Event& e, const edm::EventSetup& es) {

  cout << "In Phase2PixelQCoreNtupleOriginal::analyze" << endl;

  edm::Handle<edm::DetSetVector<QCore> > aQCoreVector;
  e.getByToken(qcore_token_, aQCoreVector);

  edm::Handle<edm::DetSetVector<ROCBitStream> > aBitStreamVector;
  e.getByToken(bitstream_token_, aBitStreamVector);

  std::cout << "RETRIEVED DETSETVECTOR QCORE : " << std::endl;

  edm::DetSetVector<QCore>::const_iterator iterDet;
  for ( iterDet = aQCoreVector->begin();
        iterDet != aQCoreVector->end();
        iterDet++ ) {

    DetId tkId = iterDet->id;

    edm::DetSet<QCore> theQCores = (*aQCoreVector)[ tkId ];

    std::cout << "QCORE DETID NEW : " << tkId.rawId() << std::endl;

    for ( auto iterQCore = theQCores.begin();
          iterQCore != theQCores.end();
          ++iterQCore ) {
      std::cout << "QCORE : " << iterQCore->rocid() << " " << iterQCore->ccol() << " " << iterQCore->qcrow() << std::endl;
    }
  }

  edm::DetSetVector<ROCBitStream>::const_iterator iterDetBitStream;
  for ( iterDetBitStream = aBitStreamVector->begin();
        iterDetBitStream != aBitStreamVector->end();
        iterDetBitStream++ ) {

    DetId tkId = iterDetBitStream->id;

    edm::DetSet<ROCBitStream> theBitStreams = (*aBitStreamVector)[ tkId ];

    std::cout << "BITSTREAM DETID : " << tkId.rawId() << std::endl;

    for ( auto iterBitStream = theBitStreams.begin();
          iterBitStream != theBitStreams.end();
          ++iterBitStream ) {
      std::cout << "BITSTREAM : " << iterBitStream->rocid() << " size = " << iterBitStream->bitstream().size() << std::endl;
    }
  }

  edm::Handle<edm::DetSetVector<PixelDigi> > pixelDigiHandle;
  e.getByToken(pixelDigi_token_, pixelDigiHandle);


  /*

  edm::DetSetVector<PixelDigi>::const_iterator iterDet;
  for ( iterDet = pixelDigiHandle->begin();
        iterDet != pixelDigiHandle->end();
        iterDet++ ) {

    DetId tkId = iterDet->id;

    edm::DetSet<PixelDigi> theDigis = (*pixelDigiHandle)[ tkId ];

    std::vector<std::pair<int,int>> hitlist;
    
    if (tkId.subdetId() == PixelSubdetector::PixelBarrel) {
      int layer_num = tTopo->pxbLayer(tkId.rawId());
      int ladder_num = tTopo->pxbLadder(tkId.rawId());
      int module_num = tTopo->pxbModule(tkId.rawId());
      //cout << "tkID: "<<tkId.subdetId()<<" Layer="<<layer_num<<" Ladder="<<ladder_num<<" Module="<<module_num<<endl;
    } else if (tkId.subdetId() == PixelSubdetector::PixelEndcap) {
      int module_num = tTopo->pxfModule(tkId());
      int disk_num = tTopo->pxfDisk(tkId());
      int blade_num = tTopo->pxfBlade(tkId());
      int panel_num = tTopo->pxfPanel(tkId());
      int side_num = tTopo->pxfSide(tkId());
      //cout << "tkID: "<<tkId.subdetId()<<" Module="<<module_num<<" Disk="<<disk_num<<" Blade="<<blade_num
      //	   <<" Panel="<<panel_num<<" Side="<<side_num<<endl;
    }


    for ( auto iterDigi = theDigis.begin();
          iterDigi != theDigis.end();
          ++iterDigi ) {
      hitlist.emplace_back(iterDigi->row(),iterDigi->column());
    }

    processHits(hitlist);

  }

  */

  if ((recHitColl.product())->dataSize() > 0) {
    std::string detname;

    evt_.init();
    fillEvt(e);

    // Loop over Detector IDs
    for (auto recHitIdIterator : *(recHitColl.product())) {
      SiPixelRecHitCollection::DetSet detset = recHitIdIterator;

      if (detset.empty())
        continue;
      DetId detId = DetId(detset.detId());  // Get the Detid object

      const GeomDet* geomDet(theGeometry->idToDet(detId));

      // Loop over rechits for this detid
      for (auto iterRecHit : detset) {
        // get matched simhit
        matched.clear();
        matched = associate.associateHit(iterRecHit);
        if (!matched.empty()) {
          float closest = 9999.9;
          LocalPoint lp = iterRecHit.localPosition();
          float rechit_x = lp.x();
          float rechit_y = lp.y();
          //loop over simhits and find closest
          for (auto const& m : matched) {
            float sim_x1 = m.entryPoint().x();
            float sim_x2 = m.exitPoint().x();
            float sim_xpos = 0.5 * (sim_x1 + sim_x2);
            float sim_y1 = m.entryPoint().y();
            float sim_y2 = m.exitPoint().y();
            float sim_ypos = 0.5 * (sim_y1 + sim_y2);

            float x_res = sim_xpos - rechit_x;
            float y_res = sim_ypos - rechit_y;
            float dist = sqrt(x_res * x_res + y_res * y_res);
            if (dist < closest) {
              closest = dist;
              closest_simhit = &m;
            }
          }  // end of simhit loop
        }    // end matched emtpy
        unsigned int subid = detId.subdetId();
        int detid_db = detId.rawId();
        int layer_num = -99, ladder_num = -99, module_num = -99, disk_num = -99, blade_num = -99, panel_num = -99,
            side_num = -99;
        if ((subid == PixelSubdetector::PixelBarrel) || (subid == PixelSubdetector::PixelEndcap)) {
          // 1 = PXB, 2 = PXF
          if (subid == PixelSubdetector::PixelBarrel) {
            layer_num = tTopo->pxbLayer(detId.rawId());
            ladder_num = tTopo->pxbLadder(detId.rawId());
            module_num = tTopo->pxbModule(detId.rawId());
#ifdef EDM_ML_DEBUG
            std::cout << "\ndetId = " << subid << " : " << tTopo->pxbLayer(detId.rawId()) << " , "
                      << tTopo->pxbLadder(detId.rawId()) << " , " << tTopo->pxbModule(detId.rawId());
#endif
          } else if (subid == PixelSubdetector::PixelEndcap) {
            module_num = tTopo->pxfModule(detId());
            disk_num = tTopo->pxfDisk(detId());
            blade_num = tTopo->pxfBlade(detId());
            panel_num = tTopo->pxfPanel(detId());
            side_num = tTopo->pxfSide(detId());
          }
          int num_simhit = matched.size();
          recHit_.init();
          // filling in on ALL track rechits
          fillPRecHit(detid_db,
                      subid,
                      layer_num,
                      ladder_num,
                      module_num,
                      disk_num,
                      blade_num,
                      panel_num,
                      side_num,
                      &iterRecHit,
                      num_simhit,
                      closest_simhit,
                      geomDet);
          pixeltree_->Fill();
        }
      }  // end of rechit loop
    }    // end of detid loop
  }      // end of loop test on recHitColl size

  // Now loop over recotracks
  edm::Handle<View<reco::Track>> trackCollection;
  e.getByToken(recoTracks_token_, trackCollection);

  // -- Track trajectory association map
  edm::Handle<TrajTrackAssociationCollection> hTTAC;
  e.getByToken(tta_token_, hTTAC);
  TrajectoryStateCombiner tsoscomb;

  if (!trackCollection.isValid()) {
    if (picky_) {
      throw cms::Exception("ProductNotValid") << "TrackCollection product not valid";
    } else {
      std::cout << "TrackCollection product not valid" << endl;
      ;
    }

  } else if (!hTTAC.isValid()) {
    if (picky_) {
      throw cms::Exception("ProductNotValid") << "TrajectoryAssociationCollection product not valid";
    } else {
      std::cout << "TrajectoryAssociationCollection product not valid" << endl;
    }

  } else {
    int rT = 0;
    const TrajTrackAssociationCollection ttac = *(hTTAC.product());
    for (TrajTrackAssociationCollection::const_iterator it = ttac.begin(); it != ttac.end(); ++it) {
      ++rT;
      const edm::Ref<std::vector<Trajectory>> refTraj = it->key;
      auto track = it->val;
      trkIsHighPurity_ = track->quality(reco::TrackBase::highPurity);
      trkPt_ = track->pt();
      trkEta_ = track->eta();
      trkTheta_ = track->theta();
      trkPhi_ = track->phi();

      int iT = 0;
#ifdef EDM_ML_DEBUG
      std::cout << " num of hits for track " << rT << " = " << track->recHitsSize() << std::endl;
#endif

      std::vector<TrajectoryMeasurement> tmeasColl = refTraj->measurements();
      for (auto const& tmeasIt : tmeasColl) {
        if (!tmeasIt.updatedState().isValid())
          continue;
        if (!tmeasIt.recHit()->isValid())
          continue;

        const TrackingRecHit* hit = tmeasIt.recHit()->hit();
        const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(hit);
        if (pixhit == nullptr)
          continue;
        if (!pixhit->isValid())
          continue;
        ++iT;
        TrajectoryStateOnSurface tsos = tsoscomb(tmeasIt.forwardPredictedState(), tmeasIt.backwardPredictedState());
        const DetId& detId = hit->geographicalId();
        const GeomDet* geomDet(theGeometry->idToDet(detId));

        if (pixhit) {
          // get matched simhit
          matched.clear();
          matched = associate.associateHit(*pixhit);

          if (!matched.empty()) {
            float closest = 9999.9;
            LocalPoint lp = pixhit->localPosition();
            float rechit_x = lp.x();
            float rechit_y = lp.y();

            //loop over simhits and find closest
            //for (std::vector<PSimHit>::const_iterator m = matched.begin(); m<matched.end(); m++)
            for (auto const& m : matched) {
              float sim_x1 = m.entryPoint().x();
              float sim_x2 = m.exitPoint().x();
              float sim_xpos = 0.5 * (sim_x1 + sim_x2);
              float sim_y1 = m.entryPoint().y();
              float sim_y2 = m.exitPoint().y();
              float sim_ypos = 0.5 * (sim_y1 + sim_y2);

              float x_res = sim_xpos - rechit_x;
              float y_res = sim_ypos - rechit_y;
              float dist = sqrt(x_res * x_res + y_res * y_res);
              if (dist < closest) {
                closest = dist;
                closest_simhit = &m;
              }
            }  // end of simhit loop
          }    // end matched emtpy

          int num_simhit = matched.size();

          int layer_num = -99, ladder_num = -99, module_num = -99, disk_num = -99, blade_num = -99, panel_num = -99,
              side_num = -99;

          unsigned int subid = detId.subdetId();
          int detid_db = detId.rawId();
          if ((subid == PixelSubdetector::PixelBarrel) || (subid == PixelSubdetector::PixelEndcap)) {
            // 1 = PXB, 2 = PXF
            if (subid == PixelSubdetector::PixelBarrel) {
              layer_num = tTopo->pxbLayer(detId.rawId());
              ladder_num = tTopo->pxbLadder(detId.rawId());
              module_num = tTopo->pxbModule(detId.rawId());
#ifdef EDM_ML_DEBUG
              std::cout << "\ndetId = " << subid << " : " << tTopo->pxbLayer(detId.rawId()) << " , "
                        << tTopo->pxbLadder(detId.rawId()) << " , " << tTopo->pxbModule(detId.rawId()) << std::endl;
#endif
            } else if (subid == PixelSubdetector::PixelEndcap) {
              module_num = tTopo->pxfModule(detId());
              disk_num = tTopo->pxfDisk(detId());
              blade_num = tTopo->pxfBlade(detId());
              panel_num = tTopo->pxfPanel(detId());
              side_num = tTopo->pxfSide(detId());
            }

            recHit_.init();
            // fill on track rechits
            fillPRecHit(detid_db,
                        subid,
                        layer_num,
                        ladder_num,
                        module_num,
                        disk_num,
                        blade_num,
                        panel_num,
                        side_num,
                        hit,     // TransientTrackingRecHit *
                        pixhit,  // SiPixelRecHit *
                        num_simhit,
                        closest_simhit,
                        geomDet,
                        tsos);
            pixeltreeOnTrack_->Fill();
          }  // if ( (subid==1)||(subid==2) )
        }    // if cast is possible to SiPixelHit
      }      //end of loop on tracking rechits
    }        // end of loop on recotracks
  }          // else track collection is valid
}  // end analyze function

// Function for filling in all the rechits
// I know it is lazy to pass everything, but I'm doing it anyway. -EB
void Phase2PixelQCoreNtupleOriginal::fillPRecHit(const int detid_db,
                                    const int subid,
                                    const int layer_num,
                                    const int ladder_num,
                                    const int module_num,
                                    const int disk_num,
                                    const int blade_num,
                                    const int panel_num,
                                    const int side_num,
                                    SiPixelRecHitCollection::DetSet::const_iterator pixeliter,
                                    const int num_simhit,
                                    const PSimHit* closest_simhit,
                                    const GeomDet* PixGeom) {
  LocalPoint lp = pixeliter->localPosition();
  LocalError le = pixeliter->localPositionError();

  recHit_.x = lp.x();
  recHit_.y = lp.y();
  recHit_.xx = le.xx();
  recHit_.xy = le.xy();
  recHit_.yy = le.yy();
  GlobalPoint GP = PixGeom->surface().toGlobal(lp);
  recHit_.gx = GP.x();
  recHit_.gy = GP.y();
  recHit_.gz = GP.z();
  GlobalPoint GP0 = PixGeom->surface().toGlobal(LocalPoint(0, 0, 0));
  recHit_.theta = GP0.theta();
  recHit_.phi = GP0.phi();

  SiPixelRecHit::ClusterRef const& Cluster = pixeliter->cluster();
  recHit_.q = Cluster->charge();
  recHit_.spreadx = Cluster->sizeX();
  recHit_.spready = Cluster->sizeY();

  recHit_.subid = subid;
  recHit_.nsimhit = num_simhit;

  recHit_.layer = layer_num;
  recHit_.ladder = ladder_num;
  recHit_.module = module_num;
  recHit_.module = module_num;
  recHit_.disk = disk_num;
  recHit_.blade = blade_num;
  recHit_.panel = panel_num;
  recHit_.side = side_num;

  /*-- module topology --*/
  const PixelGeomDetUnit* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(PixGeom);
  const PixelTopology* topol = &(theGeomDet->specificTopology());
  recHit_.nRowsInDet = topol->nrows();
  recHit_.nColsInDet = topol->ncolumns();
  recHit_.pitchx = topol->pitch().first;
  recHit_.pitchy = topol->pitch().second;
  recHit_.thickness = theGeomDet->surface().bounds().thickness();

  MeasurementPoint mp = topol->measurementPosition(LocalPoint(recHit_.x, recHit_.y));
  recHit_.row = mp.x();
  recHit_.col = mp.y();

  if (Cluster.isNonnull()) {
    // compute local angles from det position
    std::pair<float, float> local_angles = computeAnglesFromDetPosition(*Cluster, *topol, *theGeomDet);
    recHit_.cotAlphaFromDet = local_angles.first;
    recHit_.cotBetaFromDet = local_angles.second;

    // -- Get digis of this cluster
    const std::vector<SiPixelCluster::Pixel>& pixvector = Cluster->pixels();
#ifdef EDM_ML_DEBUG
    std::cout << "  Found " << pixvector.size() << " pixels for this cluster " << std::endl;
#endif
    for (unsigned int i = 0; i < pixvector.size(); ++i) {
      if (recHit_.fDgN > DGPERCLMAX - 1)
        break;
      SiPixelCluster::Pixel holdpix = pixvector[i];

      recHit_.fDgRow[recHit_.fDgN] = holdpix.x;
      recHit_.fDgCol[recHit_.fDgN] = holdpix.y;
#ifdef EDM_ML_DEBUG
      std::cout << "holdpix " << holdpix.x << " " << holdpix.y << std::endl;
#endif
      recHit_.fDgDetId[recHit_.fDgN] = detid_db;
      recHit_.fDgAdc[recHit_.fDgN] = -99.;
      recHit_.fDgCharge[recHit_.fDgN] = holdpix.adc / 1000.;
      ++recHit_.fDgN;
    }
  }  // if ( Cluster.isNonnull() )

#ifdef EDM_ML_DEBUG
  std::cout << "num_simhit = " << num_simhit << std::endl;
#endif
  if (num_simhit > 0) {
    recHit_.pdgid = closest_simhit->particleType();
    recHit_.process = closest_simhit->processType();

    float sim_x1 = closest_simhit->entryPoint().x();
    float sim_x2 = closest_simhit->exitPoint().x();
    recHit_.hx = 0.5 * (sim_x1 + sim_x2);
    float sim_y1 = closest_simhit->entryPoint().y();
    float sim_y2 = closest_simhit->exitPoint().y();
    recHit_.hy = 0.5 * (sim_y1 + sim_y2);

    float time_to_detid_ns = GP0.mag() / (CLHEP::c_light * CLHEP::ns / CLHEP::cm);  // speed of light in ns
    recHit_.ht = closest_simhit->timeOfFlight() - time_to_detid_ns;

    recHit_.tx = closest_simhit->localDirection().x();
    recHit_.ty = closest_simhit->localDirection().y();
    recHit_.tz = closest_simhit->localDirection().z();

    MeasurementPoint hmp = topol->measurementPosition(LocalPoint(recHit_.hx, recHit_.hy));
    recHit_.hrow = hmp.x();
    recHit_.hcol = hmp.y();

    // Leaving the comment below, useful for future reference
    // alpha: angle with respect to local x axis in local (x,z) plane
    // float cotalpha = sim_xdir/sim_zdir;
    // beta: angle with respect to local y axis in local (y,z) plane
    // float cotbeta = sim_ydir/sim_zdir;

#ifdef EDM_ML_DEBUG
    std::cout << "num_simhit x, y = " << 0.5 * (sim_x1 + sim_x2) << " " << 0.5 * (sim_y1 + sim_y2) << std::endl;
#endif
  }
#ifdef EDM_ML_DEBUG
  std::cout << "Found RecHit in " << subid
            << " global x/y/z : " << PixGeom->surface().toGlobal(pixeliter->localPosition()).x() << " "
            << PixGeom->surface().toGlobal(pixeliter->localPosition()).y() << " "
            << PixGeom->surface().toGlobal(pixeliter->localPosition()).z() << std::endl;
#endif
}

// Function for filling in on track rechits
void Phase2PixelQCoreNtupleOriginal::fillPRecHit(const int detid_db,
                                    const int subid,
                                    const int layer_num,
                                    const int ladder_num,
                                    const int module_num,
                                    const int disk_num,
                                    const int blade_num,
                                    const int panel_num,
                                    const int side_num,
                                    const TrackingRecHit* recHit,
                                    const SiPixelRecHit* pixHit,
                                    const int num_simhit,
                                    const PSimHit* closest_simhit,
                                    const GeomDet* PixGeom,
                                    const TrajectoryStateOnSurface tsos) {
  LocalPoint lp = recHit->localPosition();
  LocalError le = recHit->localPositionError();

  recHit_.x = lp.x();
  recHit_.y = lp.y();
  recHit_.xx = le.xx();
  recHit_.xy = le.xy();
  recHit_.yy = le.yy();

  recHit_.probQ = pixHit->probabilityQ();
  recHit_.probXY = pixHit->probabilityXY();
  //std::cout << "printing pixHit_.xxloc " << recHit_.xxloc << std::endl;

  GlobalPoint GP = PixGeom->surface().toGlobal(recHit->localPosition());
  recHit_.gx = GP.x();
  recHit_.gy = GP.y();
  recHit_.gz = GP.z();
  GlobalPoint GP0 = PixGeom->surface().toGlobal(LocalPoint(0, 0, 0));
  recHit_.theta = GP0.theta();
  recHit_.phi = GP0.phi();
  recHit_.subid = subid;

  SiPixelRecHit::ClusterRef const& Cluster = pixHit->cluster();
  recHit_.q = Cluster->charge();
  recHit_.spreadx = Cluster->sizeX();
  recHit_.spready = Cluster->sizeY();

  recHit_.nsimhit = num_simhit;

  recHit_.layer = layer_num;
  recHit_.ladder = ladder_num;
  recHit_.module = module_num;
  recHit_.module = module_num;
  recHit_.disk = disk_num;
  recHit_.blade = blade_num;
  recHit_.panel = panel_num;
  recHit_.side = side_num;

  /*-- module topology --*/
  const PixelGeomDetUnit* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(PixGeom);
  const PixelTopology* topol = &(theGeomDet->specificTopology());
  recHit_.nRowsInDet = topol->nrows();
  recHit_.nColsInDet = topol->ncolumns();
  recHit_.pitchx = topol->pitch().first;
  recHit_.pitchy = topol->pitch().second;
  recHit_.thickness = theGeomDet->surface().bounds().thickness();

  if (Cluster.isNonnull()) {
    // compute local angles from det position
    std::pair<float, float> local_angles = computeAnglesFromDetPosition(*Cluster, *topol, *theGeomDet);
    recHit_.cotAlphaFromDet = local_angles.first;
    recHit_.cotBetaFromDet = local_angles.second;

    // compute local angles from track trajectory
    recHit_.cotAlphaFromTrack = tsos.localParameters().dxdz();
    recHit_.cotBetaFromTrack = tsos.localParameters().dydz();

    // -- Get digis of this cluster
    const std::vector<SiPixelCluster::Pixel>& pixvector = Cluster->pixels();
#ifdef EDM_ML_DEBUG
    std::cout << "  Found " << pixvector.size() << " pixels for this cluster " << std::endl;
#endif
    for (unsigned int i = 0; i < pixvector.size(); ++i) {
      if (recHit_.fDgN > DGPERCLMAX - 1)
        break;
      SiPixelCluster::Pixel holdpix = pixvector[i];

      recHit_.fDgRow[recHit_.fDgN] = holdpix.x;
      recHit_.fDgCol[recHit_.fDgN] = holdpix.y;
#ifdef EDM_ML_DEBUG
      std::cout << "holdpix " << holdpix.x << " " << holdpix.y << std::endl;
#endif
      recHit_.fDgDetId[recHit_.fDgN] = detid_db;
      recHit_.fDgAdc[recHit_.fDgN] = -99.;
      recHit_.fDgCharge[recHit_.fDgN] = holdpix.adc / 1000.;
      ++recHit_.fDgN;
    }
  }  // if ( Cluster.isNonnull() )

  if (num_simhit > 0) {
    recHit_.pdgid = closest_simhit->particleType();
    recHit_.process = closest_simhit->processType();

    float sim_x1 = closest_simhit->entryPoint().x();
    float sim_x2 = closest_simhit->exitPoint().x();
    recHit_.hx = 0.5 * (sim_x1 + sim_x2);
    float sim_y1 = closest_simhit->entryPoint().y();
    float sim_y2 = closest_simhit->exitPoint().y();
    recHit_.hy = 0.5 * (sim_y1 + sim_y2);

    float time_to_detid_ns = GP0.mag() / (CLHEP::c_light * CLHEP::ns / CLHEP::cm);  // speed of light in ns
    recHit_.ht = closest_simhit->timeOfFlight() - time_to_detid_ns;

    recHit_.tx = closest_simhit->localDirection().x();
    recHit_.ty = closest_simhit->localDirection().y();
    recHit_.tz = closest_simhit->localDirection().z();

    MeasurementPoint hmp = topol->measurementPosition(LocalPoint(recHit_.hx, recHit_.hy));
    recHit_.hrow = hmp.x();
    recHit_.hcol = hmp.y();

    // Leaving the comment below, useful for future reference
    // alpha: angle with respect to local x axis in local (x,z) plane
    // float cotalpha = sim_xdir/sim_zdir;
    // beta: angle with respect to local y axis in local (y,z) plane
    // float cotbeta = sim_ydir/sim_zdir;

#ifdef EDM_ML_DEBUG
    std::cout << "num_simhit x, y = " << 0.5 * (sim_x1 + sim_x2) << " " << 0.5 * (sim_y1 + sim_y2) << std::endl;
#endif
  }
}


void Phase2PixelQCoreNtupleOriginal::fillEvt(const edm::Event& E) {
  evt_.run = E.id().run();
  evt_.evtnum = E.id().event();
}

void Phase2PixelQCoreNtupleOriginal::evt::init() {
  int dummy_int = 9999;
  run = dummy_int;
  evtnum = dummy_int;
}


void Phase2PixelQCoreNtupleOriginal::RecHit::init() {
  float dummy_float = 9999.0;

  pdgid = 0;
  process = 0;
  q = dummy_float;
  x = dummy_float;
  y = dummy_float;
  xx = dummy_float;
  xy = dummy_float;
  yy = dummy_float;

  row = dummy_float;
  col = dummy_float;
  gx = dummy_float;
  gy = dummy_float;
  gz = dummy_float;
  nsimhit = 0;
  subid = -99;
  module = -99;
  layer = -99;
  ladder = -99;
  disk = -99;
  blade = -99;
  panel = -99;
  side = -99;
  spreadx = 0;
  spready = 0;
  hx = dummy_float;
  hy = dummy_float;
  ht = dummy_float;
  tx = dummy_float;
  ty = dummy_float;
  tz = dummy_float;
  theta = dummy_float;
  phi = dummy_float;

  fDgN = DGPERCLMAX;
  for (int i = 0; i < fDgN; ++i) {
    fDgRow[i] = fDgCol[i] = -9999;
    fDgAdc[i] = fDgCharge[i] = -9999.;
    //    fDgRoc[i] = fDgRocR[i] = fDgRocC[i] = -9999;
  }
  fDgN = 0;
}
std::pair<float, float> Phase2PixelQCoreNtupleOriginal::computeAnglesFromDetPosition(const SiPixelCluster& cl,
                                                                        const PixelTopology& theTopol,
                                                                        const GeomDetUnit& theDet) const {
  // get cluster center of gravity (of charge)
  float xcenter = cl.x();
  float ycenter = cl.y();

  // get the cluster position in local coordinates (cm)

  // ggiurgiu@jhu.edu 12/09/2010 : This function is called without track info, therefore there are no track
  // angles to provide here. Call the default localPosition (without track info)
  LocalPoint lp = theTopol.localPosition(MeasurementPoint(xcenter, ycenter));
  const Local3DPoint origin = theDet.surface().toLocal(GlobalPoint(0, 0, 0));  // can be computed once...

  auto gvx = lp.x() - origin.x();
  auto gvy = lp.y() - origin.y();
  auto gvz = -1.f / origin.z();
  // normalization not required as only ratio used...

  // calculate angles
  float cotalpha_ = gvx * gvz;
  float cotbeta_ = gvy * gvz;

  return std::make_pair(cotalpha_, cotbeta_);
}

void Phase2PixelQCoreNtupleOriginal::processHits(const std::vector<std::pair<int, int> >& hitList) {

  //for(const auto& hit:hitList) {
  //  cout << "row, col : "<<hit.first<<" "<<hit.second<<endl;
  //}

}



//define this as a plug-in
DEFINE_FWK_MODULE(Phase2PixelQCoreNtupleOriginal);
