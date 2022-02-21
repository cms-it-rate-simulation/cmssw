/*
   class Phase2PixelQCoreNtupleAltered
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

class Phase2PixelQCoreNtupleAltered : public edm::one::EDAnalyzer<> {
public:
  explicit Phase2PixelQCoreNtupleAltered(const edm::ParameterSet& conf);
  virtual ~Phase2PixelQCoreNtupleAltered();
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

Phase2PixelQCoreNtupleAltered::Phase2PixelQCoreNtupleAltered(edm::ParameterSet const& conf)
  :
      qcore_token_(consumes<edm::DetSetVector<QCore> >(edm::InputTag("PixelQCore"))),
      bitstream_token_(consumes<edm::DetSetVector<ROCBitStream> >(edm::InputTag("PixelQCore"))),
      pixelDigi_token_(consumes<edm::DetSetVector<PixelDigi>>(conf.getParameter<edm::InputTag>("siPixelDigi"))),
      pixeltree_(0),
      pixeltreeOnTrack_(0) {}

Phase2PixelQCoreNtupleAltered::~Phase2PixelQCoreNtupleAltered() {}

void Phase2PixelQCoreNtupleAltered::endJob() {}

void Phase2PixelQCoreNtupleAltered::beginJob() {
  pixeltree_ = tFileService->make<TTree>("PixelNtuple", "Pixel hit analyzer ntuple");
  pixeltreeOnTrack_ = tFileService->make<TTree>("PixelNtupleOnTrack", "On-Track Pixel hit analyzer ntuple");

}

// Functions that gets called by framework every event
void Phase2PixelQCoreNtupleAltered::analyze(const edm::Event& e, const edm::EventSetup& es) {

  cout << "In Phase2PixelQCoreNtupleAltered::analyze" << endl;

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
      std::cout << "BITSTREAM : " << iterBitStream->get_rocid() << " size = " << iterBitStream->get_bitstream().size() << std::endl;
    }
  }

  edm::Handle<edm::DetSetVector<PixelDigi> > pixelDigiHandle;
  e.getByToken(pixelDigi_token_, pixelDigiHandle);
}


void Phase2PixelQCoreNtupleAltered::fillEvt(const edm::Event& E) {
  evt_.run = E.id().run();
  evt_.evtnum = E.id().event();
}

void Phase2PixelQCoreNtupleAltered::evt::init() {
  int dummy_int = 9999;
  run = dummy_int;
  evtnum = dummy_int;
}


//define this as a plug-in
DEFINE_FWK_MODULE(Phase2PixelQCoreNtupleAltered);
