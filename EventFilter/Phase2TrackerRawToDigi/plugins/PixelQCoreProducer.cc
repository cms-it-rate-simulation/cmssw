// -*- C++ -*-
// hi
// Package:    EventFilter/Phase2TrackerRawToDigi
// Class:      PixelQCoreProducer
//

/*
 Description: Make QCore objects for digis

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rohan Misra
//         Created:  Thu, 30 Sep 2021 02:04:06 GMT
//
//

// system include files
#include <memory>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Math/interface/Point3D.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/transform.h"

#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "DataFormats/Phase2TrackerDigi/interface/Hit.h"
#include "DataFormats/Phase2TrackerDigi/interface/QCore.h"
#include "DataFormats/Phase2TrackerDigi/interface/ReadoutChip.h"
#include "DataFormats/Phase2TrackerDigi/interface/ROCBitStream.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"


//
// class declaration
//

class PixelQCoreProducer : public edm::stream::EDProducer<> {
public:
  explicit PixelQCoreProducer(const edm::ParameterSet&);
  ~PixelQCoreProducer();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  
  edm::InputTag src_;
  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> pixelDigi_token_;
  typedef math::XYZPointD Point;
  typedef std::vector<Point> PointCollection;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
PixelQCoreProducer::PixelQCoreProducer(const edm::ParameterSet& iConfig)
{
  src_  = iConfig.getParameter<edm::InputTag>( "src" );
  //produces<PointCollection>( "innerPoint" ).setBranchAlias( "innerPoints");
  //produces<PointCollection>( "outerPoint" ).setBranchAlias( "outerPoints");
  //register your products
  
  pixelDigi_token_ = consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("siPixelDigi"));

  //produces<int>("qcores").setBranchAlias( "qcores" );
  produces<edm::DetSetVector<QCore> >();
  produces<edm::DetSetVector<ROCBitStream> >();

/* Examples
  produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
*/
  //now do what ever other initialization is needed
}

PixelQCoreProducer::~PixelQCoreProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

std::vector<Hit> adjustEdges(std::vector<Hit> hitList) {
        std::vector<Hit> hitListOne = {};
        std::vector<Hit> hitListTwo = {};

        for(auto hit:hitList) {
                if(hit.row() < 672) {
                        hitListOne.push_back(hit);
                } else if(672 <= hit.row() && hit.row() <= 675) {
                        hitListOne.push_back(Hit(671, hit.col(), hit.adc()));
                } else if(676 <= hit.row() && hit.row() <= 679) {
                        hitListOne.push_back(Hit(672, hit.col(), hit.adc()));
                } else if(hit.row() > 679) {
                        hitListOne.push_back(Hit(hit.row() - 8, hit.col(), hit.adc()));
                }
        }

        for(auto hit:hitListOne) {
                if(hit.col() < 216) {
                        hitListTwo.push_back(hit);
                } else if(hit.col() == 216) {
                        hitListTwo.push_back(Hit(hit.row(), 215, hit.adc()));
                } else if(hit.col() == 217) {
                        hitListTwo.push_back(Hit(hit.row(), 216, hit.adc()));
                } else if(hit.col() > 217) {
                        hitListTwo.push_back(Hit(hit.row(), hit.col() - 2, hit.adc()));
                }
        }

        return hitListTwo;
}

std::vector<ReadoutChip> splitByChip(std::vector<Hit> hitList) {
        std::vector<Hit> chip1 = {};
        std::vector<Hit> chip2 = {};
        std::vector<Hit> chip3 = {};
        std::vector<Hit> chip4 = {};

        for(auto hit:hitList) {
                if(hit.row() < 672) {
                        if(hit.col() < 216) {
                                chip1.push_back(hit);
                        } else {
                                chip2.push_back(hit);
                        }
                } else {
                        if(hit.col() < 216) {
                                chip3.push_back(hit);
                        } else {
                                chip4.push_back(hit);
                        }
                }
        }

        return {ReadoutChip(0,chip1), ReadoutChip(1, chip2), ReadoutChip(2, chip3),ReadoutChip(3, chip4)};
}

std::vector<ReadoutChip> processHits(std::vector<Hit> hitList) {
        std::vector<Hit> newHitList;

        std::cout << "Hits:" << "\n";
        for(auto& hit:hitList) {
                std::cout << "row, col : " << hit.row() << ", " << hit.col() << "\n";
        }

        newHitList = adjustEdges(hitList);

        std::vector<ReadoutChip> chips = splitByChip(newHitList);
        std::vector<bool> code;

        for(size_t i = 0; i < chips.size(); i++) {
                ReadoutChip chip = chips[i];
                code = chip.get_chip_code();

                std::cout << "number of hits: " << chip.size() << "\n";
                std::cout << "code length: " << code.size() << "\n";
                std::cout << "chip code: ";
                for(size_t j = 0; j < code.size(); j++) {
                        std::cout << code[j];
                }
                std::cout << "\n";
        }
        std::cout << "\n";
	return chips;
}

// ------------ method called to produce the data  ------------
void PixelQCoreProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm; 
  //using namespace reco; 
  using namespace std;
  // retrieve the tracks
  //Handle<TrackCollection> tracks;
  //iEvent.getByLabel( src_, tracks );
  // create the vectors. Use auto_ptr, as these pointers will automatically
  // delete when they go out of scope, a very efficient way to reduce memory leaks.
  
  unique_ptr<edm::DetSetVector<QCore> > aQCoreVector = make_unique<edm::DetSetVector<QCore> >();
  unique_ptr<edm::DetSetVector<ROCBitStream> > aBitStreamVector = make_unique<edm::DetSetVector<ROCBitStream> >();


  //Retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  const TrackerTopology* const tTopo = tTopoHandle.product();


  edm::Handle<edm::DetSetVector<PixelDigi> > pixelDigiHandle;
  iEvent.getByToken(pixelDigi_token_, pixelDigiHandle);

  edm::DetSetVector<PixelDigi>::const_iterator iterDet;
  for ( iterDet = pixelDigiHandle->begin();
        iterDet != pixelDigiHandle->end();
        iterDet++ ) {

    DetId tkId = iterDet->id;

    edm::DetSet<PixelDigi> theDigis = (*pixelDigiHandle)[ tkId ];



    //break;

    std::vector<Hit> hitlist;

    std::vector<int> id;

    if (tkId.subdetId() == PixelSubdetector::PixelBarrel) {
        int layer_num = tTopo->pxbLayer(tkId.rawId());
        int ladder_num = tTopo->pxbLadder(tkId.rawId());
        int module_num = tTopo->pxbModule(tkId.rawId());
        id = {tkId.subdetId(), layer_num, ladder_num, module_num};
        cout << "tkID: "<<tkId.subdetId()<<" Layer="<<layer_num<<" Ladder="<<ladder_num<<" Module="<<module_num<<endl;
    } else if (tkId.subdetId() == PixelSubdetector::PixelEndcap) {
        int module_num = tTopo->pxfModule(tkId());
        int disk_num = tTopo->pxfDisk(tkId());
        int blade_num = tTopo->pxfBlade(tkId());
        int panel_num = tTopo->pxfPanel(tkId());
        int side_num = tTopo->pxfSide(tkId());
        id = {tkId.subdetId(), module_num, disk_num, blade_num, panel_num, side_num};
        cout << "tkID: "<<tkId.subdetId()<<" Module="<<module_num<<" Disk="<<disk_num<<" Blade="<<blade_num<<" Panel="<<panel_num<<" Side="<<side_num<<endl;
    }

    for ( auto iterDigi = theDigis.begin();
          iterDigi != theDigis.end();
          ++iterDigi ) {
      hitlist.emplace_back(Hit(iterDigi->row(),iterDigi->column(),iterDigi->adc()));
    }

    std::vector<ReadoutChip> chips = processHits(hitlist);

    std::cout << "Make DetSet" << std::endl;

    DetSet<QCore> DetSetQCores(tkId);
    DetSet<ROCBitStream> DetSetBitStream(tkId);

    for(size_t i = 0; i < chips.size(); i++) {

      std::cout << "Retrieve chip " << i << std::endl;

      ReadoutChip chip = chips[i];

      std::cout << "Got chip " << i << std::endl;
      
      std::vector<QCore> qcores = chip.getOrganizedQCores();

      std::cout << "Got qcores " << std::endl;
      
      for (auto& qcore:qcores) {
	std::cout << "push qcore" << std::endl;
	DetSetQCores.push_back(qcore);
      }

      ROCBitStream aROCBitStream(i,chip.get_chip_code());

      DetSetBitStream.push_back(aROCBitStream);

      std::cout << "Done processing chip" << std::endl;
    }  

    aBitStreamVector->insert(DetSetBitStream);
    std::cout << "Add DetSetQCores to DetSetVector" << std::endl;
    aQCoreVector->insert(DetSetQCores);
    std::cout << "Done adding DetSetQCores to DetSetVector" << std::endl;
  }

  std::cout << "WILL STORE QCORE DETSETVECTOR IN EVENT" <<std::endl;
    
  iEvent.put( std::move(aQCoreVector) );
  iEvent.put( std::move(aBitStreamVector) );

  std::cout << "DONE STORE QCORE DETSETVECTOR IN EVENT" <<std::endl;
  
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
//void PixelQCoreProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
//}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
//void PixelQCoreProducer::endStream() {
  // please remove this method if not needed
//}

// ------------ method called when starting to processes a run  ------------
/*
void
PixelQCoreProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a job  ------------
void
PixelQCoreProducer::beginJob(edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a run  ------------
/*
void
PixelQCoreProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a job  ------------
void
PixelQCoreProducer::endJob()
{
}


// ------------ method called when starting to processes a luminosity block  ------------
/*
void
PixelQCoreProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
PixelQCoreProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//void PixelQCoreProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
//  edm::ParameterSetDescription desc;
//  desc.setUnknown();
//  descriptions.addDefault(desc);
//}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelQCoreProducer);
