// -*- C++ -*-
//
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

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetId.h"

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

  //produces<int>("integer").setBranchAlias( "integer" );
  produces<edm::DetSetVector<int> >();

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
  unique_ptr<edm::DetSetVector<int> > aIntVector = make_unique<edm::DetSetVector<int> >();
  const DetSet<int> one(1);
  const DetSet<int> three(3);
  const DetSet<int> seven(7);

  aIntVector->insert(one);
  aIntVector->insert(three);
  aIntVector->insert(seven);

  //auto_ptr<PointCollection> outerPoints( new PointCollection );
  // and already reserve some space for the new data, to control the size
  // of your executible's memory use.

  //const int size = tracks->size();
  //innerPoints->reserve( size );
  //outerPoints->reserve( size );
  // loop over the tracks:
  //for( TrackCollection::const_iterator track = tracks->begin(); 
  //     track != tracks->end(); ++ track ) {
    // fill the points in the vectors
  //  innerPoints->push_back( track->innerPosition() );
  //  outerPoints->push_back( track->outerPosition() );
  //}
  // and save the vectors


  std::cout << "WILL STORE INT DETSETVECTOR IN EVENT" <<std::endl;

  iEvent.put( std::move(aIntVector) );
  //iEvent.put( std::move(aInteger), "integer" );

  //iEvent.put( outerPoints, "outerPoint" );
   
/* This is an event example
  //Read 'ExampleData' from the Event
  ExampleData const& in = iEvent.get(inToken_);

  //Use the ExampleData to create an ExampleData2 which 
  // is put into the Event
  iEvent.put(std::make_unique<ExampleData2>(in));
*/

/* this is an EventSetup example
  //Read SetupData from the SetupRecord in the EventSetup
  SetupData& setup = iSetup.getData(setupToken_);
*/
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
