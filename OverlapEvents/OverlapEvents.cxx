#ifndef LARLITE_OVERLAPEVENTS_CXX
#define LARLITE_OVERLAPEVENTS_CXX

#include "OverlapEvents.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  bool OverlapEvents::initialize() {

    std::cout << "aaa" << std::endl;

    for (size_t i=0; i < 3456; i++)
      _laser_chans.push_back(std::vector<float>(9600,0.));

    _evt = 0;

    return true;
  }
  
  bool OverlapEvents::analyze(storage_manager* storage) {

    // read in RawDigits
    auto const ev_wf     = storage->get_data<event_rawdigit>("daq");



    // loop over waveforms
    for (size_t i=0; i < ev_wf->size(); i++){
      
      auto const wf   = ev_wf->at(i);
      auto const ADCs = wf.ADCs();
      auto const ch   = wf.Channel();
      auto const pl   = larutil::Geometry::GetME()->ChannelToPlane(ch);
      auto const wire = larutil::Geometry::GetME()->ChannelToWireID(ch);

      if (pl != 2) continue;

      for (size_t t=0; t < ADCs.size(); t++)
	_laser_chans.at(wire.Wire).at(t) += ADCs[t];

    }

    _evt += 1;
    
    if (_evt >= 180){
      auto const ev_laser  = storage->get_data<event_rawdigit>("laser");
      // set event ID through storage manager
      storage->set_id(storage->get_data<event_rawdigit>("daq")->run(),
		      storage->get_data<event_rawdigit>("daq")->subrun(),
		      storage->get_data<event_rawdigit>("daq")->event_id()-179);
      for (size_t i=0; i < ev_wf->size(); i++){
	auto const wf2   = ev_wf->at(i);
	auto ch   = wf2.Channel();
	auto pl   = larutil::Geometry::GetME()->ChannelToPlane(ch);
	auto wire = larutil::Geometry::GetME()->ChannelToWireID(ch);
	auto comp = wf2.Compression();
	auto samples = wf2.ADCs().size();
	std::vector<short> adcs;
	for (size_t s = 0; s < samples; s++)
	  adcs.push_back((short)(_laser_chans.at(wire.Wire).at(s)/_evt));
	larlite::rawdigit wf3(ch,samples,adcs,comp);
	ev_laser->push_back(wf3);
      }// for all wfs
    }
    return true;
  }

  bool OverlapEvents::finalize() {

    return true;
  }

}
#endif
