#ifndef LARLITE_MAKEHITS_CXX
#define LARLITE_MAKEHITS_CXX

#include "MakeHits.h"
#include "DataFormat/hit.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  bool MakeHits::initialize() {

    _producer = "gaushit";

    return true;
  }
  
  bool MakeHits::analyze(storage_manager* storage) {

    auto evt_hits = storage->get_data<event_hit>(_producer);
    auto new_hits = storage->get_data<event_hit>("testhit");
    
    storage->set_id(evt_hits->run(),
		    evt_hits->subrun(),
		    evt_hits->event_id());

    std::vector<int> wirenums;//{1,3,8,9,2,4,6,7,10,5};
    std::vector<int> timeticks;//{1,3,8,9,2,4,6,7,10,5};
    int m = 20;
    for (int i=1; i < 100; i++){
      int remainder = i%m;
      int division  = i/m;
      int result    = m+division*m-remainder;
      std::cout << "result is " << result << std::endl;
      wirenums.push_back(result);
      timeticks.push_back(result);
      //wirenums.push_back(i);
      //timeticks.push_back(i);

    }
    
    for (size_t i=0; i < wirenums.size(); i++){

      
      larlite::hit h(evt_hits->at(i));
      /*
      std::cout << "Hit info           : " << std::endl
		<< "Start Tick         : " << h.StartTick() << std::endl
		<< "End Tick           : " << h.EndTick() << std::endl
		<< "Peak Time          : " << h.PeakTime() << std::endl
		<< "SigmaPeakTime      : " << h.SigmaPeakTime() << std::endl
		<< "RMS                : " << h.RMS() << std::endl
		<< "PeakAmplitude      : " << h.PeakAmplitude() << std::endl
		<< "SigmaPeakAmplitude : " << h.SigmaPeakAmplitude() << std::endl
		<< "SummedADC          : " << h.SummedADC() << std::endl
		<< "Integral           : " << h.Integral() << std::endl
		<< "SigmaIntegral      : " << h.SigmaIntegral() << std::endl
		<< "Multiplicity       : " << h.Multiplicity() << std::endl
		<< "LocalIndex         : " << h.LocalIndex() << std::endl
		<< "GoodnessOfFit      : " << h.GoodnessOfFit() << std::endl
		<< "DegreesOfFreedom   : " << h.DegreesOfFreedom() << std::endl
		<< "Channel            : " << h.Channel() << std::endl
		<< "View               : " << h.View() << std::endl
		<< "SignalType         : " << h.SignalType() << std::endl << std::endl;
      */
      h.set_time_peak(3*timeticks[i],2);
      h.set_time_range(3*timeticks[i]-1,3*timeticks[i]+1);
      h.set_time_rms(1.0);
      //h.set_amplitude(100,10);
      //h.set_integral(1000,100);
      //h.set_sumq(1000);
      //h.set_rms(1);
      //h.set_multiplicity(1);
      larlite::geo::WireID wire(1,1,2,wirenums[i]);
      h.set_wire(wire);
      h.set_signal_type(larlite::geo::SigType_t::kCollection);
      h.set_view(larlite::geo::View_t::kW);
      h.set_channel(larutil::Geometry::GetME()->PlaneWireToChannel(2,wirenums[i]));
      new_hits->push_back(h);
    }// for all hits that we want to create
  
    return true;
  }

  bool MakeHits::finalize() {

  
    return true;
  }

}
#endif
