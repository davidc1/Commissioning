#ifndef LARLITE_WIREHITMAKER_CXX
#define LARLITE_WIREHITMAKER_CXX

#include "WireHitMaker.h"
#include "DataFormat/wire.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  bool WireHitMaker::initialize() {

    // producer name for RawDigits
    _producer = "daq";
    // cut on noise level (only channels with noise less than this value will be searched for hits)
    _rms_max = 3.0;
    // number of ticks to use to calculate baseline & rms
    _nticks = 100;
    // sigma cut on the noise -> once waveform goes above/below this threshold defines the start/end of a hit
    _sigmacut = 4.0;
    // minimum width for a hit (in ticks)
    _minwidth = 4;

    return true;
  }
  
  bool WireHitMaker::analyze(storage_manager* storage) {

    // read in Wires
    auto const ev_wf   = storage->get_data<event_wire>(_producer);
    // create Hits
    auto const ev_hits = storage->get_data<event_hit>("rawhit");

    // set event ID through storage manager
    storage->set_id(storage->get_data<event_wire>(_producer)->run(),
		    storage->get_data<event_wire>(_producer)->subrun(),
		    storage->get_data<event_wire>(_producer)->event_id());

    // loop over waveforms
    for (size_t i=0; i < ev_wf->size(); i++){

      auto const wf = ev_wf->at(i);

      // get vector of ROIs
      auto const& ROIs = wf.SignalROI();
      
      for (auto const& ROI : ROIs.get_ranges()){

	// for each waveform calculate a baseline
	auto const baseline = getBaseline(ROI);
	auto const rms      = getRMS(ROI,baseline);
	
	// if the rms is larger than some max amount, do not search for hits
	if (rms > _rms_max) continue;
	
	// start looking for hits
	auto hits = getHits(ROI,baseline,rms);
	
	// loop over hits and assign the channel information
	for (auto & h : hits){
	  auto ch         = wf.Channel();
	  auto const view = larutil::Geometry::GetME()->View(ch);
	  auto const pl   = larutil::Geometry::GetME()->ChannelToPlane(ch);
	  auto const sigt = larutil::Geometry::GetME()->SignalType(ch);
	  auto wire       = larutil::Geometry::GetME()->ChannelToWireID(ch);
	  
	  h.set_channel(ch);
	  h.set_view(view);
	  h.set_signal_type(sigt);
	  h.set_wire(wire);
	  // push this hit back to the hit object
	  ev_hits->push_back(h);
	}// for all hits that were found
      }// for all ROIs in the wire
    }// loop thrugh wires

    return true;
  }

  bool WireHitMaker::finalize() {

    return true;
  }

  double WireHitMaker::getBaseline(const std::vector<float> &wf){
   
    double baseline = 0.;
    int    ticks    = 0;

    size_t nticks = wf.size();

    for (size_t n=0; n < _nticks; n++){
      if (n >= nticks) break;
      baseline += wf[n];
      ticks += 1;
    }

    baseline /= (double)ticks;

    return baseline;
   }

  double WireHitMaker::getRMS(const std::vector<float> &wf, const double &baseline){
   
    double rms    = 0.;
    int    ticks  = 0;

    size_t nticks = wf.size();

    for (size_t n=0; n < _nticks; n++){
      if (n >= nticks) break;
      rms += (wf[n]-baseline)*(wf[n]-baseline);
      ticks += 1;
    }

    rms = sqrt(rms/ticks);

    return rms;
   }

  std::vector<larlite::hit> WireHitMaker::getHits(const std::vector<float> &wf,
						  const double &baseline,
						  const double &rms){
    
    // prepare a vector of hits to populate
    std::vector<larlite::hit> hits;

    // keep track of:
    double amp;   // amplitude of hit
    double area;  // area of hit
    int    peak;  // peak-time of hit
    int    start; // start time of hit
    int    end;   // end time of hit
    // are we in an active hit?
    bool active = false;

    // loop through waveform
    for (size_t i=0; i < wf.size(); i++){

      double h = wf[i]-baseline;

      // is it above the cut we want to apply?
      if (h > rms*_sigmacut){
	// if not active start the hit
	if (!active){
	  amp = h;
	  area = h;
	  start = i;
	  peak = i;
	  active = true;
	}
	// if already in an active region keep adding to the hit
	else{
	  area += h;
	  if (h > amp){
	    amp = h;
	    peak = i;
	  }
	}
      }
      // else -> we are not in an active region
      else{
	// if we were in an active region we have reached the end of a hit
	if (active){
	  end = i;
	  // now make hit
	  // if passes minimum hit width condition
	  if ( (end-start) > _minwidth){
	    double hiterr = (end-start)/2.;
	    larlite::hit hit;
	    hit.set_time_peak(peak,hiterr);
	    hit.set_time_rms(hiterr);
	    hit.set_amplitude(amp,0.);
	    hit.set_sumq(area);
	    hit.set_integral(area,0.);
	    hits.push_back(hit);
	  }
	  active = false;
	  // clear hit attributes
	  peak = 0;
	  area = 0;
	  end = 0;
	  start = 0;
	  amp = 0;
	}
      }// if not in an active region
    }// looping over waveform

    return hits;
  }

}
#endif
