#ifndef LARLITE_SIGNALPROCAREANORMALIZATION_CXX
#define LARLITE_SIGNALPROCAREANORMALIZATION_CXX

#include "SignalProcAreaNormalization.h"
#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"
#include "DataFormat/track.h"
#include "DataFormat/wire.h"
#include "DataFormat/rawdigit.h"

namespace larlite {

  bool SignalProcAreaNormalization::initialize() {

    _tree = new TTree("tree","Signal Processing Hit Normalization");
    _tree->Branch("_wire",&_wire,"wire/I");
    _tree->Branch("_chan",&_chan,"chan/I");
    _tree->Branch("_raw_area",&_raw_area,"raw_area/D");
    _tree->Branch("_reco_ara",&_reco_area,"reco_area/D");
    _tree->Branch("_tick",&_tick,"tick/I");
    _tree->Branch("_trk_size",&_trk_size,"trk_size/I");
    _tree->Branch("_hit_multiplicity",&_hit_multiplicity,"hit_multiplicity/I");

    return true;
  }
  
  bool SignalProcAreaNormalization::analyze(storage_manager* storage) {

    // load clusters
    auto ev_track = storage->get_data<event_track>("pandoraCosmic");

    event_hit* ev_recohit = nullptr;
    auto const& ass_hit_v = storage->find_one_ass(ev_track->id(), ev_recohit, ev_track->name());

    if ( !ev_recohit or (ev_recohit->size() == 0) ){
      std::cout << "no hits found. exiting" << std::endl;
      return false;
    }

    // load raw hits (produced starting fromraw digits)
    auto ev_rawhit  = storage->get_data<event_hit>("rawhit");

    // make a map from reco hit index to wire number
    std::map<size_t, int> _HitIdx_to_WireNum;
    // make a map from wire number to vector of raw hit indices
    std::map<int, std::vector<size_t> > _WireNum_to_RawHitIdx_v;

    std::cout << "we found " << ass_hit_v.size() << " tracks in this event" << std::endl;

    /*
    // add reco hits to map
    for (size_t i=0; i < ass_hit_v.size(); i++){

      //auto const& trk = ev_track->at(i);
      //if (trk.NumberTrajectoryPoints() < 20)
      //continue;
      //_trk_size = trk.NumberTrajectoryPoints();

      auto ass_hits = ass_hit_v.at(i);
      //if (ass_hits.size() < 20)
      //continue;

      for (size_t j=0; j < ass_hits.size(); j++){

	auto const& hit = ev_recohit->at( ass_hits[j] );

	// ignore non-collection plane
	if (hit.WireID().Plane != 2)
	  continue;

	// get wire number
	auto const& wire = hit.WireID().Wire;
	
	_HitIdx_to_WireNum[i] = wire;
      }// for all hits associated to a large track
    }// for all tracks in the event
    */

    // add reco hits to map
    for (size_t i = 0; i < ev_recohit->size(); i++){

      auto const& hit = ev_recohit->at( i );
      
      // ignore non-collection plane
      if (hit.WireID().Plane != 2)
	continue;
      
      // get wire number
      auto const& wire = hit.WireID().Wire;
      
      _HitIdx_to_WireNum[i] = wire;
    }// for all hits associated to a large track


    // add raw hits to map
    for (size_t i=0; i < ev_rawhit->size(); i++){

      auto const& hit = ev_rawhit->at(i);

      // if plane != collection -> ignore
      if (hit.WireID().Plane != 2)
	continue;

      // get wire number
      auto const& wire = hit.WireID().Wire;

      if (_WireNum_to_RawHitIdx_v.find(wire) == _WireNum_to_RawHitIdx_v.end() ){
	std::vector<size_t> hitidx_v = {i};
	_WireNum_to_RawHitIdx_v[wire] = hitidx_v;
      }
      else
	_WireNum_to_RawHitIdx_v[wire].push_back(i);

    }// for all raw hits -> fill wire map
    
    // for every reconstructed hit, find a matching raw hit in time (if possible)
    // and compare their charge
    typedef std::map<size_t, int>::iterator itt;
    for (itt it = _HitIdx_to_WireNum.begin(); it != _HitIdx_to_WireNum.end(); it++){

      auto const& reco_hit = ev_recohit->at (it->first );

      // get hit time-tick
      _tick = reco_hit.PeakTime();

      // get the wire
      _wire = reco_hit.WireID().Wire;

      // get list of raw hits @ this wire
      auto raw_hits = _WireNum_to_RawHitIdx_v[_wire];

      _hit_multiplicity = raw_hits.size();

      std::cout << "@ wire " << _wire << " found " << raw_hits.size() << " that match reco hit" << std::endl;

      // go through the raw hits, find the one that has ~ the same
      // time tick. This is the correct match -> compare areas
      for (auto& hit_index : raw_hits){

	auto const& raw_hit = ev_rawhit->at(hit_index);

	auto tick = raw_hit.PeakTime();

	if ( ( (tick - _tick) < 6) && ( (tick - _tick) > -6) ){

	// we have a match!
	// compare the areas
	_raw_area  = raw_hit.Integral();
	_reco_area = reco_hit.Integral();

	std::cout << "reco hit time = " << _tick << ", raw hit time = " << tick << "\t reco area = " << _reco_area << "\t raw area = " << _raw_area << std::endl;
	
	_tree->Fill();
	
	}// if we have a match
	
      }// for all raw hits on the same wire
      
    }// for all reconstructed hits

    return true;
  }

  bool SignalProcAreaNormalization::finalize() {

    _tree->Write();

    return true;
  }

}
#endif
