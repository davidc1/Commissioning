#ifndef LARLITE_SIGNALPROCAREANORMALIZATION_CXX
#define LARLITE_SIGNALPROCAREANORMALIZATION_CXX

#include "SignalProcAreaNormalization.h"
#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"
#include "DataFormat/track.h"
#include "DataFormat/wire.h"
#include "DataFormat/simch.h"
#include "DataFormat/rawdigit.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  bool SignalProcAreaNormalization::initialize() {

    _tree = new TTree("tree","Signal Processing Hit Normalization");
    _tree->Branch("_pl",&_pl,"pl/I");
    _tree->Branch("_wire",&_wire,"wire/I");
    _tree->Branch("_chan",&_chan,"chan/I");
    _tree->Branch("_raw_area",&_raw_area,"raw_area/D");
    _tree->Branch("_reco_ara",&_reco_area,"reco_area/D");
    _tree->Branch("_q",&_q,"q/D");
    _tree->Branch("_tick",&_tick,"tick/I");
    _tree->Branch("_trk_size",&_trk_size,"trk_size/I");
    _tree->Branch("_hit_multiplicity",&_hit_multiplicity,"hit_multiplicity/I");

    return true;
  }
  
  bool SignalProcAreaNormalization::analyze(storage_manager* storage) {

    auto geom    = ::larutil::Geometry::GetME();

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

    /*
    // load simchannels
    auto ev_simch   = storage->get_data<event_simch>("largeant");
    // make a map from channel to position in ev_simch
    std::map<unsigned int, size_t> Ch2SimchIdx_map;
    for (size_t i=0; i < ev_simch->size(); i++)
      Ch2SimchIdx_map[ ev_simch->at(i).Channel() ] = i;
    */
    
    // make a map from reco hit index to wire number
    std::map<size_t, int> _HitIdx_to_WireNum;
    // make a map from wire number to vector of raw hit indices
    std::map<int, std::vector<size_t> > _WireNum_to_RawHitIdx_v;

    std::cout << "we found " << ass_hit_v.size() << " tracks in this event" << std::endl;



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

      std::cout << "tick " << _tick << std::endl;

      // get the wire
      _wire = reco_hit.WireID().Wire;
      // and the plane
      _pl   = reco_hit.WireID().Plane;

      // get list of raw hits @ this wire
      auto raw_hits = _WireNum_to_RawHitIdx_v[_wire];

      _hit_multiplicity = raw_hits.size();

      //std::cout << "@ wire " << _wire << " found " << raw_hits.size() << " that match reco hit" << std::endl;

      // go through the raw hits, find the one that has ~ the same
      // time tick. This is the correct match -> compare areas
      for (auto& hit_index : raw_hits){

	auto const& raw_hit = ev_rawhit->at(hit_index);

	auto tick = raw_hit.PeakTime() - 2405;

	std::cout << "\traw tick : " << tick << std::endl;

	if ( ( (tick - _tick) < 6) && ( (tick - _tick) > -6) ){

	// we have a match!
	// compare the areas
	_raw_area  = raw_hit.Integral();
	_reco_area = reco_hit.Integral();

	std::cout << "reco hit time = " << _tick << ", raw hit time = " << tick << "\t reco area = " << _reco_area << "\t raw area = " << _raw_area << std::endl;

	/*
	// find the simchannel info for this channel
	auto chan = geom->PlaneWireToChannel(_pl,_wire);

	// find the simch associated with this 
	// vector of LArLite IDEs
	auto const& simch = ev_simch->at( Ch2SimchIdx_map[chan] );

	// get the charge deposited in the appropriate time-range
	auto const& ide_v = simch.TrackIDsAndEnergies( reco_hit.PeakTime() - 3*reco_hit.RMS()+3050,
						       reco_hit.PeakTime() + 3*reco_hit.RMS()+3050);

	_q = 0;
	for (auto const&  ide : ide_v)
	  _q += ide.numElectrons;
	*/
	
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
