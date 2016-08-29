#ifndef LARLITE_ELECGAINEXTRACTION_CXX
#define LARLITE_ELECGAINEXTRACTION_CXX

#include "ElecGainExtraction.h"
#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"
#include "DataFormat/track.h"
#include "DataFormat/wire.h"
#include "DataFormat/simch.h"
#include "DataFormat/rawdigit.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  bool ElecGainExtraction::initialize() {

    _tree = new TTree("tree","Signal Processing Hit Normalization");
    _tree->Branch("_pl",&_pl,"pl/I");
    _tree->Branch("_wire",&_wire,"wire/I");
    _tree->Branch("_chan",&_chan,"chan/I");
    _tree->Branch("_reco_area",&_reco_area,"reco_area/D");
    _tree->Branch("_reco_ampl",&_reco_ampl,"reco_ampl/D");
    _tree->Branch("_q",&_q,"q/D");
    _tree->Branch("_tick",&_tick,"tick/I");

    return true;
  }
  
  bool ElecGainExtraction::analyze(storage_manager* storage) {

    auto geom    = ::larutil::Geometry::GetME();


    auto const& ev_clus = storage->get_data<event_cluster>(_clus_producer);

    if ( !ev_clus or (ev_clus->size() == 0) ){
      std::cout << "no hits found. exiting" << std::endl;
      return false;
    }
    
    event_hit* ev_recohit = nullptr;
    auto const& ass_hit_v = storage->find_one_ass(ev_clus->id(), ev_recohit, ev_clus->name());
    
    if ( !ev_recohit or (ev_recohit->size() == 0) ){
      std::cout << "no hits found. exiting" << std::endl;
      return false;
    }

    /*
    auto const& ev_hit = storage->get_data<event_hit>(_hit_producer);

    if ( !ev_hit or (ev_hit->size() == 0) ){
      std::cout << "no hits found. exiting" << std::endl;
      return false;
    }
    */

    // load simchannels
    auto ev_simch   = storage->get_data<event_simch>("largeant");
    // make a map from channel to position in ev_simch
    std::map<unsigned int, size_t> Ch2SimchIdx_map;
    for (size_t i=0; i < ev_simch->size(); i++)
      Ch2SimchIdx_map[ ev_simch->at(i).Channel() ] = i;
    
    // make a map from reco hit index to wire number
    std::map<size_t, int> _HitIdx_to_WireNum;
    // make a map from wire number to vector of raw hit indices
    std::map<int, std::vector<size_t> > _WireNum_to_RawHitIdx_v;

    //    for (size_t i=0; i < ev_hit->size(); i++){
    
    for (auto const& ass_clus_hit : ass_hit_v){

      if (ass_clus_hit.size() < 100)
	continue;

      for (auto const& i : ass_clus_hit){
	
	auto const& reco_hit = ev_recohit->at ( i );
	
	// get hit time-tick
	_tick = reco_hit.PeakTime();
	
	// get the wire
	_wire = reco_hit.WireID().Wire;
	// and the plane
	_pl   = reco_hit.WireID().Plane;
	
	_reco_area = reco_hit.Integral();
	_reco_ampl = reco_hit.PeakAmplitude();
	
	// find the simchannel info for this channel
	auto chan = geom->PlaneWireToChannel(_pl,_wire);
	
	// are there simch @ this channel?
	if (Ch2SimchIdx_map.find(chan) == Ch2SimchIdx_map.end())
	  continue;
	
	// find the simch associated with this 
	// vector of LArLite IDEs
	auto const& simch = ev_simch->at( Ch2SimchIdx_map[chan] );
	
	std::cout << "Hit time : " << reco_hit.PeakTime() << std::endl;
	
	// get the charge deposited in the appropriate time-range
	auto const& ide_v = simch.TrackIDsAndEnergies( reco_hit.PeakTime() - 4*reco_hit.RMS() + 7298,
						       reco_hit.PeakTime() + 4*reco_hit.RMS() + 7298);
	
	_q = 0;
	for (auto const&  ide : ide_v)
	  _q += ide.numElectrons;
	
	_tree->Fill();
	
      }// for all reconstructed hits
    }// for all clusters

    return true;
  }

  bool ElecGainExtraction::finalize() {

    _tree->Write();

    return true;
  }

}
#endif
