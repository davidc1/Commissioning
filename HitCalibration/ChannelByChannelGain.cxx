#ifndef LARLITE_CHANNELBYCHANNELGAIN_CXX
#define LARLITE_CHANNELBYCHANNELGAIN_CXX

#include "ChannelByChannelGain.h"
#include "DataFormat/hit.h"

namespace larlite {

  bool ChannelByChannelGain::initialize() {

    if (_tree) { delete _tree; }

    _tree = new TTree("_tree","tree");
    _tree->Branch("_ch",&_ch,"ch/I");
    _tree->Branch("_pl",&_pl,"pl/I");
    _tree->Branch("_run",&_run,"run/I");
    _tree->Branch("_subrun",&_subrun,"subrun/I");
    _tree->Branch("_event",&_event,"event/I");
    _tree->Branch("_area",&_area,"area/D");
    _tree->Branch("_chi",&_chi,"chi/D");

    return true;
  }
  
  bool ChannelByChannelGain::analyze(storage_manager* storage) {

    _run    = storage->run_id();
    _subrun = storage->subrun_id();
    _event  = storage->event_id();

    auto *ev_hit = storage->get_data<event_hit>(_hit_producer);
    
    if (!ev_hit or (ev_hit->size() == 0) ){
      return false;
    }

    for (size_t i =0; i < ev_hit->size(); i++){

      auto const& hit = ev_hit->at(i);

      if ( (hit.WireID().Wire < 3000) or (hit.WireID().Wire > 3100) )
	continue;

      if ( (hit.Multiplicity() != 1) or (hit.GoodnessOfFit() > 1) )
	continue;

      _area = hit.Integral();
      _chi  = hit.GoodnessOfFit();
      _ch   = hit.WireID().Wire;
      _pl   = hit.WireID().Plane;

      _tree->Fill();
      
    }
  
    return true;
  }

  bool ChannelByChannelGain::finalize() {

    _tree->Write();

    return true;
  }

}
#endif
