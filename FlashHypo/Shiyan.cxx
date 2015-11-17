#ifndef LARLITE_SHIYAN_CXX
#define LARLITE_SHIYAN_CXX

#include "Shiyan.h"
#include "FlashHypo.h"

namespace larlite {

  bool Shiyan::initialize() {

    if (_tree) {delete _tree;}
    _tree = new TTree("Paddletree", "Paddle Tree");
    _tree->Branch("PE","std::vector<double>",&_pe);
    
    return true;
  }
  
  bool Shiyan::analyze(storage_manager* storage) {
  
    auto ev_reco = storage->get_data<event_track>("trackkalmanhit");
    if(!ev_reco)
      {std::cout<<"wait.....\n";}
    for(size_t i = 0;i <ev_reco->size();i++ ){
      auto const& trk = ev_reco->at(i);
      if (trk.NumberTrajectoryPoints() > 1) {
	::geoalgo::Trajectory trj;
	for (size_t pt = 0; pt < trk.NumberTrajectoryPoints(); pt++) {
	  auto const& pos = trk.LocationAtPoint(pt);
	  trj.push_back(::geoalgo::Vector(pos[0], pos[1], pos[2]));
	}
	///debug
	FlashHypo Shiyan;
	std::vector<double> pe;
	//Shiyan.Set_Gap(0);
	//Shiyan.TrackStart(false);
	//Shiyan.TrackEnd(false);
	pe = Shiyan.FlashHypothesis(trj);
	_pe = pe;
	//std::cout<<pe[0]<<std::endl;
	_tree->Fill();
      }
    }
    return true;
  }

  bool Shiyan::finalize() {

    if(_fout){
      _fout->cd();
      if(_tree) _tree->Write();
    }

    
    return true;
  }

}
#endif
