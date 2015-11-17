#ifndef LARLITE_SHIYAN_CXX
#define LARLITE_SHIYAN_CXX

#include "Shiyan.h"
#include "FlashHypo.h"

namespace larlite {

  bool Shiyan::initialize() {

    if (_tree) {delete _tree;}
    _tree = new TTree("Paddletree", "Paddle Tree");
    _tree->Branch("PE","std::vector<double>",&_pe);
    
    _tree->Branch("run",&_run,"_run/I");
    _tree->Branch("subrun",&_subrun,"_subrun/I");
    _tree->Branch("event",&_event,"_event/I");
    _tree->Branch("trk_id",&_trk_id,"_trk_id/I");
    
    _length_xfiducial = larutil::Geometry::GetME()->DetHalfWidth();
    _length_yfiducial = larutil::Geometry::GetME()->DetHalfHeight();
    _length_zfiducial = larutil::Geometry::GetME()->DetLength();
    
    _vfiducial = ::geoalgo::AABox(0, -_length_yfiducial, 0, 2 * _length_xfiducial, _length_yfiducial,_length_zfiducial);
    _vmucs_top = ::geoalgo::AABox(-71.795, 393.941, 531.45, -23.795, 398.451, 579.45);
    _vmucs_bottom = ::geoalgo::AABox(-19.6948, 316.041, 533.25, 28.3052, 320.551, 581.25);

    _n_evt = 0;
    
    return true;
  }
  
  bool Shiyan::analyze(storage_manager* storage) {
    std::vector<double> pe;
    pe.resize(32,0.0);
    _pe.resize(32,0.0);
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
	::geoalgo::HalfLine trj_prj(trj[0], trj[0]-trj[trj.size()/2]);
	auto const&  intersection_trj_prj_top    = _geoAlgo.Intersection(_vmucs_top,trj_prj);
        auto const&  intersection_trj_prj_bottom = _geoAlgo.Intersection(_vmucs_bottom,trj_prj);
	if(intersection_trj_prj_top.size()>0 && intersection_trj_prj_bottom.size()>0){
	  ///debug
	  _n_evt++;
	  FlashHypo Shiyan;
	  //Shiyan.Set_Gap(0);
	  Shiyan.TrackStart(true);
	  Shiyan.TrackEnd(true);
	  std::cout<<Shiyan.FlashHypothesis(trj).size()<<std::endl;
	  std::cout<<trj.size()<<std::endl;
	  /*
	    for(int pmt_id=0; pmt_id<32;pmt_id++){
	    pe[pmt_id] += Shiyan.FlashHypothesis(trj)[pmt_id];
	    _pe[pmt_id] = pe[pmt_id];
	    }
	  */
	}
      }//For track with at least 2 hits
    }//Loop over each track in this event
    _run    = storage->get_data<event_track>("trackkalmanhit")->run();
    _subrun = storage->get_data<event_track>("trackkalmanhit")->subrun();
    _event  = storage->get_data<event_track>("trackkalmanhit")->event_id();
    //_trk_id = trk.ID();
    _tree->Fill();
    return true;
  }

  bool Shiyan::finalize() {

    if(_fout){
      _fout->cd();
      if(_tree) _tree->Write();
    }

    std::cout<<"hey yo:"<<_n_evt<<std::endl;
    return true;
  }

}
#endif
