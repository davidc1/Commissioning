#ifndef LARLITE_PADDLETRACKOPFLASH_CXX
#define LARLITE_PADDLETRACKOPFLASH_CXX

#include "PaddleTrackOpflash.h"
#include <algorithm>
#include <functional>

namespace larlite {
  
  bool PaddleTrackOpflash::initialize() {
    
    if (_tree) {delete _tree;}
    _tree = new TTree("Paddle Tree", "Paddle Tree");
    _tree->Branch("t_opflash","std::vector<double>",&_t_opflash);
    
    _length_xfiducial = larutil::Geometry::GetME()->DetHalfWidth();
    _length_yfiducial = larutil::Geometry::GetME()->DetHalfHeight();
    _length_zfiducial = larutil::Geometry::GetME()->DetLength();

    _vfiducial = ::geoalgo::AABox(0, -_length_yfiducial, 0, 2 * _length_xfiducial, _length_yfiducial,_length_zfiducial);
    
    _vmucs_top = ::geoalgo::AABox(-71.795, 393.941, 531.45, -23.795, 398.451, 579.45);
    _vmucs_bottom = ::geoalgo::AABox(-19.6948, 316.041, 533.25, 28.3052, 320.551, 581.25);

    _n_evt = 0;
    
    return true;
  }
  
  bool PaddleTrackOpflash::analyze(storage_manager* storage) {
    
    _t_opflash.clear();
    
    auto ev_reco = storage->get_data<event_track>("trackkalmanhit");
    if (!ev_reco) {
      std::cout<<"........Couldn't find reco track data product in this event...... "<<std::endl;
    }

    auto ev_ophit = storage->get_data<event_ophit>("opFlash");
    if (!ev_ophit) {
      std::cout<<"........Couldn't find ophit data product in this event...... "<<std::endl;
    }
    
    auto ev_opflash = storage->get_data<event_opflash>("opFlash");
    if (!ev_opflash) {
      std::cout<<"........Couldn't find opflash data product in this event...... "<<std::endl;
    }
    
    for(size_t i = 0; i <ev_reco->size(); i++ ){
      
      auto const& trk = ev_reco->at(i);

      if (trk.NumberTrajectoryPoints() > 1) {

        //a geoalgo::Trajectory            
        ::geoalgo::Trajectory trj;

        for (size_t pt = 0; pt < trk.NumberTrajectoryPoints(); pt++) {

          auto const& pos = trk.LocationAtPoint(pt);

          trj.push_back(::geoalgo::Vector(pos[0], pos[1], pos[2]));
        }// for all points in track

	::geoalgo::HalfLine trj_prj(trj[0], trj[0]-trj[trj.size()/2]);
        ::geoalgo::HalfLine trj_prj_neg(trj[trj.size()-1], trj[trj.size()/2]-trj[trj.size()-1]);

        auto const&  intersection_trj_prj_top    = _geoAlgo.Intersection(_vmucs_top,trj_prj);
        auto const&  intersection_trj_prj_bottom = _geoAlgo.Intersection(_vmucs_bottom,trj_prj);
	
        auto const&  intersection_trj_prj_top_neg    = _geoAlgo.Intersection(_vmucs_top,trj_prj_neg);
        auto const&  intersection_trj_prj_bottom_neg = _geoAlgo.Intersection(_vmucs_bottom,trj_prj_neg);
	
	if(intersection_trj_prj_top.size()>0 && intersection_trj_prj_bottom.size()>0){
          _n_evt++;
	  //std::cout<<intersection_trj_prj_top.size()<<",";
	  _MuCS_ints_x_top = intersection_trj_prj_top.at(0).at(0);
          _MuCS_ints_z_top = intersection_trj_prj_top.at(0).at(2);
          _MuCS_ints_x_bottom =  intersection_trj_prj_bottom.at(0).at(0);
          _MuCS_ints_z_bottom =  intersection_trj_prj_bottom.at(0).at(2);
	  
	  //
	  {
	    std::cout<<ev_opflash->size()<<std::endl;
	    for(size_t j = 0; j < ev_opflash->size(); j++){
	      auto const& opflash = ev_opflash->at(j);
	      _t_opflash.push_back(opflash.Time());
	      
	    }
	    std::cout<<_t_opflash.size()<<std::endl;
	    std::cout<<_t_opflash.at(0)<<std::endl;
	  }//for all opflashes
	  
	  
	  _tree->Fill();
        }
        //correct for wrong recon direction
        if(intersection_trj_prj_top_neg.size()>0 && intersection_trj_prj_bottom_neg.size()>0){
          _n_evt++;
	  //std::cout<<intersection_trj_prj_top_neg.size();
          _MuCS_ints_x_top = intersection_trj_prj_top_neg.at(0).at(0);
          _MuCS_ints_z_top = intersection_trj_prj_top_neg.at(0).at(2);
          _MuCS_ints_x_bottom =  intersection_trj_prj_bottom_neg.at(0).at(0);
          _MuCS_ints_z_bottom =  intersection_trj_prj_bottom_neg.at(0).at(2);
	  _tree->Fill();
        }

	
      }
    }///for all tracks
   
    return true;
  }
  
  bool PaddleTrackOpflash::finalize() {
    if(_fout){
      _fout->cd();
      if(_tree) _tree->Write();
      std::cout<<"\nThrough-going muon found in "<<_n_evt<<" events!!!???"<<std::endl;
    }
    
    
    return true;
  }

}
#endif
