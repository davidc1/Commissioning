#ifndef LARLITE_PADDLETRACKANA_CXX
#define LARLITE_PADDLETRACKANA_CXX

#include "PaddleTrackAna.h"
#include "DataFormat/track.h"

namespace larlite {

  bool PaddleTrackAna::initialize() {
    if(!_x_z_ints_top)
      _x_z_ints_top = new TH2F("_x_z_ints_top","Intersection(x,c)_Top_Dectecor",49,-72,-23,49,531,580);
    if(!_x_z_ints_bottom)
      _x_z_ints_bottom = new TH2F("_x_z_ints_top","Intersection(x,c)_Bottome_Dectecor",49,-20,-29,49,533,582);
    
    
    if (_tree) {delete _tree;}
    _tree = new TTree("Paddletree", "Paddle Tree");
    
    _tree->Branch("x_top",&_MuCS_ints_x_top,"_MuCS_ints_x_top/D");
    _tree->Branch("z_top",&_MuCS_ints_z_top,"_MuCS_ints_z_top/D");
    _tree->Branch("x_bottom",&_MuCS_ints_x_bottom,"_MuCS_ints_x_bottom/D");
    _tree->Branch("z_bottom",&_MuCS_ints_z_bottom,"_MuCS_ints_z_bottom/D");
    

    _length_xfiducial = larutil::Geometry::GetME()->DetHalfWidth();
    _length_yfiducial = larutil::Geometry::GetME()->DetHalfHeight();
    _length_zfiducial = larutil::Geometry::GetME()->DetLength();

    _vfiducial = ::geoalgo::AABox(0, -_length_yfiducial, 0, 2 * _length_xfiducial, _length_yfiducial, _length_zfiducial);
    _vmucs_top = ::geoalgo::AABox(-71.795, 398.451, 531.45, -23.695, 398.551, 579.45);
    _vmucs_bottom = ::geoalgo::AABox(-19.6948, 320.551, 533.25, 28.3052, 320.552, 581.25);

    _n_intersections_FV = 0;
    _n_intersections_mucs_top = 0;
    _n_intersections_mucs_bottom = 0;
    _n_evt =0;
    
    return true;
  }
  
  bool PaddleTrackAna::analyze(storage_manager* storage) {
    
    _trj.clear();
    _trj_con.clear();
    _trj_mucs.clear();
    _prj_lineseg.clear();
    
    auto ev_reco = storage->get_data<event_track>("trackkalmanhit");
    if (!ev_reco) {
      std::cout<<"........Couldn't find reco track data product in this event...... "<<std::endl;    
    }
    for(size_t i=0; i<ev_reco->size();i++){
      
      auto const& trk = ev_reco->at(i);
      
      if (trk.NumberTrajectoryPoints() > 1) {
	
	//a geoalgo::Trajectory
	::geoalgo::Trajectory trj;
	
	for (size_t pt = 0; pt < trk.NumberTrajectoryPoints(); pt++) {
	  
	  auto const& pos = trk.LocationAtPoint(pt);
	  
	  trj.push_back(::geoalgo::Vector(pos[0], pos[1], pos[2]));
	}// for all points in track 
      
	::geoalgo::HalfLine trj_prj(trj[0], trj[0]-trj[trj.size()/2]);
	auto const&  intersection_trj_prj_top    = _geoAlgo.Intersection(_vmucs_top,trj_prj);
	auto const&  intersection_trj_prj_bottom = _geoAlgo.Intersection(_vmucs_bottom,trj_prj);
	
	

	if(intersection_trj_prj_top.size()>0 && intersection_trj_prj_bottom.size()>0){
	  _n_evt++;
	  //if(intersection_trj_prj_top.size()>0 && intersection_trj_prj_bottom.size()>0){
	  std::cout<<intersection_trj_prj_top.size();
	  //std::cout<<intersection_trj_prj_top.at(0)<<std::endl;
	  //std::cout<<intersection_trj_prj_top.at(0).at(0)<<std::endl;
	  //std::cout<<intersection_trj_prj_top.at(0).at(1)<<std::endl; 
	  _MuCS_ints_x_top = intersection_trj_prj_top.at(0).at(0);
          _MuCS_ints_z_top = intersection_trj_prj_top.at(0).at(2);	  
	  _MuCS_ints_x_bottom =  intersection_trj_prj_bottom.at(0).at(0);
	  _MuCS_ints_z_bottom =  intersection_trj_prj_bottom.at(0).at(2);
	  //std::cout<<_MuCS_ints_x_top<<std::endl;
	  //std::cout<<_MuCS_ints_y_top<<std::endl;	  
	  _x_z_ints_top->Fill(intersection_trj_prj_top.at(0).at(0),intersection_trj_prj_top.at(0).at(2));
	  _x_z_ints_bottom->Fill(intersection_trj_prj_bottom.at(0).at(0),intersection_trj_prj_bottom.at(0).at(2));
	  _tree->Fill();
	}
	
      }
    }//for all tracks
    
    return true;
  }

  bool PaddleTrackAna::finalize() {
    if(_fout){
      _fout->cd();
      
      if(_tree) _tree->Write();
      std::cout<<"Through-going muon found in "<<_n_evt<<" events!!!???"<<std::endl;
    }

    if(_save_histos){
      _x_z_ints_top->Write();
      _x_z_ints_bottom->Write();
    }
    return true;
  }

}
#endif
