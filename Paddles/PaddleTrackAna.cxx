#ifndef LARLITE_PADDLETRACKANA_CXX
#define LARLITE_PADDLETRACKANA_CXX

#include "PaddleTrackAna.h"
#include "DataFormat/track.h"

namespace larlite {

  bool PaddleTrackAna::initialize() {
    if(!_x_z_ints_top)
      _x_z_ints_top = new TH2F("_x_z_ints_top","Intersection(x,z)_Top_Dectecor",12,-71.795,-23.795,12,531.45,579.45);
    if(!_x_z_ints_bottom)
      _x_z_ints_bottom = new TH2F("_x_z_ints_bottom","Intersection(x,z)_Bottom_Dectecor",12,-19.6948,28.3052,12,533.25,581.25);
    
    
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

    //Postion of MuCS https://docs.google.com/spreadsheets/d/11XsZOU9kNe363-j1mTPTD-sKoNshvGotC5yMcyzXv18/edit#gid=1930333484
    //_vmucs_top = ::geoalgo::AABox(-71.795, 393.941, 531.45, -23.795, 398.451, 579.45);
    _vmucs_bottom = ::geoalgo::AABox(-19.6948, 316.041, 533.25, 28.3052, 320.551, 581.25);
    
    //Double the area the MuCS, correction for not-well construted track
    //_vmucs_top = ::geoalgo::AABox(-95.795, 393.941, 507.45, 0.205, 398.451, 603.45);
    //_vmucs_bottom = ::geoalgo::AABox(-43.6948, 316.041, 509.25, 52.3052, 320.551, 605.25);
    
    //Plane y=316.041 and 398.451    
    _vmucs_top = ::geoalgo::AABox(-256, 393.941, 0, 256, 393.551, 1037);
    //_vmucs_bottom = ::geoalgo::AABox(-256, 319.551, 0, 256, 320.551, 1037);
    
    
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
	::geoalgo::HalfLine trj_prj_neg(trj[trj.size()-1], trj[trj.size()/2]-trj[trj.size()-1]);
	
	auto const&  intersection_trj_prj_top    = _geoAlgo.Intersection(_vmucs_top,trj_prj);
	auto const&  intersection_trj_prj_bottom = _geoAlgo.Intersection(_vmucs_bottom,trj_prj);
	
	auto const&  intersection_trj_prj_top_neg    = _geoAlgo.Intersection(_vmucs_top,trj_prj_neg);
        auto const&  intersection_trj_prj_bottom_neg = _geoAlgo.Intersection(_vmucs_bottom,trj_prj_neg);
	

	if(intersection_trj_prj_top.size()>0 && intersection_trj_prj_bottom.size()>0){
	  _n_evt++;
	  //if(intersection_trj_prj_top.size()>0 && intersection_trj_prj_bottom.size()>0){
	  
	  //std::cout<<intersection_trj_prj_top.size();
	  
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

	  // get ID information
	  _run = storage->get_data<event_track>("trackkalmanhit")->run();
	  _subrun = storage->get_data<event_track>("trackkalmanhit")->subrun();
	  _event = storage->get_data<event_track>("trackkalmanhit")->event_id();
	  _trk_id = trk.ID();
	  
	  //std::cout<<"Run: "<<_run<<" , Subrun: "<<_subrun<<" , Event: "<<_event<<" , Track_ID: "<<_trk_id<<std::endl;
	  
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
          _x_z_ints_top->Fill(intersection_trj_prj_top_neg.at(0).at(0),intersection_trj_prj_top_neg.at(0).at(2));
          _x_z_ints_bottom->Fill(intersection_trj_prj_bottom_neg.at(0).at(0),intersection_trj_prj_bottom_neg.at(0).at(2));
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
      std::cout<<"\nThrough-going muon found in "<<_n_evt<<" events!!!???"<<std::endl;
    }

    if(_save_histos){
      _x_z_ints_top->Write();
      _x_z_ints_bottom->Write();
    }
    return true;
  }

}
#endif
