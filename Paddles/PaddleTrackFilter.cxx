#ifndef LARLITE_PADDLETRACKFILTER_CXX
#define LARLITE_PADDLETRACKFILTER_CXX

#include "PaddleTrackFilter.h"
#include "DataFormat/track.h"

namespace larlite {

  bool PaddleTrackFilter::initialize() {
    
    if(_tree){delete _tree;}
    _tree = new TTree("paddletree","Paddle Tree");

    _length_xfiducial = larutil::Geometry::GetME()->DetHalfWidth();
    _length_yfiducial = larutil::Geometry::GetME()->DetHalfHeight();
    _length_zfiducial = larutil::Geometry::GetME()->DetLength();
    
    _vfiducial = ::geoalgo::AABox(0,-_length_yfiducial,0,2*_length_xfiducial,_length_yfiducial,_length_zfiducial);
    _vmucs_top = ::geoalgo::AABox(-71.795, 398.451, 579.45, -23.695, 398.351, 531.45);
    _vmucs_bottom= ::geoalgo::AABox(-19.6948, 320.05, 581.25, 28.3052, 320.15, 533.25);

    _n_intersections_FV = 0;
    _n_intersections_mucs_top = 0;
    _n_intersections_mucs_bottom = 0;
    
    return true;
  }
  
  bool PaddleTrackFilter::analyze(storage_manager* storage) {
    
    auto ev_reco = storage->get_data<event_track>("trackkalmanhit");
    if(!ev_reco){
      std::cout<<"........Couldn't find reco track data product in this event...... "<<std::endl;
    }

    for (size_t i=0; i < ev_reco->size(); i++){
      
      auto const& trk = ev_reco->at(i);
      
      std::cout << "this track has " << trk.NumberTrajectoryPoints() << " steps" << std::endl;
      
      if (trk.NumberTrajectoryPoints() > 1){

	// a geoalgo::Trajectory
	::geoalgo::Trajectory trj;
	
	for (size_t pt=0; pt < trk.NumberTrajectoryPoints(); pt++){

	  auto const& pos = trk.LocationAtPoint(pt);
	  
	  trj.push_back(::geoalgo::Vector(pos[0],pos[1],pos[2]));

	}// for all points in track
	
	std::cout << "the total track length is " << trj.Length() << std::endl;
	
	// make a HalfLine that starts at the beginning of the track and points
	// backwards in the track direction:
	
	::geoalgo::HalfLine trj_prj(trj[0], trj[0]-trj[trj.size()-1]);

	auto const& intersections_trj     = _geoAlgo.Intersection(_vfiducial,trj);
	auto const& intersections_trj_prj_bottom = _geoAlgo.Intersection(_vmucs_bottom,trj_prj);
	auto const& intersections_trj_prj_top = _geoAlgo.Intersection(_vmucs_bottom,trj_prj);
	
	if (intersections_trj.size() > 0){
	  _n_intersections_FV++;
	  std::cout << "this track intersects the Fiducial Volume!" << std::endl;}
	if (intersections_trj_prj_bottom.size() > 0){
	  _n_intersections_mucs_top++;
	  std::cout << "this track's projection backwards intersects the MuCS Bottom" << std::endl;}
	if (intersections_trj_prj_top.size() > 0){
	  _n_intersections_mucs_bottom++;
	  std::cout << "this track's projection backwards intersects the MuCS Top" << std::endl;}
	
	std::cout << std::endl;
	
	
      }// if there are at least 2 points in the track

    }// for all tracks
    
    return true;
  }

  bool PaddleTrackFilter::finalize() {

    if(_tree)
      _tree->Write();
    std::cout<<"track intersections with FV: "<<_n_intersections_FV<<std::endl;
    std::cout<<"track intersections with MuCS Top Detector: "<<_n_intersections_mucs_top<<std::endl;
    std::cout<<"track intersections with MuCS Bottome Detector: "<<_n_intersections_mucs_bottom<<std::endl;

    return true;
  }

}
#endif
