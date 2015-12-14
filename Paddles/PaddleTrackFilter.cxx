#ifndef LARLITE_PADDLETRACKFILTER_CXX
#define LARLITE_PADDLETRACKFILTER_CXX

#include "PaddleTrackFilter.h"
#include "DataFormat/track.h"

namespace larlite {

bool PaddleTrackFilter::initialize() {

  if (_tree) {delete _tree;}
  _tree = new TTree("paddletree", "Paddle Tree");

  _length_xfiducial = larutil::Geometry::GetME()->DetHalfWidth();
  _length_yfiducial = larutil::Geometry::GetME()->DetHalfHeight();
  _length_zfiducial = larutil::Geometry::GetME()->DetLength();

  _vfiducial = ::geoalgo::AABox(0, -_length_yfiducial, 0, 2 * _length_xfiducial, _length_yfiducial, _length_zfiducial);
  //Postion of MuCS https://docs.google.com/spreadsheets/d/11XsZOU9kNe363-j1mTPTD-sKoNshvGotC5yMcyzXv18/edit#gid=1930333484
  _vmucs_top = ::geoalgo::AABox(-71.795, 393.941, 531.45, -23.795, 398.451, 579.45);
  _vmucs_bottom = ::geoalgo::AABox(-19.6948, 316.041, 533.25, 28.3052, 320.551, 581.25);

  _n_intersections_FV = 0;
  _n_intersections_mucs_top = 0;
  _n_intersections_mucs_bottom = 0;

  return true;
}

bool PaddleTrackFilter::analyze(storage_manager* storage) {

  _n_evt++;
  
  _trj.clear();
  _trj_con.clear();
  _trj_mucs.clear();
  _prj_lineseg.clear();
  
  auto ev_reco = storage->get_data<event_track>("trackkalmanhit");
  if (!ev_reco) {
    std::cout<<"........Couldn't find reco track data product in this event...... "<<std::endl;
  }

  for (size_t i = 0; i < ev_reco->size(); i++) {

    auto const& trk = ev_reco->at(i);
    //std::cout << "this track has " << trk.NumberTrajectoryPoints() << " steps" << std::endl;
    _n_tracks++;
    
    if (trk.NumberTrajectoryPoints() > 1) {

      // a geoalgo::Trajectory
      ::geoalgo::Trajectory trj;

      for (size_t pt = 0; pt < trk.NumberTrajectoryPoints(); pt++) {

        auto const& pos = trk.LocationAtPoint(pt);

        trj.push_back(::geoalgo::Vector(pos[0], pos[1], pos[2]));
      }// for all points in track
      
      //std::cout << "the total track length is " << trj.Length() << std::endl;
      
      // make a HalfLine that starts at the beginning of the track and points
      // backwards in the track direction:
      
      ::geoalgo::HalfLine trj_prj_start(trj[0], trk.DirectionAtPoint(0));
      ::geoalgo::HalfLine trj_prj_Negstart(trj[0], -trk.DirectionAtPoint(0));
      ::geoalgo::HalfLine trj_prj_end(trj.back(), trk.EndDirection());
      ::geoalgo::HalfLine trj_prj_Negend(trj.back(), -trk.EndDirection());

      //std::cout << trj[0] << std::endl;

      auto const& intersections_trj_start     = _geoAlgo.Intersection(_vfiducial, trj);
      auto const& intersections_trj_prj_bottom_start = _geoAlgo.Intersection(_vmucs_bottom, trj_prj_start);
      auto const& intersections_trj_prj_top_start = _geoAlgo.Intersection(_vmucs_top, trj_prj_start);
      
      auto const& intersections_trj_Negstart     = _geoAlgo.Intersection(_vfiducial, trj);
      auto const& intersections_trj_prj_bottom_Negstart = _geoAlgo.Intersection(_vmucs_bottom, trj_prj_Negstart);
      auto const& intersections_trj_prj_top_Negstart = _geoAlgo.Intersection(_vmucs_top, trj_prj_Negstart);

      auto const& intersections_trj_end     = _geoAlgo.Intersection(_vfiducial, trj);
      auto const& intersections_trj_prj_bottom_end = _geoAlgo.Intersection(_vmucs_bottom, trj_prj_end);
      auto const& intersections_trj_prj_top_end = _geoAlgo.Intersection(_vmucs_top, trj_prj_end);
      
      auto const& intersections_trj_Negend     = _geoAlgo.Intersection(_vfiducial, trj);
      auto const& intersections_trj_prj_bottom_Negend = _geoAlgo.Intersection(_vmucs_bottom, trj_prj_Negend);
      auto const& intersections_trj_prj_top_Negend = _geoAlgo.Intersection(_vmucs_top, trj_prj_Negend);

      
      if (intersections_trj_start.size() > 0 || intersections_trj_end.size() > 0 || intersections_trj_Negstart.size() > 0 ||
           intersections_trj_Negend.size() > 0) {
        //_n_intersections_FV++;
      }
      if (intersections_trj_prj_bottom_start.size() > 0 || intersections_trj_prj_bottom_end.size() > 0 ||
          intersections_trj_prj_bottom_Negstart.size() > 0 || intersections_trj_prj_bottom_Negend.size() > 0) {
        _n_intersections_mucs_top++;
      }
      
      if (intersections_trj_prj_top_start.size() > 0 || intersections_trj_prj_top_end.size() > 0 ||
          intersections_trj_prj_top_Negstart.size() > 0 || intersections_trj_prj_top_Negend.size() > 0) {
        _n_intersections_mucs_bottom++;
	}
      //store tracks proj by lineseg
      ::geoalgo::LineSegment prj_lineseg((trj[0]+(trj[0]-trj[trj.size()/2-1])*7), trj[trj.size()/2-1]);
       _prj_lineseg.push_back(prj_lineseg);
      //store tracks in an event contained by TPCFV
      if(_vfiducial.Contain(trj.at(0))==1){_trj_con.push_back(trj);}
      //store mucs tracks in an event
      //::geoalgo::HalfLine trj_prj(trj[0], trj[0]-trj[trj.size()/2]);
      ::geoalgo::HalfLine trj_prj(trj[0], trj[0]-trj[1]);
      _trj_prj.push_back(trj_prj);
      auto const& intersections_trj_prj_top = _geoAlgo.Intersection(_vmucs_top, trj_prj);
      auto const& intersections_trj_prj_bottom = _geoAlgo.Intersection(_vmucs_bottom, trj_prj);
      
      if(intersections_trj_prj_top.size()>0 || intersections_trj_prj_bottom.size()>0){
	_trj_mucs.push_back(trj);
	_n_intersections_FV++;

	_run    = storage->get_data<event_track>("trackkalmanhit")->run();
	_subrun = storage->get_data<event_track>("trackkalmanhit")->subrun();
	_event  = storage->get_data<event_track>("trackkalmanhit")->event_id();
	_trk_id = trk.ID();
	
      }
      //store all  tracks in an event
      _trj.push_back(trj);  
    }// if there are at least 2 points in the track
    
  }// for all tracks

  return true;
}


bool PaddleTrackFilter::finalize() {

  if (_tree)
    _tree->Write();
  std::cout << "track intersections with FV: " << _n_intersections_FV << std::endl;
  std::cout << "track intersections with MuCS Bottome Detector: " << _n_intersections_mucs_bottom << std::endl;
  std::cout << "track intersections with MuCS Top Detector: " << _n_intersections_mucs_top << std::endl;
  return true;
}

}
#endif
