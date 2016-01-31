#ifndef LARLITE_MUCSTAGGER_CXX
#define LARLITE_MUCSTAGGER_CXX

#include "MuCSTagger.h"
#include "FhiclLite/ConfigManager.h"
#include "DataFormat/track.h"
#include "DataFormat/cosmictag.h"
namespace larlite {

  MuCSTagger::MuCSTagger()
    : ana_base()
    , _mucs_upper_box(0,0,0,1,1,1)
    , _mucs_lower_box(0,0,0,1,1,1)
    , _mucs_dir(0,0,0,1,1,1)
  { _name="MuCSTagger"; _fout=0; _configured=false;}

  void MuCSTagger::configure(const std::string config_file) {

    ::fcllite::ConfigManager cfg_mgr(_name);

    cfg_mgr.AddCfgFile(config_file);

    auto const& main_cfg = cfg_mgr.Config().get_pset(_name);

    _store_match = main_cfg.get<bool>("StoreMatch");

    std::vector<double> upper_box = main_cfg.get<std::vector<double> >("UpperBox");
    if(upper_box.size()!=6) {
      print(msg::kERROR,__FUNCTION__,"UpperBox needs length 6 floating point...");
      throw std::exception();
    }
    std::vector<double> lower_box = main_cfg.get<std::vector<double> >("LowerBox");
    if(lower_box.size()!=6) {
      print(msg::kERROR,__FUNCTION__,"LowerBox needs length 6 floating point...");
      throw std::exception();
    }

    _mucs_upper_box.Min(upper_box[0],upper_box[1],upper_box[2]);
    _mucs_upper_box.Max(upper_box[3],upper_box[4],upper_box[5]);

    _mucs_lower_box.Min(lower_box[0],lower_box[1],lower_box[2]);
    _mucs_lower_box.Max(lower_box[3],lower_box[4],lower_box[5]);

    _xmin = main_cfg.get<double>("XMin",-50.);
    _xmax = main_cfg.get<double>("XMax",300.);

    _producer = main_cfg.get<std::string>("Producer","trackkalmanhit");

    _hit_both_box = main_cfg.get<bool>("HitBoth");

    _min_track_length = main_cfg.get<double>("MinTrackLength",50.);

    _scan_length = main_cfg.get<double>("ScanLength",0.);

    _segment_length = main_cfg.get<double>("SegmentLength",10.);

    if(_segment_length>_min_track_length) {
      print(msg::kERROR,__FUNCTION__,"Segment length cannot be larger than MinTracklength!");
      throw std::exception();
    }

    _allow_flip_direction = main_cfg.get<bool>("AllowFlip",true);

    _configured = true;
  }

  bool MuCSTagger::initialize() {

    if(!_configured) {
      print(msg::kERROR,__FUNCTION__,"Must be configured before running!");
      return false;
    }
    
    _hNumTracks = new TH1D("hNumTracks","# Tracks above min length; # Tracks; Events",
			   20,0,100);

    _hNumTagged = new TH1D("hNumTagged","# Tracks tagged; # tagged; Events",
			   20,-0.5,19.5);

    return true;
  }

  bool MuCSTagger::Intersect(const TVector3& start, const TVector3& end) {

    ::geoalgo::GeoAlgo alg;

    _mucs_dir.Start(start[0],start[1],start[2]);

    _mucs_dir.Dir(end[0]-start[0], end[1]-start[1], end[2]-start[2]);

    bool res = alg.Intersection(_mucs_upper_box,_mucs_dir).empty();

    if(!res && _hit_both_box) return false;
    
    if(res && !_hit_both_box) return true;

    return !(alg.Intersection(_mucs_lower_box,_mucs_dir).empty());

  }

  bool MuCSTagger::IntersectDumb(const TVector3& start, const TVector3& end) {

    TVector3 dir = end - start;
    double mag = dir.Mag();
    dir[0] /= mag;
    dir[1] /= mag;
    dir[2] /= mag;
    if(dir[1]>0) return false;
    std::cout<<"start "<<start[0]<<" "<<start[1]<<" "<<start[2]<<std::endl;
    std::cout<<"end   "<<end[0]<<" "<<end[1]<<" "<<end[2]<<std::endl;
    std::cout<<"dir   "<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<std::endl;

    auto const& upper_min = _mucs_upper_box.Min();
    auto const& upper_max = _mucs_upper_box.Max();
    auto const& lower_min = _mucs_lower_box.Min();
    auto const& lower_max = _mucs_lower_box.Max();

    double upper_x = start[0] + dir[0] * upper_min[1];
    double upper_z = start[2] + dir[2] * upper_min[1];

    double lower_x = start[0] + dir[0] * lower_min[1];
    double lower_z = start[2] + dir[2] * lower_min[1];

    std::cout<<"upper xs: "<<upper_x<<" "<<upper_z<<std::endl;
    std::cout<<"lower xs: "<<lower_x<<" "<<lower_z<<std::endl;

    bool upper_hit = ( upper_x > upper_min[0] && upper_x < upper_max[0] &&
		       upper_z > upper_min[2] && upper_z < upper_max[2] );

    bool lower_hit = ( lower_x > lower_min[0] && lower_x < lower_max[0] &&
		       lower_z > lower_min[2] && lower_z < lower_max[2] );

    std::cout << "upper: " << (upper_hit ? "yes" : "no") << " ... "
	      << "lower: " << (upper_hit ? "yes" : "no") << std::endl << std::endl;
    
    if(_hit_both_box) return (upper_hit && lower_hit);
    return (upper_hit || lower_hit);

  }
  
  bool MuCSTagger::analyze(storage_manager* storage) {

    _matched_trj_v.clear();
    _matched_dir_v.clear();

    auto ev_track = storage->get_data<event_track>(_producer);
    if(!ev_track) {
      print(msg::kERROR,__FUNCTION__,"No matching data product found for specified track producer name!");
      throw DataFormatException("No data found");
    }

    ::geoalgo::HalfLine dir(0,0,0,1,1,1);

    auto ev_ctag = storage->get_data<event_cosmictag>("MuCSTagger");
    if(ev_ctag) ev_ctag->resize(ev_track->size());
    //if(ev_ctag) ev_ctag->reserve(ev_track->size());

    cosmictag ctag;
    size_t num_tracks = 0;
    size_t num_tagged = 0;
    bool tagged = false;
    for(size_t track_index=0; track_index<ev_track->size(); ++track_index) {

      ctag.fCosmicScore = -1;
      auto const& trk = (*ev_track)[track_index];

      if(trk.NumberTrajectoryPoints()<2) {
	if(ev_ctag) (*ev_ctag)[track_index].fCosmicScore = -1;
	continue;
      }
      
      auto const& start = trk.LocationAtPoint(0);
      auto const& end   = trk.LocationAtPoint(trk.NumberTrajectoryPoints()-1);

      if( (start[0] < _xmin || start[0] > _xmax) || (end[0]   < _xmin || end[0]   > _xmax) ) {
	//if(ev_ctag) ev_ctag->push_back(ctag);
	if(ev_ctag) (*ev_ctag)[track_index].fCosmicScore = -1;
	continue;
      }

      ctag.fCosmicScore = -1;
	
      double length = 0;
      size_t ref_index = 0;
      for(size_t i=0; i<(trk.NumberTrajectoryPoints()-1) && ctag.fCosmicScore<0; ++i) {
	auto const& pt1 = trk.LocationAtPoint(i);
	auto const& pt2 = trk.LocationAtPoint(i+1);
	length += (pt2 - pt1).Mag();
	if(length >= _segment_length && length < _segment_length + _scan_length) {
	  auto const& ref = trk.LocationAtPoint(ref_index);
	  if(Intersect(pt2,ref)) ctag.fCosmicScore = 1.;
	  ++ref_index;
	}
	if(length > _min_track_length) break;
      }

      if(length < _min_track_length) {
	//if(ev_ctag) ev_ctag->push_back(ctag);
	if(ev_ctag) (*ev_ctag)[track_index].fCosmicScore = -1;
	continue;
      }
      
      ++num_tracks;
      
      if(ctag.fCosmicScore<0. && _allow_flip_direction) {
	length = 0;
	ref_index = trk.NumberTrajectoryPoints()-1;
	for(int i=trk.NumberTrajectoryPoints()-1; i>0 && ctag.fCosmicScore<0.; --i) {
	  auto const& pt1 = trk.LocationAtPoint(i);
	  auto const& pt2 = trk.LocationAtPoint(i-1);
	  length += (pt2 - pt1).Mag();
	  if(length >= _segment_length && length < _segment_length + _scan_length) {
	    auto const& ref = trk.LocationAtPoint(ref_index);
	    if(Intersect(pt2,ref)) ctag.fCosmicScore = 0.5;
	    --ref_index;
	  }
	  if(length > (_scan_length + _segment_length)) break;
	}
      }

      if(ctag.fCosmicScore>0) {
	tagged = true;
	if(_store_match) {
	  _matched_dir_v.push_back(_mucs_dir);
	  ::geoalgo::Trajectory trj;
	  trj.resize(trk.NumberTrajectoryPoints());
	  for(size_t i=0; i<trk.NumberTrajectoryPoints(); ++i) {

	    auto const& pt = trk.LocationAtPoint(i);

	    trj.emplace_back(::geoalgo::Vector(pt));

	  }
	  _matched_trj_v.emplace_back(trj);
	}	  
	++num_tagged;
      }
      //if(ctag.fCosmicScore>0) std::cout<<"\033[93m"<<track_index<<"\033[00m"<<std::endl;
      //if(ev_ctag) ev_ctag->push_back(ctag);
      if(ev_ctag) (*ev_ctag)[track_index].fCosmicScore = ctag.fCosmicScore;
    }

    _hNumTracks->Fill(num_tracks);
    _hNumTagged->Fill(num_tagged);

    storage->set_id(storage->run_id(),
		    storage->subrun_id(),
		    storage->event_id());

    if(ev_ctag && ev_ctag->size()!=ev_track->size()) {
      print(msg::kERROR,__FUNCTION__,"Mismatch in ctag & track length!");
      throw std::exception();
    }

    return tagged;
  }

  bool MuCSTagger::finalize() {

    if(_fout) {
      _fout->cd();
      _hNumTracks->Write();
      _hNumTagged->Write();
    }
  
    return true;
  }

}
#endif
