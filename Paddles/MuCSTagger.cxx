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

    _hit_upper_box = main_cfg.get<bool>("HitUpperBox");

    _hit_lower_box = main_cfg.get<bool>("HitLowerBox");

    if(!_hit_upper_box && !_hit_lower_box) {
      print(msg::kERROR,__FUNCTION__,"Both HitUpperBox and HitLowerBox cannot be false!");
      throw std::exception();
    }

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

    bool res = true;

    ::geoalgo::GeoAlgo alg;

    _mucs_dir.Start(start[0],start[1],start[2]);

    _mucs_dir.Dir(end[0]-start[0], end[1]-start[2], end[1]-start[2]);
    
    if(_hit_upper_box && alg.Intersection(_mucs_upper_box,_mucs_dir).empty()) res=false;

    if(res && _hit_lower_box && alg.Intersection(_mucs_lower_box,_mucs_dir).empty()) res=false;

    return res;

  }
  
  bool MuCSTagger::analyze(storage_manager* storage) {

    auto ev_track = storage->get_data<event_track>(_producer);
    if(!ev_track) {
      print(msg::kERROR,__FUNCTION__,"No matching data product found for specified track producer name!");
      throw DataFormatException("No data found");
    }

    ::geoalgo::HalfLine dir(0,0,0,1,1,1);
    
    auto ev_ctag = storage->get_data<event_cosmictag>("MuCSTagger");

    size_t num_tracks = 0;
    size_t num_tagged = 0;
    for(auto const& trk : *ev_track) {

      if(trk.NumberTrajectoryPoints()<2) continue;

      auto const& start = trk.LocationAtPoint(0);
      auto const& end   = trk.LocationAtPoint(trk.NumberTrajectoryPoints()-1);

      if(start[0] < _xmin || start[0] > _xmax) continue;
      if(end[0]   < _xmin || end[0]   > _xmax) continue;

      cosmictag ctag;
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

      if(length < _min_track_length) continue;

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

	if(ev_ctag) ev_ctag->emplace_back(ctag);
	
	++num_tagged;
      }
    }

    _hNumTracks->Fill(num_tracks);
    _hNumTagged->Fill(num_tagged);

    storage->set_id(storage->run_id(),
		    storage->subrun_id(),
		    storage->event_id());
    
    return true;
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
