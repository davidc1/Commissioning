#ifndef LARLITE_MUCSTAGGER_CXX
#define LARLITE_MUCSTAGGER_CXX

#include "MuCSTagger.h"
#include "FhiclLite/ConfigManager.h"
#include "DataFormat/track.h"
#include "DataFormat/cosmictag.h"
namespace larlite {

  MuCSTagger* MuCSTagger::_me = nullptr;
  
  MuCSTagger::MuCSTagger()
    : ana_base()
    , _mucs_upper_box(0,0,0,1,1,1)
    , _mucs_lower_box(0,0,0,1,1,1)
    , _tpc_av(0,0,0,1,1,1)
    , _mucs_dir(0,0,0,1,1,1)
    , _ctag_score(-1)
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

    std::vector<double> shift(6,0.);
    shift = main_cfg.get<std::vector<double> >("UpperShift",shift);
    if(shift.size()!=6) {
      print(msg::kERROR,__FUNCTION__,"Shift needs length 6 floating point...");
      throw std::exception();
    }
    for(size_t i=0;i<upper_box.size();++i) upper_box[i] += shift[i];

    shift = main_cfg.get<std::vector<double> >("LowerShift",shift);
    if(shift.size()!=6) {
      print(msg::kERROR,__FUNCTION__,"Shift needs length 6 floating point...");
      throw std::exception();
    }
    for(size_t i=0;i<lower_box.size();++i) lower_box[i] += shift[i];

    std::vector<double> tpc_box = main_cfg.get<std::vector<double> >("TPC_Active_Volume");
    if(tpc_box.size()!=6) {
      print(msg::kERROR,__FUNCTION__,"TPC AV needs length 6 floating point...");
      throw std::exception();
    }
    _mucs_upper_box.Min(upper_box[0],upper_box[1],upper_box[2]);
    _mucs_upper_box.Max(upper_box[3],upper_box[4],upper_box[5]);

    _mucs_lower_box.Min(lower_box[0],lower_box[1],lower_box[2]);
    _mucs_lower_box.Max(lower_box[3],lower_box[4],lower_box[5]);

    _tpc_av.Min(tpc_box[0],tpc_box[1],tpc_box[2]);
    _tpc_av.Max(tpc_box[3],tpc_box[4],tpc_box[5]);
    
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
			   20,-0.5,19.5);

    _hNumTagged = new TH1D("hNumTagged","# Tracks tagged; # tagged; Events",
			   10,-0.5,9.5);

    return true;
  }

  bool MuCSTagger::Intersect(const ::geoalgo::HalfLine& mucs_dir) {

    ::geoalgo::GeoAlgo alg;

    auto const upper_pt_v = alg.Intersection(_mucs_upper_box,mucs_dir);
    
    bool res = !(upper_pt_v.empty());

    // no intersection w/ upper box, both box hit required => return false;
    if(!res && _hit_both_box) return res;

    // intersection w/ upper box, only one box need to be hit => return true
    if(res && !_hit_both_box){
      _upper_pt.push_back(upper_pt_v.front());

      return res;
    }

    auto const lower_pt_v = alg.Intersection(_mucs_lower_box,mucs_dir);

    res = !(lower_pt_v.empty());

    if(res) _lower_pt.push_back(lower_pt_v.front());

    return res;
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

    _rui_v.clear();
    _matched_trj_v.clear();
    _matched_dir_v.clear();
    _upper_pt.clear();
    _lower_pt.clear();
    _temps1.clear();
    _temps2.clear();
    _ctag_score = -1;
    auto ev_track = storage->get_data<event_track>(_producer);
    if(!ev_track) {
      print(msg::kERROR,__FUNCTION__,"No matching data product found for specified track producer name!");
      throw DataFormatException("No data found");
    }

    _rui_v.resize(ev_track->size());
    
    ::geoalgo::HalfLine mucs_dir(0,0,0,1,1,1);

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
      ::geoalgo::HalfLine temp1(trk.LocationAtPoint(0),trk.LocationAtPoint(0)-trk.LocationAtPoint(1));
      ::geoalgo::HalfLine temp2(trk.LocationAtPoint(trk.NumberTrajectoryPoints()-1),trk.LocationAtPoint(trk.NumberTrajectoryPoints()-1)-trk.LocationAtPoint(trk.NumberTrajectoryPoints()-2));
      ::geoalgo::Trajectory trj_start;
      ::geoalgo::Trajectory trj_end;
      
      double length = 0;
      //size_t ref_index = 1;
      size_t ref_index = 0;
      for(size_t i=0; i<(trk.NumberTrajectoryPoints()-1) && length<_min_track_length; ++i) {
	//auto const& ptS = trk.LocationAtPoint(0);
	auto const& pt1 = trk.LocationAtPoint(i);
	auto const& pt2 = trk.LocationAtPoint(i+1);

	length += (pt2 - pt1).Mag();
	if(ctag.fCosmicScore<0 && length >= _segment_length && length < _segment_length + _scan_length) {
	  auto const& ref = trk.LocationAtPoint(ref_index);

	  //mucs_dir.Start(ref); mucs_dir.Dir(ptS-ref);
	  mucs_dir.Start(pt2); mucs_dir.Dir(ref-pt2);
	
	  if(Intersect(mucs_dir)){
	    ctag.fCosmicScore = 1.;
	    auto const _pt_av_up = _alg.Intersection(_tpc_av,mucs_dir);
	    trj_start.emplace_back(pt2);
	    if(!_pt_av_up.empty())_pt_up = _pt_av_up[0];
	    mucs_dir.Dir(-ref+pt2);
	    auto const _pt_av_up1 = _alg.Intersection(_tpc_av,mucs_dir);
	    if(!_pt_av_up1.empty())_pt_low = _pt_av_up1[0];
	    
	  }

	  ++ref_index;
	}
      }
      

      if(length < _min_track_length) {
	if(ev_ctag) (*ev_ctag)[track_index].fCosmicScore = -1;
	continue;
      }
      
      ++num_tracks;
      
      if(ctag.fCosmicScore<0. && _allow_flip_direction) {
	
	length = 0;
	//ref_index = trk.NumberTrajectoryPoints()-2;
	ref_index = trk.NumberTrajectoryPoints()-1;
	for(int i=trk.NumberTrajectoryPoints()-1; i>0 && length<_min_track_length; --i) {
	  //auto const& ptE = trk.LocationAtPoint(trk.NumberTrajectoryPoints()-1);
	  auto const& pt1 = trk.LocationAtPoint(i);
	  auto const& pt2 = trk.LocationAtPoint(i-1);

	  length += (pt2 - pt1).Mag();
	  if(ctag.fCosmicScore<0 && length >= _segment_length && length < _segment_length + _scan_length) {
	    auto const& ref = trk.LocationAtPoint(ref_index);

	    //mucs_dir.Start(ref); mucs_dir.Dir(ptE-ref);
	    mucs_dir.Start(pt2); mucs_dir.Dir(ref-pt2);

	    if(Intersect(mucs_dir)){
	      
	      auto const _pt_av_low = _alg.Intersection(_tpc_av,mucs_dir);
	      if(!_pt_av_low.empty())_pt_up = _pt_av_low[0];
	      mucs_dir.Dir(-ref+pt2);
	      auto const _pt_av_low1 = _alg.Intersection(_tpc_av,mucs_dir);
	      if(!_pt_av_low1.empty())_pt_low = _pt_av_low1[0];
	      
	    }
	    
	    if(Intersect(mucs_dir)) ctag.fCosmicScore = 0.5;

	    --ref_index;
	  }
	}
      }

      if(ctag.fCosmicScore>0) {

	tagged = true;
	if(_store_match) {
	  _matched_dir_v.push_back(_mucs_dir);
	  ::geoalgo::Trajectory trj;
	  trj.resize(trk.NumberTrajectoryPoints());
	  for(size_t i=0; i<trk.NumberTrajectoryPoints(); ++i) {

	    //auto const& pt = trk.LocationAtPoint(i);
	    //auto const& con = _tpc_av.Contain(pt);
	    //if (con == 0) std::cout<<"this point is outside of TPC"<<std::endl;
	    //trj.emplace_back(::geoalgo::Vector(pt));
	    trj[i] = trk.LocationAtPoint(i);
	    
	  }
	  _matched_trj_v.emplace_back(trj);
	  _ctag_score = ctag.fCosmicScore;
	}	  
	++num_tagged;
	//
	// Store trajectory here
	//


	auto& trj = _rui_v[track_index];

	trj.reserve(trk.NumberTrajectoryPoints()+2);

	if(ctag.fCosmicScore == 1.0 &&_pt_up.IsValid())trj.emplace_back(_pt_up);
	if(ctag.fCosmicScore == 0.5 &&_pt_low.IsValid())trj.emplace_back(_pt_low);
	
	
	for(size_t i=0; i<trk.NumberTrajectoryPoints(); ++i) {
	  
	  auto const& pt = trk.LocationAtPoint(i);
	  trj.emplace_back(::geoalgo::Vector(pt[0],pt[1],pt[2]));
	  
	}
	
	if(ctag.fCosmicScore == 1.0 &&_pt_low.IsValid())trj.emplace_back(_pt_low);
	if(ctag.fCosmicScore == 0.5 &&_pt_up.IsValid())trj.emplace_back(_pt_up);
	
      }

      _temps1.push_back(temp1);
      _temps2.push_back(temp2);
      //if(ctag.fCosmicScore>0) std::cout<<"\033[93m"<<track_index<<"\033[00m"<<std::endl;
      //if(ev_ctag) ev_ctag->push_back(ctag);
      if(ev_ctag) (*ev_ctag)[track_index].fCosmicScore = ctag.fCosmicScore;
      _length = length;


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
