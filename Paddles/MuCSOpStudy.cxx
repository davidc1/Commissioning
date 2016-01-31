#ifndef LARLITE_MUCSOPSTUDY_CXX
#define LARLITE_MUCSOPSTUDY_CXX

#include "MuCSOpStudy.h"
#include "FhiclLite/ConfigManager.h"
#include "DataFormat/track.h"
#include "DataFormat/ophit.h"
#include "DataFormat/opflash.h"
#include "DataFormat/cosmictag.h"
namespace larlite {

  MuCSOpStudy::MuCSOpStudy() : ana_base()
  {
    _name="MuCSOpStudy";
    _fout=nullptr;
    _mgr.AddCustomAlgo(&_lpath);
    _mgr.SetAlgo(&_fhypo);
    _mgr.SetAlgo(&_qll);
  }
  
  void MuCSOpStudy::configure(const std::string config_file) {

    ::fcllite::ConfigManager cfg_mgr(_name);
    
    cfg_mgr.AddCfgFile(config_file);
    
    auto const& main_cfg = cfg_mgr.Config().get_pset(_name);
    
    _ophit_tmin = main_cfg.get<double>("OpHitTMin");
    _ophit_tmax = main_cfg.get<double>("OpHitTMax");

    _ophit_producer   = main_cfg.get<std::string>("OpHitProducer"    );
    _opflash_producer = main_cfg.get<std::string>("OpFlashProducer"  );
    _track_producer   = main_cfg.get<std::string>("TrackProducer"    );
    _ctag_producer    = main_cfg.get<std::string>("CosmicTagProducer");

    _run_match = main_cfg.get<bool>("RunMatch");
    _use_ophit_flash = main_cfg.get<bool>("UseOpHitFlash");
  }
  
  bool MuCSOpStudy::initialize() {

    _hRatioMap = new TH2D("hRatioMap","PE Assym. per PMT;Channel;Assym",
			  32,-0.5,31.5,
			  100,0,2);
    _hHitFlashScore = new TH1D("hHitFlashScore","Cheat Flash Score; Score; Cheat Flash",
			       100,0,1);
    _hMatchTime = new TH1D("hMatchTime","Matched Time;Time;Match",
			   100, -100,100);
    _hMatchScore = new TH1D("hMatchScore","Matched Score;Score;Match",
			    100,0,1);
    _hMatchScorePE = new TH2D("hMatchScorePE","Match PE vs Score;PE;Score",
			      100,0,2000,100,0,1);
    _hMatchScoreTime = new TH2D("hMatchScoreTime","Match Time vs Score:Time;Score",
				100,-100,100,100,0,1);
			      
    return true;
  }
  
  bool MuCSOpStudy::analyze(storage_manager* storage) {
  
    auto ev_ophit   = storage->get_data<event_ophit>(_ophit_producer);
    auto ev_opflash = storage->get_data<event_opflash>(_opflash_producer);
    auto ev_ctag    = storage->get_data<event_cosmictag>(_ctag_producer);
    auto ev_track   = storage->get_data<event_track>(_track_producer);

    if(!ev_ophit || !ev_ctag || !ev_track || !ev_opflash) {
      print(msg::kERROR,__FUNCTION__,"Missing data product!");
      std::cout<<"ophit:   "<<ev_ophit<<std::endl
	       <<"opflash: "<<ev_opflash<<std::endl
	       <<"ctag:    "<<ev_ctag<<std::endl
	       <<"track:   "<<ev_track<<std::endl<<std::endl;
      throw std::exception();
    }
    if(ev_ctag->size() != ev_track->size()) {
      print(msg::kERROR,__FUNCTION__,"Un-matched cosmic tag & track count!");
      throw std::exception();
    }

    _flash_v.clear();
    for(size_t flash_index=0; flash_index<ev_opflash->size(); ++flash_index) {

      auto const& opf = (*ev_opflash)[flash_index];
      ::flashana::Flash_t f;
      f.x = f.x_err = 0;
      f.y = opf.YCenter();
      f.z = opf.ZCenter();
      f.y_err = opf.YWidth();
      f.z_err = opf.ZWidth();
      f.pe_v.resize(32);
      for(size_t ch=0; ch<f.pe_v.size(); ++ch) 
	f.pe_v[ch] = opf.PE(ch);
      f.idx = flash_index;
      if(_run_match) _mgr.Add(f);
      _flash_v.emplace_back(f);
    }

    _ophit_flash = ::flashana::Flash_t();
    _ophit_flash.pe_v.resize(32);
    for(auto const& oph : *ev_ophit) {
      
      if(oph.OpChannel() >= _ophit_flash.pe_v.size()) continue;

      if(oph.PeakTime() < _ophit_tmin || oph.PeakTime() > _ophit_tmax) continue;
      
      _ophit_flash.pe_v[oph.OpChannel()] += oph.PE();

      _ophit_flash.time += oph.PeakTime() * oph.PE();
    }
    double ophit_qtot=0;
    for(auto const& v : _ophit_flash.pe_v) ophit_qtot += v;
    _ophit_flash.time /= ophit_qtot;
    _ophit_flash.idx = _flash_v.size();
    _flash_v.push_back(_ophit_flash);
    if(_run_match && _use_ophit_flash) _mgr.Add(_ophit_flash);

    _qcluster_v.clear();
    _cand_trj_v.clear();
    for(size_t track_index=0; track_index<ev_track->size(); ++track_index) {

      auto const& trk  = (*ev_track)[track_index];
      auto const& ctag = (*ev_ctag)[track_index];

      if(ctag.fCosmicScore<=0) continue;

      ::geoalgo::Trajectory trj;
      trj.reserve(trk.NumberTrajectoryPoints());

      for(size_t i=0; i<trk.NumberTrajectoryPoints(); ++i) {

	auto const& pt = trk.LocationAtPoint(i);
	trj.emplace_back(pt[0],pt[1],pt[2]);

      }
      
      auto qc = _lpath.FlashHypothesis(trj);

      ::flashana::Flash_t hypo;
      hypo.pe_v.resize(32);
      _fhypo.FillEstimate(qc,hypo);
      _hHitFlashScore->Fill(_qll.QLL(hypo,_ophit_flash));

      for(size_t ch=0; ch<hypo.pe_v.size(); ++ch)

	_hRatioMap->Fill(ch, (_ophit_flash.pe_v[ch] - hypo.pe_v[ch]) / (_ophit_flash.pe_v[ch] + hypo.pe_v[ch]) * 2.);

      if(_run_match) _mgr.Add(qc);
      _qcluster_v.emplace_back(qc);
      _cand_trj_v.emplace_back(trj);
    }

    if(_run_match) {
      auto res = _mgr.Match();
      for(auto const& match : res) {

	auto const& f = _flash_v[match.flash_id];

	_hMatchTime->Fill(f.time);
	_hMatchScore->Fill(match.score);
	double pe_tot=0;
	for(auto const& v : f.pe_v) pe_tot += v;
	_hMatchScorePE->Fill(pe_tot,match.score);
	_hMatchScoreTime->Fill(f.time,match.score);
	
      }
    }
      
    return true;
  }

  bool MuCSOpStudy::finalize() {

    if(_fout) {
      _hRatioMap->Write();
      _hHitFlashScore->Write();
      _hMatchTime->Write();
      _hMatchScore->Write();
      _hMatchScorePE->Write();
      _hMatchScoreTime->Write();
    }

    return true;
  }

}
#endif
