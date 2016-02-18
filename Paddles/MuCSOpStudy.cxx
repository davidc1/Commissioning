#ifndef LARLITE_MUCSOPSTUDY_CXX
#define LARLITE_MUCSOPSTUDY_CXX

#include "MuCSOpStudy.h"
#include "FhiclLite/ConfigManager.h"
#include "DataFormat/track.h"
#include "DataFormat/ophit.h"
#include "DataFormat/opflash.h"
#include "DataFormat/cosmictag.h"
//#include "mytimer.h"
namespace larlite {

  MuCSOpStudy::MuCSOpStudy() : ana_base()
  {
    _name="MuCSOpStudy";
    _fout=nullptr;
    _mgr.AddCustomAlgo(&_lpath);
    _mgr.SetAlgo(&_qll);
    _mgr.SetAlgo(&_fhypo);
  }
  
  void MuCSOpStudy::configure(const std::string config_file) {

    ::fcllite::ConfigManager cfg_mgr(_name);
    
    cfg_mgr.AddCfgFile(config_file);
    
    auto const& main_cfg = cfg_mgr.Config().get_pset(_name);

    _if_gain    = main_cfg.get<bool>("IfGain");
    _num_ch     = main_cfg.get<size_t>("NumChannels");
    _ophit_veto = main_cfg.get<double>("OpHitVetoPeriod");
    _ophit_tmin = main_cfg.get<double>("OpHitTMin");
    _ophit_tmax = main_cfg.get<double>("OpHitTMax");
    
    _disc_threshold   = main_cfg.get<double>("DiscThreshold");
    _ophit_producer   = main_cfg.get<std::string>("OpHitProducer"    );
    _opflash_producer = main_cfg.get<std::string>("OpFlashProducer"  );
    _track_producer   = main_cfg.get<std::string>("TrackProducer"    );
    _ctag_producer    = main_cfg.get<std::string>("CosmicTagProducer");
    _gain_correction  = main_cfg.get<std::vector<double>>("GainsCorrection");
        
    _run_match = main_cfg.get<bool>("RunMatch");
    _use_ophit_flash = main_cfg.get<bool>("UseOpHitFlash");
  }
  
  bool MuCSOpStudy::initialize() {

    _hRatioMap = new TH2D("hRatioMap","Hypothesis / Optical Reco. per PMT;Channel;Ratio",
			  _num_ch,-0.5,_num_ch-0.5,
			  100,0,3);
    _hHitFlashScore = new TH1D("hHitFlashScore","Cheat Flash Score; Score; Cheat Flash",
			       100,0,1);
    _hRatioPLOP = new TH1D("hRatioPLOP","Hypothesis / Optical Reconstructionn",50,0,2);
    
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
    Watch my_watch;
    my_watch.Start();
    _mgr.Reset();
    _ctag_score = -1;
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
    if(ev_track->size() != MuCSTagger::GetME()->Trajectories().size()) {
      print(msg::kERROR,__FUNCTION__,"Un-matched trajectory & track count!");
      throw std::exception();
    }

    _flash_v.clear();
    //int i =0;
    for(size_t flash_index=0; flash_index<ev_opflash->size(); ++flash_index) {

      auto const& opf = (*ev_opflash)[flash_index];
      ::flashana::Flash_t f;
      f.x = f.x_err = 0;
      f.y = opf.YCenter();
      f.z = opf.ZCenter();
      f.y_err = opf.YWidth();
      f.z_err = opf.ZWidth();
      f.pe_v.resize(_num_ch);
      f.time = opf.Time();
      for(size_t ch=0; ch<f.pe_v.size(); ++ch) 
	f.pe_v[ch] = opf.PE(ch);
      f.idx = flash_index;
      if(_run_match && f.TotalPE()> 10) {_mgr.Add(f);} //++i;}
      _flash_v.emplace_back(f);
    }
    //std::cout<<"found "<<i<<" falshes !\n";
    //i=0;

    

    //
    // Figure out veto channels
    //
    _veto_v.resize(_num_ch,false);
    for(auto v : _veto_v) v = false;
    
    const double veto_start = _ophit_tmin - _ophit_veto;
    const double veto_end   = _ophit_tmin;
    for(auto const& oph : *ev_ophit) {
      
      if(oph.OpChannel() >= _veto_v.size()) continue;
      
      if(oph.PeakTime() < veto_start || oph.PeakTime() > veto_end) continue;

      _veto_v[oph.OpChannel()] = true;

    }

    _ophit_flash = ::flashana::Flash_t();
    _ophit_flash.time = 0;
    _ophit_flash.pe_v.resize(_num_ch);
    for(auto const& oph : *ev_ophit) {
            
      if(oph.OpChannel() >= _ophit_flash.pe_v.size()) continue;

      if(_veto_v[oph.OpChannel()]) continue;
      
      if(oph.PeakTime() < _ophit_tmin || oph.PeakTime() > _ophit_tmax) continue;      
      
      if(_if_gain)_ophit_flash.pe_v[oph.OpChannel()] += oph.PE()/_gain_correction[oph.OpChannel()];
      else _ophit_flash.pe_v[oph.OpChannel()] += oph.PE();
      //_ophit_flash.pe_v[oph.OpChannel()] += oph.Amplitude()/20;
      
      _ophit_flash.time += oph.PeakTime() * oph.PE();
    }
    
    double ophit_qtot=0;
    for(auto& v : _ophit_flash.pe_v) {
      v *= (1./0.23);
      ophit_qtot += v;
    }
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

      /*
      ::geoalgo::Trajectory trj;
      trj.reserve(trk.NumberTrajectoryPoints());

      for(size_t i=0; i<trk.NumberTrajectoryPoints(); ++i) {

	auto const& pt = trk.LocationAtPoint(i);
	trj.emplace_back(::geoalgo::Vector(pt[0],pt[1],pt[2]));

      }
      */
      auto const& trj = MuCSTagger::GetME()->Trajectories()[track_index];
      
      _ctag_score = ctag.fCosmicScore;
      //std::cout<<"time used before LP one event is "<<my_watch.WallTime()<<"s"<<std::endl;      
      auto qc = _lpath.FlashHypothesis(trj);
      //std::cout<<"time used by lp is "<<my_watch.WallTime()<<"s"<<std::endl;
      
      ::flashana::Flash_t hypo;
      hypo.pe_v.resize(_num_ch);
      _fhypo.FillEstimate(qc,hypo);
      _ophit_hypo = hypo;
      double pl_qtot=0;
      for(size_t ch=0; ch<hypo.pe_v.size(); ++ch) {
	if(_veto_v[ch]) continue;
	auto const& v = hypo.pe_v[ch];
	if(v>=_disc_threshold) pl_qtot += v;
      }

      _hRatioPLOP->Fill(pl_qtot/ophit_qtot);
      //_hRatioPLOP->Fill(ophit_qtot);
      //std::cout<<track_index<<std::endl;
      //double score = _qll.QLL(hypo,_ophit_flash);
      //_hHitFlashScore->Fill(score);
      //std::cout<<score<<std::endl;

      auto match = _qll.Match(qc,_ophit_flash);
      _hHitFlashScore->Fill(match.score);
      //std::cout<<match.score<<std::endl;

      for(size_t ch=0; ch<hypo.pe_v.size(); ++ch) {
	//std::cout<<_ophit_flash.pe_v[ch]<< " : "<<hypo.pe_v[ch]<<std::endl;
	if(_ophit_flash.pe_v[ch]<(5/0.23) || hypo.pe_v[ch]<(5/0.23)) continue;
	//_hRatioMap->Fill(ch, (_ophit_flash.pe_v[ch] - hypo.pe_v[ch]) / (_ophit_flash.pe_v[ch] + hypo.pe_v[ch]) * 2.);
	_hRatioMap->Fill(ch, (( hypo.pe_v[ch]) / _ophit_flash.pe_v[ch] ) );
      }
      if(_run_match) _mgr.Add(qc);
      _qcluster_v.emplace_back(qc);
      _cand_trj_v.emplace_back(trj);
    }

    //std::cout<<"time for lp on all tracks and match score "<<my_watch.WallTime()<<"s"<<std::endl;
    
    
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
      
      //std::cout<<"time for 4 other histos on match  "<<my_watch.WallTime()<<"s"<<std::endl;
      
      //std::cout<<"time for used event is "<<my_watch.WallTime()<<"s"
      //<<"\n"<<"===================================="<<std::endl;
      
    return true;
  }

  bool MuCSOpStudy::finalize() {

    if(_fout) {
      _hRatioMap->Write();
      _hRatioPLOP->Write();
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
