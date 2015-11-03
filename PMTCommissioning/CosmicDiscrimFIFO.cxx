#ifndef LARLITE_COSMICDISCRIMFIFO_CXX
#define LARLITE_COSMICDISCRIMFIFO_CXX

#include "CosmicDiscrimFIFO.h"
#include "DataFormat/opdetwaveform.h"
#include "DataFormat/fifo.h"
#include "DataFormat/trigger.h"

#include <climits>
#include <limits>
#include <cmath>

namespace larlite {

  CosmicDiscrimFIFO::CosmicDiscrimFIFO()
    : _trig_tree(nullptr)
    , _tree(nullptr)
    , _rate_tree(nullptr)
  {
    _name="CosmicDiscrimFIFO";
    _fout=0;
    _use_trig=false;
    _save_wf=false;
    _verbose=false;
    _trig_producer = "";
    _fifo_producer = "";
    _adc_thresh = 0;
    _hg = true;
  }

  bool CosmicDiscrimFIFO::initialize() {

    if (_tree) delete _tree;
    if (_rate_tree) delete _rate_tree;
    if (_trig_tree) delete _trig_tree;

    _ev = 0;

    _baselines.clear();
    _rms.clear();
    for (int i=0; i < 32; i++){
      _baselines.push_back(kMAX_DOUBLE);
      _rms.push_back(kMAX_DOUBLE);
    }

    _last_trig_time = 0;
    _last_frame_num = 0;

    _tree = new TTree("_tree","opdetwf tree");
    _tree->Branch("_ch",&_ch,"ch/I");
    _tree->Branch("_ev",&_ev,"ev/I");
    _tree->Branch("_frame",&_frame,"frame/I");
    _tree->Branch("_sample",&_sample,"sample/I");
    _tree->Branch("_disc",&_disc,"disc/I");
    _tree->Branch("_adc_v","std::vector<unsigned short>",&_adc_v);
    _tree->Branch("_adcs",&_adcs,"adcs/I");
    _tree->Branch("_amp",&_amp,"amp/I");
    _tree->Branch("_rms_ch",&_rms_ch,"rms_ch/D");
    _tree->Branch("_avg_ch",&_avg_ch,"avg_ch/D");
    _tree->Branch("_first",&_first,"first/I");
    _tree->Branch("_ev_frame",&_ev_frame,"ev_frame/I");
    _tree->Branch("_frame_diff",&_frame_diff,"frame_diff/I");

    _trig_tree = new TTree("_trig_tree","Trigger Tree");
    _trig_tree->Branch("_trig_num",&_trig_num,"trig_num/I");
    _trig_tree->Branch("_trig_time",&_trig_time,"trig_time/D");
    _trig_tree->Branch("_delta_t",&_delta_t,"delta_t/D");
    _trig_tree->Branch("_frame_diff",&_frame_diff,"frame_diff/I");

    _rate_tree = new TTree("_rate_tree","Rateger Tree");
    _rate_tree->Branch("_n20_windows_h","std::vector<unsigned short>",&_n20_windows_h);
    _rate_tree->Branch("_n20_windows_l","std::vector<unsigned short>",&_n20_windows_l);
    _rate_tree->Branch("_n1k_windows","std::vector<unsigned short>",&_n1k_windows);
    _rate_tree->Branch("_rms_v","std::vector<double>",&_rms_v);
    _rate_tree->Branch("_baseline_v","std::vector<double>",&_baseline_v);
    _rate_tree->Branch("_charge_v","std::vector<double>",&_charge_v);
    _rate_tree->Branch("_pulses_03","std::vector<unsigned short>",&_pulses_03);
    _rate_tree->Branch("_pulses_06","std::vector<unsigned short>",&_pulses_06);
    _rate_tree->Branch("_pulses_09","std::vector<unsigned short>",&_pulses_09);
    _rate_tree->Branch("_areas_09","std::vector<double>",&_areas_09);
    _rate_tree->Branch("_pulses_12","std::vector<unsigned short>",&_pulses_12);
    _rate_tree->Branch("_pulses_15","std::vector<unsigned short>",&_pulses_15);
    _rate_tree->Branch("_pulses_18","std::vector<unsigned short>",&_pulses_18);
    _rate_tree->Branch("_pulses_21","std::vector<unsigned short>",&_pulses_21);
    _rate_tree->Branch("_pulses_24","std::vector<unsigned short>",&_pulses_24);
    _rate_tree->Branch("_pulses_27","std::vector<unsigned short>",&_pulses_27);
    _rate_tree->Branch("_pulses_30","std::vector<unsigned short>",&_pulses_30);
    _rate_tree->Branch("_pulses_33","std::vector<unsigned short>",&_pulses_33);
    _rate_tree->Branch("_pulses_36","std::vector<unsigned short>",&_pulses_36);
    _rate_tree->Branch("_pulses_39","std::vector<unsigned short>",&_pulses_39);
    _rate_tree->Branch("_pulses_42","std::vector<unsigned short>",&_pulses_42);
    _rate_tree->Branch("_pulses_45","std::vector<unsigned short>",&_pulses_45);
    _rate_tree->Branch("_pulses_60","std::vector<unsigned short>",&_pulses_60);
    _rate_tree->Branch("_pulses_100","std::vector<unsigned short>",&_pulses_100);
    _rate_tree->Branch("_ev",&_ev,"ev/I");

    return true;
  }
  
  bool CosmicDiscrimFIFO::analyze(storage_manager* storage) {



    auto const ev_fifo = storage->get_data<event_fifo>(_fifo_producer);

    if (!ev_fifo){
      std::cout << "nothing here..." << std::endl;
      return true;
    }

    if (ev_fifo->size() == 0){
      std::cout << "no data in fifo..." << std::endl;
      return true;
    }

    // resize windows for the number of PMTs we want
    ClearVectors();

    // if we want to use trigger
    if (_use_trig){
      auto const trig    = storage->get_data<trigger>(_trig_producer);
      if (!trig){
	std::cout << "nothing here..." << std::endl;
	return true;
      }
      _trig_num = trig->TriggerNumber();
      _trig_time = trig->TriggerTime();
      _delta_t   = _trig_time-_last_trig_time;
      _last_trig_time = _trig_time;
      _trig_tree->Fill();
    }// if we should use the trigger

    // get event information
    _ev_frame = ev_fifo->event_frame_number();
    _frame_diff = _ev_frame - _last_frame_num;
    _last_frame_num = _ev_frame;

    // find the 1500-sample window and use it to measure a baseline & rms value
    for (size_t i=0; i < ev_fifo->size(); i++){

      auto const& wf = ev_fifo->at(i);
      
      if ( (_hg == true)  &&  (wf.module_address() != 5) ) continue;
      if ( (_hg == false) &&  (wf.module_address() == 5) ) continue;

      _ch     = wf.channel_number();
      if (_ch >= 32) continue;

      if (wf.size() != 1000) continue;

      _n1k_windows[_ch] += 1;

      auto pedestal = GetBaselineRms(wf);
      _baseline_v[_ch] = pedestal.first;
      _rms_v[_ch] = pedestal.second;

      // find and keep track of SPE pulses in beam-gate window
      _pulses_03[_ch] = FindSPE(wf,pedestal.first,pedestal.second, 3);
      _pulses_06[_ch] = FindSPE(wf,pedestal.first,pedestal.second, 6);
      _pulses_09[_ch] = FindSPE(wf,pedestal.first,pedestal.second, 9);
      _pulses_12[_ch] = FindSPE(wf,pedestal.first,pedestal.second,12);
      _pulses_15[_ch] = FindSPE(wf,pedestal.first,pedestal.second,15);
      _pulses_18[_ch] = FindSPE(wf,pedestal.first,pedestal.second,18);
      _pulses_21[_ch] = FindSPE(wf,pedestal.first,pedestal.second,21);
      _pulses_24[_ch] = FindSPE(wf,pedestal.first,pedestal.second,24);
      _pulses_27[_ch] = FindSPE(wf,pedestal.first,pedestal.second,27);
      _pulses_30[_ch] = FindSPE(wf,pedestal.first,pedestal.second,30);
      _pulses_33[_ch] = FindSPE(wf,pedestal.first,pedestal.second,33);
      _pulses_36[_ch] = FindSPE(wf,pedestal.first,pedestal.second,36);
      _pulses_39[_ch] = FindSPE(wf,pedestal.first,pedestal.second,39);
      _pulses_42[_ch] = FindSPE(wf,pedestal.first,pedestal.second,42);
      _pulses_45[_ch] = FindSPE(wf,pedestal.first,pedestal.second,45);
      _pulses_60[_ch] = FindSPE(wf,pedestal.first,pedestal.second,60);
      _pulses_100[_ch] = FindSPE(wf,pedestal.first,pedestal.second,100,281);

      if (pedestal.second < _rms[_ch]){
	_baselines[_ch] = pedestal.first;
	_rms[_ch]       = pedestal.second;
      }// if we should update pedestal and rms
    }// for all waveforms in the envet


    if (_verbose){
      std::cout << "Baselines: " << std::endl
		<< "[";
      for (auto pmt : _baseline_v)
	std::cout << pmt << ", ";
      std::cout << "]" << std::endl;
      std::cout << "RMS: " << std::endl
		<< "[";
      for (auto pmt : _rms_v)
	std::cout << pmt << ", ";
      std::cout << "]" << std::endl;
    }

    for (size_t i=0; i < ev_fifo->size(); i++){

      auto const& wf = ev_fifo->at(i);

      if ( (_hg == true)  &&  (wf.module_address() != 5) ) continue;
      if ( (_hg == false) &&  (wf.module_address() == 5) ) continue;

      _adc_v.clear();
      _ch     = wf.channel_number();
      if (_ch >= 32) continue;
      _disc   = wf.disc_id();
      _adcs   = wf.size();

      //if (_adcs != 20) continue;

      _frame  = wf.readout_frame_number();
      _sample = wf.readout_sample_number_RAW();
      _frame_diff = _frame-_ev_frame;

      _rms_ch = _rms_v[_ch];
      _avg_ch = _baseline_v[_ch];

      _first = (int)wf[0];

      _amp = _first;
      for (auto& adc : wf)
	if (adc > _amp) { _amp = (int)adc; }

      // cut on amplitude of first tick
      if ( fabs(_first-_baselines[_ch]) < 3*_rms[_ch] ){
	if ((_amp-_baselines[_ch]) > _adc_thresh){
	  // baseline subtraction
	  _n20_windows_h[_ch] += 1;
	  _charge_v[_ch] += _amp-_baselines[_ch];
	}
	else
	  _n20_windows_l[_ch] += 1;
      }

      if (_save_wf)
	_adc_v  = wf;

      _tree->Fill();

    }// for all waveforms in event

    _ev += 1;

    _rate_tree->Fill();

    return true;
  }

  bool CosmicDiscrimFIFO::finalize() {

    if (_fout){
      if (_tree)
	_tree->Write();
      if (_trig_tree)
	_trig_tree->Write();
      if (_rate_tree)
	_rate_tree->Write();
    }
  
    return true;
  }

  std::pair<double,double> CosmicDiscrimFIFO::GetBaselineRms(const std::vector<unsigned short>& wf)
  {

    if (wf.size() == 0)
      return std::pair<double,double>(kMAX_DOUBLE,kMAX_DOUBLE);

    double avg = 0;
    double rms = 0;
    
    for (auto adc : wf)
      avg += adc;
    avg /= (double)wf.size();

    for (auto adc : wf)
      rms += (adc-avg)*(adc-avg);
    rms = sqrt( rms / ((double)wf.size() - 1) );
    
    return std::pair<double,double>(avg,rms);
  }
  
  unsigned short CosmicDiscrimFIFO::FindSPE(const std::vector<unsigned short>& wf,
					    const double& avg, const double& rms,
					    const double& thresh,
					    const int& deadtime)
  {
    
    bool active = false;
    double max  = kMIN_DOUBLE;
    
    unsigned short pulses = 0;
    int time = 0;
    int lasttime = 0;

    int cntr = 0;
    for (auto adc : wf){

      cntr += 1;
      
      if ( (adc-avg) > thresh){
	active = true;
	if ( (adc-avg) > max){
	  max = (adc-avg);
	  time = cntr;
	}
      }// if in pulse region
      else{
	if (active == true){
	  if ( (time-lasttime) > deadtime){
	    pulses += 1;
	    lasttime = time;
	  }
	  max = kMIN_DOUBLE;
	  active = false;
	}// if previously in active region
      }// if not in pulse region

    }// for all ADCs

    //std::cout << "found " << pulses << " pulses w/ " << thresh << "threshold!" << std::endl;
    return pulses;
  }

  void CosmicDiscrimFIFO::ClearVectors()
  {
    _n20_windows_h.clear();
    _n20_windows_h.resize(32);
    _n20_windows_l.clear();
    _n20_windows_l.resize(32);
    _n1k_windows.clear();
    _n1k_windows.resize(32);
    _rms_v.clear();
    _rms_v.resize(32);
    _baseline_v.clear();
    _baseline_v.resize(32);
    _charge_v.clear();
    _charge_v.resize(32);
    _pulses_03.clear();
    _pulses_03.resize(32);
    _pulses_06.clear();
    _pulses_06.resize(32);
    _pulses_09.clear();
    _pulses_09.resize(32);
    _pulses_12.clear();
    _pulses_12.resize(32);
    _pulses_15.clear();
    _pulses_15.resize(32);
    _pulses_18.clear();
    _pulses_18.resize(32);
    _pulses_21.clear();
    _pulses_21.resize(32);
    _pulses_24.clear();
    _pulses_24.resize(32);
    _pulses_27.clear();
    _pulses_27.resize(32);
    _pulses_30.clear();
    _pulses_30.resize(32);
    _pulses_33.clear();
    _pulses_33.resize(32);
    _pulses_36.clear();
    _pulses_36.resize(32);
    _pulses_39.clear();
    _pulses_39.resize(32);
    _pulses_42.clear();
    _pulses_42.resize(32);
    _pulses_45.clear();
    _pulses_45.resize(32);
    _pulses_60.clear();
    _pulses_60.resize(32);
    _pulses_100.clear();
    _pulses_100.resize(32);

    return;
  }


}
#endif
