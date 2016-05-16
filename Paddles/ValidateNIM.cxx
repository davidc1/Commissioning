#ifndef LARLITE_VALIDATENIM_CXX
#define LARLITE_VALIDATENIM_CXX

#include "ValidateNIM.h"
#include "DataFormat/opdetwaveform.h"
#include "DataFormat/trigger.h"

int bitConversion(int bit)
{
  int powers = 0;
  while (bit != 1){
    bit /= 2;
    powers += 1;
  }
  return powers;
}


namespace larlite {

  ValidateNIM::ValidateNIM()
    : _tree(nullptr)
  {
    _name = "ValidateNIM";
    _fout = 0;

    std::cout << "1  = 2^" << bitConversion(1)  << std::endl;
    std::cout << "2  = 2^" << bitConversion(2)  << std::endl;
    std::cout << "4  = 2^" << bitConversion(4)  << std::endl;
    std::cout << "8  = 2^" << bitConversion(8)  << std::endl;
    std::cout << "16 = 2^" << bitConversion(16) << std::endl;

  }

  bool ValidateNIM::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("trigger_tree","Trigger Tree");
    _tree->Branch("_run",&_run,"run/I");
    _tree->Branch("_subrun",&_subrun,"subrun/I");
    _tree->Branch("_event",&_event,"event/I");
    _tree->Branch("_trig_num",&_trig_num,"trig_num/I");
    _tree->Branch("_trig_time",&_trig_time,"trig_time/D");
    _tree->Branch("_beam_time",&_beam_time,"beam_time/D");
    _tree->Branch("_trig_bits",&_trig_bits,"trig_bits/I");
    _tree->Branch("_strb1",&_strb1,"strb1/I");
    _tree->Branch("_strb2",&_strb2,"strb2/I");
    _tree->Branch("_numi",&_numi,"numi/I");
    _tree->Branch("_bnb",&_bnb,"bnb/I");
    _tree->Branch("_rwm",&_rwm,"rwm/I");
    _tree->Branch("_pdls",&_pdls,"pdls/I");
    _tree->Branch("_led",&_led,"led/I");
    _tree->Branch("_flash",&_flash,"flash/I");
    _tree->Branch("_dt",&_dt,"dt/D");
    _tree->Branch("_dt_same",&_dt_same,"dt_same/D");

    _dt = 0.;
    _t_rwm = _t_strb1 = _t_strb2 = _t_numi = _t_bnb = _t_led = _t_flash = _t_pdls = 0.;

    return true;
  }
  
  bool ValidateNIM::analyze(storage_manager* storage) {
  

    auto trig = storage->get_data<trigger>("daq");

    if (!trig){
      std::cout << "No trigger data product!" << std::endl;
      return false;
    }

    auto ev_wf = storage->get_data<event_opdetwaveform>("pmtreadout");
    
    if ( (!ev_wf) || (ev_wf->size() == 0) ){
      std::cout << "No waveforms saved!" << std::endl;
      return false;
    } 

    _run    = storage->run_id();
    _subrun = storage->subrun_id();
    _event  = storage->event_id();

    // save trigger information
    _trig_num  = trig->TriggerNumber();
    _dt        = trig->TriggerTime() - _trig_time;
    _trig_time = trig->TriggerTime();
    _beam_time = trig->BeamGateTime();
    _trig_bits = bitConversion(trig->TriggerBits());

    if (_verbose){
      std::cout << "Trigger Num  : " << _trig_num << std::endl
		<< "Trigger Time : " << _trig_time << std::endl
		<< "Beam Time    : " << _beam_time << std::endl
		<< "Trigger bits : " << _trig_bits << std::endl
		<< std::endl;
    }


    // reset trigger info
    resetTriggers();


    // now look at the waveforms saved
    for (size_t n=0; n < ev_wf->size(); n++){
      
      auto const& wf = ev_wf->at(n);

      auto const& ch = wf.ChannelNumber();

      if ( ( (ch >= 36)  && (ch < 100) ) ||
	   ( (ch >= 136) && (ch < 200) ) )
	{
	  if (_verbose)
	    std::cout << "looking at channel " << ch << "\t wf size : " << wf.size() << std::endl;
	  
	  for (size_t i=0; i < wf.size(); i++){
	    if ( wf[i] > 2200 ){
	      if (_verbose)
		std::cout << "found a pulse!" << std::endl;
	      if ( ch ==  39 ) 
		{ 
		  _rwm   = 1; 
		  _dt_same = _trig_time - _t_rwm;;
		  _t_rwm = _trig_time;
		}
	      if ( ch ==  38 )
		{
		  _strb1 = 1;
		  _dt_same = _trig_time - _t_strb1;
		  _t_strb1 = _trig_time;
		}
	      if ( ch ==  37 )
		{
		  _numi  = 1;
		  _dt_same = _trig_time - _t_numi;
		  _t_numi = _trig_time;
		}
	      if ( ch ==  36 )
		{
		  _bnb   = 1; 
		  _dt_same = _trig_time - _t_bnb;
		  _t_bnb = _trig_time;
		}
	      if ( ch == 139 )
		{
		  _flash = 1; 
		  _dt_same = _trig_time - _t_flash;
		  _t_flash = _trig_time;
		}
	      if ( ch == 138 )
		{
		  _strb2 = 1; 
		  _dt_same = _trig_time - _t_strb2;
		  _t_strb2 = _trig_time;
		}
	      if ( ch == 137 )
		{
		  _pdls  = 1; 
		  _dt_same = _trig_time - _t_pdls;
		  _t_pdls = _trig_time;
		}
	      if ( ch == 136 )
		{
		  _led   = 1; 
		  _dt_same = _trig_time - _t_led;
		  _t_led = _trig_time;
		}
	      break;
	    }
	  }

	}

    }// for all channels

    _tree->Fill();

    return true;
  }

  bool ValidateNIM::finalize() {

    if (_fout){
      if (_tree) { _tree->Write(); }
    }

    return true;
  }

  void ValidateNIM::resetTriggers() {

    _rwm = 0;
    _strb1 = 0;
    _strb2 = 0;
    _bnb = 0;
    _numi = 0;
    _led = 0;
    _pdls = 0;
    _flash = 0;

  }

}
#endif
