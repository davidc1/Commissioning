#ifndef LARLITE_FILLTREE_CXX
#define LARLITE_FILLTREE_CXX

#include "FillTree.h"
#include "DataFormat/trigger.h"
#include "DataFormat/opflash.h"

namespace larlite {

  FillTree::FillTree()
    : _tree(nullptr)
  {
    _name = "FillTree";
    _fout = 0;
    _flashProducer = "opflashSat";
    _pe_min = 50;
    _minT   = -10;
    _maxT   =  30;
  }

  bool FillTree::initialize() {

    if (_tree) { delete _tree; }
    _tree = new TTree("_tree","Flash Tree");
    _tree->Branch("_event",&_event,"event/I");
    _tree->Branch("_run",&_run,"run/I");
    _tree->Branch("_subrun",&_subrun,"subrun/I");
    _tree->Branch("_dt",&_dt,"dt/D");
    _tree->Branch("_flash_time",&_flash_time,"flash_time/D");
    _tree->Branch("_flash_abstime",&_flash_abstime,"flash_abstime/D");
    _tree->Branch("_trig_time",&_trig_time,"trig_time/D");
    _tree->Branch("_pe_total",&_pe_total,"pe_total/D");
    _tree->Branch("_trig_word",&_trig_word,"trig_word/I");

    return true;
  }
  
  bool FillTree::analyze(storage_manager* storage) {

    auto trig = storage->get_data<trigger>("daq");

    auto ev_flash = storage->get_data<event_opflash>(_flashProducer);

    if (!trig){
      print(msg::kERROR,__FUNCTION__,"Missing Trigger data product!");
      throw std::exception();
    }

    if (!ev_flash){
      print(msg::kERROR,__FUNCTION__,"Missing OpFlash data product!");
      throw std::exception();
    }

    _run    = storage->run_id();
    _subrun = storage->subrun_id();
    _event  = storage->event_id();

    _trig_word = trig->TriggerBits();

    _trig_time = trig->TriggerTime();

    // choose only max PE flash per event
    double max_pe = 0;
    double time = 0;

    for (size_t i=0; i < ev_flash->size(); i++){

      auto const& flash = ev_flash->at(i);

      _pe_total = flash.TotalPE();

      _pe_total = 0;
      for (size_t pmt=0; pmt < 32; pmt++)
	_pe_total += flash.PE(pmt);

      if (_pe_total < _pe_min)
	continue;

      _flash_time    = flash.Time();
      _flash_abstime = flash.AbsTime();

      _dt = flash.AbsTime() - _trig_time;

      _tree->Fill();

      if ( (_dt < _minT) or (_dt > _maxT) )
	continue;

      if (flash.TotalPE() > max_pe){
	max_pe = flash.TotalPE();
	time   = _flash_abstime;
      }

    }// for all flashes

    /*
    if (max_pe != 0){
      _flash_abstime = time;
      _dt = _flash_abstime - _trig_time;
      _pe_total = max_pe;
      _tree->Fill();
    }
    */

    //}// for all flashes
    
    return true;
  }

  bool FillTree::finalize() {

    if (_fout){
      _fout->cd();
      _tree->Write();
    }

    return true;
  }

}
#endif
