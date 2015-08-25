#ifndef LARLITE_COSMICDISCRIMFIFO_CXX
#define LARLITE_COSMICDISCRIMFIFO_CXX

#include "CosmicDiscrimFIFO.h"
#include "DataFormat/opdetwaveform.h"
#include "DataFormat/fifo.h"
#include "DataFormat/trigger.h"

namespace larlite {

  bool CosmicDiscrimFIFO::initialize() {

    if (_tree) delete _tree;
    
    _ev = 0;

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
    _tree->Branch("_ev_frame",&_ev_frame,"ev_frame/I");
    _tree->Branch("_frame_diff",&_frame_diff,"frame_diff/I");

    _trig_tree = new TTree("_trig_tree","Trigger Tree");
    _trig_tree->Branch("_trig_num",&_trig_num,"trig_num/I");
    _trig_tree->Branch("_trig_time",&_trig_time,"trig_time/D");
    _trig_tree->Branch("_delta_t",&_delta_t,"delta_t/D");
    _trig_tree->Branch("_frame_diff",&_frame_diff,"frame_diff/I");

    return true;
  }
  
  bool CosmicDiscrimFIFO::analyze(storage_manager* storage) {


    auto const trig    = storage->get_data<trigger>("daq");
    auto const ev_fifo = storage->get_data<event_fifo>("pmtreadout");

    if (!ev_fifo or !trig){
      std::cout << "nothing here..." << std::endl;
      return true;
    }

    if (ev_fifo->size() == 0){
      std::cout << "no data in fifo..." << std::endl;
      return true;
    }

    _ev_frame = ev_fifo->event_frame_number();
    _frame_diff = _ev_frame - _last_frame_num;
    _last_frame_num = _ev_frame;

    _trig_num = trig->TriggerNumber();
    _trig_time = trig->TriggerTime();
    _delta_t   = _trig_time-_last_trig_time;
    _last_trig_time = _trig_time;
    _trig_tree->Fill();




    for (size_t i=0; i < ev_fifo->size(); i++){

      auto const& wf = ev_fifo->at(i);

      _ch     = wf.channel_number();
      _disc   = wf.disc_id();
      _adcs   = wf.size();
      _adc_v  = wf;
      _frame  = wf.readout_frame_number();
      _sample = wf.readout_sample_number_RAW();
      _frame_diff = _frame-_ev_frame;
      _tree->Fill();
      
    }// for all waveforms in event

    _ev += 1;

    std::cout << "event: " << _ev << std::endl;
  
    return true;
  }

  bool CosmicDiscrimFIFO::finalize() {

    if (_fout and _tree)
      _tree->Write();
    if (_fout and _trig_tree)
      _trig_tree->Write();
  
    return true;
  }

}
#endif
