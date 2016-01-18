#ifndef LARLITE_MAKETRIGINFOTREE_CXX
#define LARLITE_MAKETRIGINFOTREE_CXX

#include "MakeTrigInfoTree.h"
#include "DataFormat/fifo.h"
#include "DataFormat/trigger.h"

namespace larlite {

  MakeTrigInfoTree::MakeTrigInfoTree()
    : _tree(nullptr)
  {
    _name = "MakeTrigInfoTree";
    _fout = 0;
  }

  bool MakeTrigInfoTree::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("tree","tree");
    _tree->Branch("_event",&_event,"event/I");
    _tree->Branch("_event_frame_num",&_event_frame_num,"event_frame_num/I");
    _tree->Branch("_fem_trig_sample_number_RAW",&_fem_trig_sample_number_RAW,"fem_trig_sample_number_RAW/I");
    _tree->Branch("ch0","std::vector<unsigned short>",&_ch0);
    _tree->Branch("ch1","std::vector<unsigned short>",&_ch1);
    _tree->Branch("_trig_time",&_trig_time,"trig_time/D");
    
    return true;
  }
  
  bool MakeTrigInfoTree::analyze(storage_manager* storage) {


    auto trig = storage->get_data<trigger>("daq");

    if (!trig){
      std::cout << "No trigger data product!" << std::endl;
      return false;
    }

    auto ev_fifo = storage->get_data<event_fifo>("pmt_xmit");

    if (!ev_fifo){
      std::cout << "No fifo data product!" << std::endl;
      return false;
    }

    _event = storage->event_id();

    _trig_time = trig->TriggerTime();

    _event_frame_num            = ev_fifo->event_frame_number();
    _fem_trig_sample_number_RAW = ev_fifo->fem_trig_sample_number_RAW();
    

    // find channels 0 and 1
    for (size_t i=0; i < ev_fifo->size(); i++){
      
      auto const& fifo = ev_fifo->at(i);
      
      if (fifo.channel_number() > 1)
	continue;
      
      if (fifo.channel_number() == 0)
	_ch0 = fifo;
      if (fifo.channel_number() == 1)
	_ch1 = fifo;
    }      

    _tree->Fill();
    
    return true;
  }

  bool MakeTrigInfoTree::finalize() {


    if (_fout)
      if (_tree) _tree->Write();
    
    return true;
  }

}
#endif
