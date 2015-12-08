#ifndef LARLITE_MAKETRIGINFOTREE_CXX
#define LARLITE_MAKETRIGINFOTREE_CXX

#include "MakeTrigInfoTree.h"
#include "DataFormat/fifo.h"
#include "DataFormat/trigger.h"

namespace larlite {

  bool MakeTrigInfoTree::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("tree","tree");
    _tree->Branch("_event",&_event,"event/I");
    _tree->Branch("ch0","std::vector<short>",&_ch0);
    _tree->Branch("ch1","std::vector<short>",&_ch1);
    _tree->Branch("_trig_time",&_trig_time,"trig_time/D");
    
    return true;
  }
  
  bool MakeTrigInfoTree::analyze(storage_manager* storage) {


    auto trig = storage->get_data<trigger>("daq");

    if (!trig){
      std::cout << "No trigger data product!" << std::endl;
      return false;
    }

    auto fifo = storage->get_data<event_fifo>("pmt_xmit");

    if (!fifo){
      std::cout << "No fifo data product!" << std::endl;
      return false;
    }

    _event = storage->event_id();
    
  
    return true;
  }

  bool MakeTrigInfoTree::finalize() {


    if (_fout)
      if (_tree) _tree->Write();
    
    return true;
  }

}
#endif
