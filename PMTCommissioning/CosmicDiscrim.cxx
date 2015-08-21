#ifndef LARLITE_COSMICDISCRIM_CXX
#define LARLITE_COSMICDISCRIM_CXX

#include "CosmicDiscrim.h"
#include "DataFormat/opdetwaveform.h"

namespace larlite {

  bool CosmicDiscrim::initialize() {

    if (_tree) delete _tree;
    
    _ev = 0;

    _tree = new TTree("_tree","opdetwf tree");
    _tree->Branch("_ch",&_ch,"ch/I");
    _tree->Branch("_ev",&_ev,"ev/I");
    _tree->Branch("_adc_v","std::vector<short>",&_adc_v);
    _tree->Branch("_adcs",&_adcs,"adcs/I");
    _tree->Branch("_tstart",&_tstart,"tstart/D");
    _tree->Branch("_tend",&_tend,"tend/D");


    return true;
  }
  
  bool CosmicDiscrim::analyze(storage_manager* storage) {

    auto const ev_opdetwf = storage->get_data<event_opdetwaveform>("pmtreadout");

    for (size_t i=0; i < ev_opdetwf->size(); i++){

      auto const& wf = ev_opdetwf->at(i);

      _ch     = wf.ChannelNumber();
      _adcs   = wf.size();
      _adc_v  = wf;
      _tstart = wf.TimeStamp();
      _tend   = _tstart + _adcs;
      _tree->Fill();
      
    }// for all waveforms in event
    
    
    _ev += 1;
    
    return true;
  }

  bool CosmicDiscrim::finalize() {

    if (_fout and _tree)
      _tree->Write();

    return true;
  }

}
#endif
