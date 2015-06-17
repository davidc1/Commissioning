#ifndef SIMPLEWFANA_CXX
#define SIMPLEWFANA_CXX

#include "SimpleWFAna.h"

namespace larlite {

  bool SimpleWFAna::initialize() {

    _evtN = 0;

    _t_ch = new TTree("ch_tree","");
    _t_ch->Branch("larch",&_larch,"larch/i");
    _t_ch->Branch("rms",&_rms,"rms/F");
    _t_ch->Branch("evt",&_evtN,"evt/i");
    _t_ch->Branch("mean",&_mean,"mean/F");
    _t_ch->SetDirectory(0);

    return true;
  }
  
  bool SimpleWFAna::analyze(storage_manager* storage) {

    auto wfs = storage->get_data<event_rawdigit>("daq");

    // Make a map from channel index to channel number
    for (size_t i=0; i < wfs->size(); i++){

      auto const& wf = (*wfs).at(i);

      int chnum = wf.Channel();
      _larch = chnum;
      
      auto const& adcs = wf.ADCs();

      
      // Mean
      for(size_t j=0; j<adcs.size(); ++j)
	_mean += adcs[j];
      _mean /= ((float)adcs.size());
      // RMS
      for(size_t k=0; k<adcs.size(); ++k)
	_rms += (adcs[k]-_mean)*(adcs[k]-_mean);
      _rms = sqrt( _rms / ((float)adcs.size()));
      
      _t_ch->Fill();	

    }

    _evtN += 1;
    
    return true;
  }

  bool SimpleWFAna::finalize() {
    
    if(_fout){
      _fout->cd();
      std::cout << "writing ch tree" << std::endl;
      _t_ch->Write();
    }
    return true;
  }

}
#endif
