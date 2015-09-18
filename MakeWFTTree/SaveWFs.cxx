#ifndef LARLITE_SAVEWFS_CXX
#define LARLITE_SAVEWFS_CXX

#include "SaveWFs.h"
#include <algorithm>

namespace larlite {

  SaveWFs::SaveWFs()
    : _tree(nullptr)
  {
    _name = "SaveWFs";
    _fout = 0;
  }

  bool SaveWFs::initialize() {

    if (_tree) { delete _tree; }
    _tree = new TTree("tree","Waveform Tree");
    _tree->Branch("_ch",&_ch,"ch/I");
    _tree->Branch("_pl",&_pl,"pl/I");
    _tree->Branch("_evt",&_evt,"evt/I");
    _tree->Branch("_wire",&_wire,"wire/I");
    _tree->Branch("_ADCs","std::vector<short>",&_ADCs);

    std::cout << "Saving waveform for the following channels:" << std::endl;
    for (auto const& ch : _chList)
      std::cout << ch << std::endl;
    std::cout << std::endl;

    return true;
  }
  
  bool SaveWFs::analyze(storage_manager* storage) {
  
    _evt += 1;

    // read in RawDigits
    auto const ev_wf   = storage->get_data<event_rawdigit>("daq");

    for (size_t n=0; n < ev_wf->size(); n++){

      auto const wf = ev_wf->at(n);
      
      auto ch         = wf.Channel();
      auto const view = larutil::Geometry::GetME()->View(ch);
      auto const pl   = larutil::Geometry::GetME()->ChannelToPlane(ch);
      auto wire       = larutil::Geometry::GetME()->ChannelToWireID(ch);

      // is the channel in the channel-list?
      if ( find(_chList.begin(),_chList.end(),ch) != _chList.end() ){
	_ADCs = wf.ADCs();
	_ch   = ch;
	_pl   = pl;
	_wire = wire.Wire;
	_tree->Fill();
      }

    }// for all waveforms
	

    return true;
  }

  bool SaveWFs::finalize() {

    if (_fout)
      _tree->Write();

    return true;
  }

}
#endif
