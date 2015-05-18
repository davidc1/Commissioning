#ifndef ZIGZAG_CXX
#define ZIGZAG_CXX

#include "ZigZag.h"
#include <TStopwatch.h>
namespace larlite {

  void ZigZag::SetMap(unsigned int lar_ch,
		      short crate,
		      short slot,
		      short femch)
  {
    if(_ch_to_crate.size() <= lar_ch) {
      _ch_to_crate.resize(lar_ch+1,-1);
      _ch_to_slot.resize(lar_ch+1,-1);
      _ch_to_femch.resize(lar_ch+1,-1);
    }
    _ch_to_crate[lar_ch] = crate;
    _ch_to_slot[lar_ch]  = slot;
    _ch_to_femch[lar_ch] = femch;
  }

  bool ZigZag::initialize() {

    _tree_v.resize((*_ref_channel_s.rbegin())+1,nullptr);
    _diff_v.resize(_nsamples);

    _evtN = 0;

    int counter_tree = 0;

    _t_ch = new TTree("ch_tree","");
    _t_ch->Branch("larch",&_larch,"larch/i");
    _t_ch->Branch("crate",&_crate,"crate/i");
    _t_ch->Branch("slot",&_slot,"slot/i");
    _t_ch->Branch("femch",&_femch,"femch/i");
    _t_ch->Branch("rms",&_rms,"rms/F");
    _t_ch->Branch("evt",&_evtN,"evt/i");
    _t_ch->Branch("mean",&_mean,"mean/F");
    _t_ch->Branch("block",&_block,"block/i");
    _t_ch->Branch("adc_v","std::vector<short>",&_adc_v);
    _t_ch->Branch("zigzag_start","std::vector<short>",&_zigzag_start);
    _t_ch->Branch("zigzag_end","std::vector<short>",&_zigzag_end);
    _t_ch->SetDirectory(0);

    _t_zig = new TTree("zig_tree","");
    _t_zig->Branch("evt",&_evtN,"evt/i");
    _t_zig->Branch("bursts","std::vector<int>",&_bursts);
    _t_zig->Branch("frac",&_frac,"frac/F");

    _t_z = new TTree("z_tree","");
    _t_z->Branch("evt",&_evtN,"evt/i");
    _t_z->Branch("start",&_start,"start/i");
    _t_z->Branch("end",&_end,"end/i");
    _t_z->Branch("len",&_len,"len/i");
    _t_z->Branch("amp",&_amp,"amp/f");
    _t_z->Branch("larch",&_larch,"larch/i");
    _t_z->Branch("crate",&_crate,"crate/i");
    _t_z->Branch("slot",&_slot,"slot/i");
    _t_z->Branch("femch",&_femch,"femch/i");

    _t_corr = new TTree("corr_tree","");
    _t_corr->Branch("ref_larch",&_ref_larch,"ref_larch/i");
    _t_corr->Branch("ref_crate",&_ref_crate,"ref_crate/i");
    _t_corr->Branch("ref_slot",&_ref_slot,"ref_slot/i");
    _t_corr->Branch("ref_femch",&_ref_femch,"ref_femch/i");
    _t_corr->Branch("sub_larch",&_sub_larch,"sub_larch/i");
    _t_corr->Branch("sub_crate",&_sub_crate,"sub_crate/i");
    _t_corr->Branch("sub_slot",&_sub_slot,"sub_slot/i");
    _t_corr->Branch("sub_femch",&_sub_femch,"sub_femch/i");
    _t_corr->Branch("block",&_block,"block/i");
    _t_corr->Branch("rms_ref",&_rms_ref,"rms_ref/F");
    _t_corr->Branch("rms_sub",&_rms_sub,"rms_sub/F");
    _t_corr->Branch("rms_diff",&_rms_diff,"rms_diff/F");
    _t_corr->Branch("mean_ref",&_mean_ref,"mean_ref/F");
    _t_corr->Branch("mean_sub",&_mean_sub,"mean_sub/F");
    _t_corr->Branch("mean_diff",&_mean_diff,"mean_diff/F");
    _t_corr->Branch("corr",&_corr,"corr/F");
    _t_corr->Branch("corrall",&_corrall,"corrall/F");
    _t_corr->Branch("corrzig",&_corrzig,"corrzig/F");
    //_t_corr->Branch("corr2",&_corr2,"corr2/F");
    //_t_corr->Branch("diff_v","std::vector<short>",&_diff_v);
    _t_corr->SetDirectory(0);
    //_tree_v[counter_tree] = t;
    //counter_tree += 1;
    //    }

    return true;
  }
  
  bool ZigZag::analyze(storage_manager* storage) {

    auto wfs = storage->get_data<event_rawdigit>("daq");
    size_t nwfs = wfs->size();

    TStopwatch fWatch;
    fWatch.Start();
    
    std::map<int,int> chmap;

    // Prepare vectors
    _ch_to_index.resize(nwfs);;

    _zzMap.clear();
    
    // per-event quantities
    _frac = 0;
    _bursts.clear();
    _bursts.resize(9600);

    // Make a map from channel index to channel number
    for (size_t i=0; i < wfs->size(); i++){
      auto const& wf = (*wfs).at(i);
      int chnum = wf.Channel();
      chmap[chnum] = i;

      _larch = chnum;
      _crate = _ch_to_crate.at(chnum);
      _slot  = _ch_to_slot.at(chnum);
      _femch = _ch_to_femch.at(chnum);
      
      _adc_v.clear();
      
      auto const& adcs = wf.ADCs();

      
      // Mean
      for(size_t j=0; j<adcs.size(); ++j){
	_adc_v.push_back(adcs[j]);
	_mean += adcs[j];
      }
      _mean /= ((float)adcs.size());
      // RMS
      for(size_t k=0; k<adcs.size(); ++k)
	_rms += (adcs[k]-_mean)*(adcs[k]-_mean);
      _rms = sqrt( _rms / ((float)adcs.size()));
      //if (_ref_channel_s.find(chnum) == _ref_channel_s.end())
      //_adc_v.clear();
      
      // Find ZigZag region
      _zzMap[chnum] = findZigZag(adcs,_mean);
      if (_verbose) { std::cout << "found zig-zag for Chan: " << chnum << " with " << _zzMap[chnum].size() << " segments" << std::endl; }
      std::vector<std::pair<short,short> > thiszz = _zzMap[chnum];
      _zigzag_start.clear();
      _zigzag_end.clear();
      for (size_t y=0; y < thiszz.size(); y++){
	_zigzag_start.push_back(thiszz[y].first);
	_zigzag_end.push_back(thiszz[y].second);
	int num = thiszz[y].second - thiszz[y].first;
	_frac += num;
	float maxAmp = 0;
	for (int k =0; k < num; k++){
	  _bursts[thiszz[y].first+k] += 1;
	  if ((adcs[thiszz[y].first+k]-_mean) > maxAmp)
	    maxAmp = (adcs[thiszz[y].first+k]-_mean);
	}
	_amp = maxAmp;
	_start = thiszz[y].first;
	_end  = thiszz[y].second;
	_len = _end-_start;
	_t_z->Fill();
      }
      _t_ch->Fill();	
    }

    _frac /= (9600.*(wfs->size()));
    _t_zig->Fill();

    /*
    // Loop over blocks
    //size_t loops = wf.size()/((float)_nsamples);
    size_t loops = 1;//9500./_nsamples;
    for (size_t l=0; l < loops; l++){

      _block = l+(_evtN*loops);

      _ch_mean_v.clear();
      _ch_rms_v.clear();
      _ch_mean_v.resize(nwfs);
      _ch_rms_v.resize(nwfs);

      size_t ctr=0;

      // Loop over waveforms in list
      //for(auto const& ch1 : _ref_channel_s) {
      //auto const& wf = (*wfs).at(chmap[ch1]);
      for (size_t chn = 0; chn < wfs->size(); chn++){
	auto const& wf = (*wfs).at(chn);
	auto const& ch = wf.Channel();      

	auto const& adcs = wf.fADC;

	int idx = chmap[ch];

	_adc_v.clear();

	_larch = ch;

	_mean=0;
	_rms =0;
	
	_crate = _ch_to_crate.at(ch);
	_slot  = _ch_to_slot.at(ch);
	_femch = _ch_to_femch.at(ch);

	// Mean
	for(size_t j=0; j<_nsamples; ++j){
	  _adc_v.push_back(adcs[j+l*_nsamples]);
	  _mean += adcs[j+l*_nsamples];
	}
	_mean /= ((float)_nsamples);
	// RMS
	for(size_t k=0; k<_nsamples; ++k)
	  _rms += (adcs[k+l*_nsamples]-_mean)*(adcs[k+l*_nsamples]-_mean);
	_rms = sqrt( _rms / ((float)_nsamples));

	_ch_mean_v[idx] = _mean;
	_ch_rms_v[idx] = _rms;
	ctr++;
      }// for all channels
      
      std::cout<<"Finished computing baseline for "<<ctr<<" channels... "<<fWatch.RealTime()<<" [s]"<<std::endl;
      fWatch.Start();
      
      // For "reference" channels in list
      int counter_ch = 0;
      for(auto const& ch : _ref_channel_s) {
      auto const& wf_ref = (*wfs).at(chmap[ch]);
      //for (size_t chn = 0; chn < wfs->size(); chn++){
	//auto const& wf_ref = (*wfs).at(chn);
	auto const& ch_ref = wf_ref.Channel();
	_ref_larch = ch_ref;
	_ref_crate = _ch_to_crate.at(ch_ref);
	_ref_slot  = _ch_to_slot.at(ch_ref);
	_ref_femch = _ch_to_femch.at(ch_ref);
	
	_diff_v.clear();
	_mean_ref = _ch_mean_v.at(chmap[ch_ref]);
	_rms_ref  = _ch_rms_v.at(chmap[ch_ref]);
	
	fWatch.Start();
	//for(auto const& ch2 : _ref_channel_s) {
	//auto const& wf_sub = (*wfs).at(chmap[ch2]);

	for (size_t chn = 0; chn < wfs->size(); chn++){
	  auto const& wf_sub = (*wfs).at(chn);

	  auto const& ch_sub = wf_sub.Channel();
	  _sub_larch = ch_sub;
	  _sub_crate = _ch_to_crate.at(ch_sub);
	  _sub_slot  = _ch_to_slot.at(ch_sub);
	  _sub_femch = _ch_to_femch.at(ch_sub);
	  _mean_sub = _ch_mean_v.at(chmap[ch_sub]);
	  _rms_sub  = _ch_rms_v.at(chmap[ch_sub]);
	  
	  _mean_diff = 0;
	  _rms_diff  = 0;
	  auto const& adc_ref = wf_ref.fADC;
	  auto const& adc_sub = wf_sub.fADC;
	  for(size_t i=0; i<_nsamples; ++i) {
	    _diff_v[i] = adc_ref[i+l*_nsamples] - adc_sub[i+l*_nsamples];
	    _mean_diff += _diff_v[i+l*_nsamples];
	  }
	  _mean_diff /= ((float)_nsamples);
	  
	  for(size_t i=0; i<_nsamples; ++i)
	    _rms_diff += pow((float)(_diff_v[i+l*_nsamples]) - _mean_diff,2);
	  _rms_diff = sqrt( _rms_diff / ((float)_nsamples));
	  
	  _corr = 0;
	  _corrall = 0;
	  _corrzig = 0;
	  int counter = 0;
	  int zigsamples = 0;
	  // markers keeping track of the position in the vector of zig-zags that we have reached
	  size_t idx1 = 0;
	  size_t idx2 = 0;
	  size_t zig1 = 0;
	  size_t zig2 = 0;
	  for(size_t j=0; j<_nsamples; ++j){
	    // need to make sure this tick is in neither waveform's zig-zag region
	    //std::cout << "using: ch: [" << ch_ref << ", " << ch_sub << "]. Ticks: " 
	    //	      << j+l*_nsamples << " Idx: [" << idx1 << ", " << idx2 << "]" << std::endl;
	    double corrthistick = ((float)(adc_ref[j+l*_nsamples])-_mean_ref)*((float)(adc_sub[j+l*_nsamples])-_mean_sub);
	    if ( isNotZigZagRegion(ch_ref,ch_sub,j+l*_nsamples,idx1,idx2) ){
	      _corr += corrthistick;
	      counter += 1;
	    }
	    if ( isZigZagRegion(ch_ref,ch_sub,j+l*_nsamples,zig1,zig2) ){
	      _corrzig += corrthistick;
	      zigsamples += 1;
	    }	
	    _corrall += corrthistick;
	  }
	  if (_verbose) {std::cout << "Chans: [" << ch_ref << ", " << ch_sub << "]. Tick fraction: " << float(counter)/_nsamples << std::endl; }
	  if ( (_rms_sub != 0) and (_rms_ref != 0) ){
	    if (_verbose) { std::cout << "RMS ref: " << _rms_ref << "\tRMS sub: " << _rms_sub << "\t Ticks: " << _nsamples << std::endl; } 
	    if (counter != 0)
	      _corr /= (float)(_rms_sub*_rms_ref*counter);
	    else
	      _corr = -2;
	    _corrall /= (float)(_rms_sub*_rms_ref*_nsamples);
	    if (zigsamples != 0)
	      _corrzig /= (float)(_rms_sub*_rms_ref*zigsamples);
	    else
	      _corrzig = -2;
	  }
	  if (_verbose) {std::cout << "Corr: " << _corr << "\tCorr All: " << _corrall << "\tCorr Zig: " << _corrzig << std::endl << std::endl; }
	  
	  //_tree_v.at(counter_ch)->Fill();
	  _t_corr->Fill();
	}// for all "sub" waveforms
	counter_ch += 1;
      }// for all "ref" waveforms
    }// for all loops
    */
    _evtN += 1;
    
    return true;
  }

  bool ZigZag::finalize() {
    
    if(_fout){
      _fout->cd();
      //      for(auto& t : _tree_v)
      //	if(t) t->Write();
      std::cout << "writing corr tree" << std::endl;
      _t_corr->Write();
      std::cout << "writing ch tree" << std::endl;
      _t_ch->Write();
      std::cout << "writing zigzag tree" << std::endl;
      _t_zig->Write();
      std::cout << "writing zigzag tree" << std::endl;
      _t_z->Write();
      std::cout << "done writing trees" << std::endl;
    }
    return true;
  }

  bool ZigZag::isZigZag(float const& a, float const& b, float const& c){
    if ( (((a > 0) and (b < 0) and (c > 0)) or ((a < 0) and (b > 0) and (c < 0))) and
	 ( fabs(a-b) > 3 ) and ( fabs(b-c) > 3 ) ){
      return true;
    }
    else
    return false;
  }

  std::vector<std::pair<short,short> > ZigZag::findZigZag(std::vector<short> const& adcs, float const& baseline) {

    short start = 0;
    short end   = 0;
    size_t ticksOn = 0;
    std::vector<std::pair<short,short> > ROIs;
    std::pair<short,short> ROI;
    for (size_t i=0; i < 8300; i++){//adcs.size()-2; i++){
      // if 3 consecutive ticks show zig-zag
      if ( isZigZag(float(adcs[i])-baseline,float(adcs[i+1])-baseline,float(adcs[i+2])-baseline) ){
	// if not in ROI restart
	if (ticksOn == 0)
	  start = i;
	ticksOn += 1;
      }
      // if not in ROI (anymore)
      else{
	// if we were in ROI -> close segment
	if (ticksOn != 0){
	  // if ROI long enough -> make pair
	  if (ticksOn >= _minLength){
	    end = i+1;
	    ROI = std::make_pair(start,end);
	    ROIs.push_back(ROI);
	  }
	}
	ticksOn = 0;
      }
    }
    
    return ROIs;
  }
  
  bool ZigZag::isNotZigZagRegion(unsigned int wf1, unsigned int wf2, short tick, size_t &idx1, size_t &idx2){

    // boolean to keep track if in zigzag od one wf
    bool zz1 = false;
    bool zz2 = false;

    if (_zzMap[wf1].size() > idx1){
      auto const& reg1 = _zzMap[wf1][idx1];
      short start1 = reg1.first;
      short end1   = reg1.second;
      //if (_verbose) { std::cout << "Wf1: [" << start1 << ", " << end1 << "] " << tick << std::endl; }
      if (end1 < tick)
	idx1 += 1;
      if ( (tick > start1) and (tick < end1) )
	zz1 = true;
    }
    /*
    else
      if (_verbose) { std::cout << "Out of bounds for WF1. segments: " << _zzMap[wf1].size() << " Idx: " << idx1 << std::endl; }
    */
    if (_zzMap[wf2].size() > idx2){
      auto const& reg2 = _zzMap[wf2][idx2];
      short start2 = reg2.first;
      short end2   = reg2.second;
      //if (_verbose) { std::cout << "Wf2: [" << start2 << ", " << end2 << "] " << tick << std::endl; }
      if (end2 < tick)
	idx2 += 1;
      if ( (tick > start2) and (tick < end2) )
	zz2 = true;
    }
    /*
    else
      if (_verbose) { std::cout << "Out of bounds for WF1. segments: " << _zzMap[wf1].size() << " Idx: " << idx1 << std::endl; }
    */
    if ( (zz1 == false) and (zz2 == false) )
      return true;
    
    return false;
  }


  bool ZigZag::isZigZagRegion(unsigned int wf1, unsigned int wf2, short tick, size_t &idx1, size_t &idx2){

    // boolean to keep track if in zigzag od one wf
    bool zz1 = false;
    bool zz2 = false;

    if (_zzMap[wf1].size() > idx1){
      auto const& reg1 = _zzMap[wf1][idx1];
      short start1 = reg1.first;
      short end1   = reg1.second;
      //if (_verbose) { std::cout << "Wf1: [" << start1 << ", " << end1 << "] " << tick << std::endl; }
      if (end1 < tick)
	idx1 += 1;
      if ( (tick > start1) and (tick < end1) )
	zz1 = true;
    }
    /*
    else
      if (_verbose) { std::cout << "Out of bounds for WF1. segments: " << _zzMap[wf1].size() << " Idx: " << idx1 << std::endl; }
    */
    if (_zzMap[wf2].size() > idx2){
      auto const& reg2 = _zzMap[wf2][idx2];
      short start2 = reg2.first;
      short end2   = reg2.second;
      //if (_verbose) { std::cout << "Wf2: [" << start2 << ", " << end2 << "] " << tick << std::endl; }
      if (end2 < tick)
	idx2 += 1;
      if ( (tick > start2) and (tick < end2) )
	zz2 = true;
    }
    /*
    else
      if (_verbose) { std::cout << "Out of bounds for WF1. segments: " << _zzMap[wf1].size() << " Idx: " << idx1 << std::endl; }
    */
    if ( (zz1 == true) and (zz2 == true) )
      return true;
    
    return false;
  }

  
}
#endif
