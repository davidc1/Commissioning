/**
 * \file ZigZag.h
 *
 * \ingroup SimpleWFAna
 * 
 * \brief Class def header for a class ZigZag
 *
 * @author davidc1
 */

/** \addtogroup SimpleWFAna

    @{*/

#ifndef ZIGZAG_H
#define ZIGZAG_H

#include "Analysis/ana_base.h"
#include <map>
#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include "DataFormat/rawdigit.h"

namespace larlite {
  /**
     \class ZigZag
     User custom analysis class made by kterao
   */
  class ZigZag : public ana_base{
  
  public:

    /// Default constructor
    ZigZag(){ _name="ZigZag"; _fout=0; _nsamples=9550; _minLength=4; _verbose=false; };

    /// Default destructor
    virtual ~ZigZag(){};

    /** IMPLEMENT in ZigZag.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ZigZag.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ZigZag.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void add(unsigned int n) { _ref_channel_s.insert(n); }

    void setNSamples(size_t n) { _nsamples = n; }

    void SetMap(unsigned int lar_ch,
		short crate,
		short slot,
		short femch);

    void setZigZagMinLength(size_t l) {_minLength = l; }

    void setVerbose(bool on) { _verbose = on; }

  protected:

    bool _verbose;

    bool isZigZag(float const& a, float const& b, float const& c);
    std::vector<std::pair<short,short> > findZigZag(std::vector<short> const& adcs, float const& baseline);
    bool isNotZigZagRegion(unsigned int wf1, unsigned int wf2, short tick, size_t &idx1, size_t &idx2);
    bool isZigZagRegion(unsigned int wf1, unsigned int wf2, short tick, size_t &idx1, size_t &idx2);

    // nimumum number of ticks to create ZigZag region
    size_t _minLength;

    // map linking channel number to zigzag regions
    std::map<unsigned int, std::vector<std::pair<short,short> > > _zzMap;

    // vector holding number of zig-zags in wf
    std::vector<short> _nZigZags;

    size_t _nsamples;    
    std::set<unsigned int> _ref_channel_s;
    std::vector<float>  _ch_mean_v;
    std::vector<float>  _ch_rms_v;

    int   _block;
    float _rms, _mean;
    std::vector<short> _zigzag_start; 
    std::vector<short> _zigzag_end; 
    float _rms_ref;
    float _rms_sub;
    float _rms_diff;
    float _mean_ref;
    float _mean_sub;
    float _mean_diff;
    float _corr, _corr2,_corrall, _corrzig;
    std::vector<short> _diff_v;
    std::vector<short> _adc_v;

    std::vector<int> _bursts;
    float _frac;

    unsigned int _larch, _crate, _slot, _femch;
    unsigned int _ref_larch, _ref_crate, _ref_slot, _ref_femch;
    unsigned int _sub_larch, _sub_crate, _sub_slot, _sub_femch;

    std::vector<int> _ch_to_index;
    std::vector<short> _ch_to_crate;
    std::vector<short> _ch_to_slot;
    std::vector<short> _ch_to_femch;

    std::vector<TTree*> _tree_v;
    TTree* _t_corr;
    TTree* _t_ch;
    TTree* _t_zig;
    TTree* _t_z;
    float _amp;
    int _start;
    int _end;
    int _len;

    int _evtN;

  };
}



#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
