/**
 * \file CosmicDiscrimFIFO.h
 *
 * \ingroup PMTCommissioning
 * 
 * \brief Class def header for a class CosmicDiscrimFIFO
 *
 * @author david caratelli
 */

/** \addtogroup PMTCommissioning

    @{*/

#ifndef LARLITE_COSMICDISCRIMFIFO_H
#define LARLITE_COSMICDISCRIMFIFO_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {

  
  static const double kINVALID_DOUBLE = std::numeric_limits<double>::max();
  
  static const double kMAX_DOUBLE = std::numeric_limits<double>::max();
  
  static const double kMIN_DOUBLE = std::numeric_limits<double>::min();

  /**
     \class CosmicDiscrimFIFO
     User custom analysis class made by SHELL_USER_NAME
   */
  class CosmicDiscrimFIFO : public ana_base{
  
  public:

    /// Default constructor
    CosmicDiscrimFIFO();

    /// Default destructor
    virtual ~CosmicDiscrimFIFO(){}

    /** IMPLEMENT in CosmicDiscrimFIFO.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CosmicDiscrimFIFO.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in CosmicDiscrimFIFO.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    // setter for whether to use trigger info
    void UseTrigger(bool on) { _use_trig = on; }
    // setter for whether to save waveforms or not
    void SaveWF(bool on) { _save_wf = on; }
    // trigger producer
    void SetTrigProducer(std::string s) { _trig_producer = s; }
    // fifo producer
    void SetFifoProducer(std::string s) { _fifo_producer = s; }
    // set ADC threshold (wf counted only if within the 20-sample window the ADC goes above that range)
    void SetADCThresh(short adc) { _adc_thresh = adc; }

    // verbosity
    void SetVerbose(bool on) { _verbose = on; }

  protected:
    
    // bool verbose
    bool _verbose;

    // boolean for whether to use trigger info
    bool _use_trig;
    // boolean : save waveforms?
    bool _save_wf;
    // trigger producer name
    std::string _trig_producer;
    // fifo producer name
    std::string _fifo_producer;
    // threshold for ADC
    short _adc_thresh;
    // baseline vector -> holds a baseline for each PMT
    std::vector<double> _baselines;
    // rms vector -> holds best RMS found per PMT
    std::vector<double> _rms;

    // function to get baseline and rms
    std::pair<double,double> GetBaselineRms(const std::vector<unsigned short>& wf);

    TTree* _trig_tree;
    int _trig_num;
    double _trig_time;
    double _delta_t;
    double _last_trig_time;
    int _last_frame_num;

    TTree* _tree;
    int _ch;
    int _ev;
    int _adcs;
    int _amp;
    int _first; // ADC-baseline for first tick value in 20-sample window
    int _disc;
    int _frame;
    int _sample;
    int _ev_frame;
    int _frame_diff; // wf frame # - event frame #
    double _avg_ch;
    double _rms_ch;
    std::vector<unsigned short> _adc_v;

    TTree* _rate_tree;
    std::vector<unsigned short> _n20_windows_l;
    std::vector<unsigned short> _n20_windows_h;
    std::vector<unsigned short> _n1k_windows;
    std::vector<double> _rms_v;
    std::vector<double> _baseline_v;
    
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
