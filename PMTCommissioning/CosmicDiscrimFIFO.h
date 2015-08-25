/**
 * \file CosmicDiscrimFIFO.h
 *
 * \ingroup PMTCommissioning
 * 
 * \brief Class def header for a class CosmicDiscrimFIFO
 *
 * @author davidc1
 */

/** \addtogroup PMTCommissioning

    @{*/

#ifndef LARLITE_COSMICDISCRIMFIFO_H
#define LARLITE_COSMICDISCRIMFIFO_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class CosmicDiscrimFIFO
     User custom analysis class made by SHELL_USER_NAME
   */
  class CosmicDiscrimFIFO : public ana_base{
  
  public:

    /// Default constructor
    CosmicDiscrimFIFO()
      : _tree(nullptr)
      , _trig_tree(nullptr)
    { _name="CosmicDiscrimFIFO"; _fout=0;}

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

  protected:

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
    int _disc;
    int _frame;
    int _sample;
    int _ev_frame;
    int _frame_diff; // wf frame # - event frame #
    std::vector<unsigned short> _adc_v;
    
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
