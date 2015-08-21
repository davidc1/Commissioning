/**
 * \file CosmicDiscrim.h
 *
 * \ingroup PMTCommissioning
 * 
 * \brief Class def header for a class CosmicDiscrim
 *
 * @author david
 */

/** \addtogroup PMTCommissioning

    @{*/

#ifndef LARLITE_COSMICDISCRIM_H
#define LARLITE_COSMICDISCRIM_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class CosmicDiscrim
     User custom analysis class made by SHELL_USER_NAME
   */
  class CosmicDiscrim : public ana_base{
  
  public:

    /// Default constructor
    CosmicDiscrim()
      : _tree(nullptr)
    { _name="CosmicDiscrim"; _fout=0;}

    /// Default destructor
    virtual ~CosmicDiscrim(){}

    /** IMPLEMENT in CosmicDiscrim.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CosmicDiscrim.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:

    TTree* _tree;
    int _ch;
    int _ev;
    int _adcs;
    double _tstart;
    double _tend;
    std::vector<short> _adc_v;
    
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
