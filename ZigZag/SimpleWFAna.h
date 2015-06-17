/**
 * \file SimpleWFAna.h
 *
 * \ingroup SimpleWFAna
 * 
 * \brief Class def header for a class SimpleWFAna
 *
 * @author davidc1
 */

/** \addtogroup SimpleWFAna

    @{*/

#ifndef SIMPLEWFANA_H
#define SIMPLEWFANA_H

#include "Analysis/ana_base.h"
#include <map>
#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include "DataFormat/rawdigit.h"

namespace larlite {
  /**
     \class SimpleWFAna
     User custom analysis class made by kterao
   */
  class SimpleWFAna : public ana_base{
  
  public:

    /// Default constructor
    SimpleWFAna(){ _name="SimpleWFAna"; _fout=0;  _verbose=false; };

    /// Default destructor
    virtual ~SimpleWFAna(){};

    /** IMPLEMENT in SimpleWFAna.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SimpleWFAna.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SimpleWFAna.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void setVerbose(bool on) { _verbose = on; }

  protected:

    bool _verbose;

    size_t _nsamples;    

    float _rms, _mean;

    unsigned int _larch;

    TTree* _t_ch;
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
