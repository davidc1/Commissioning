/**
 * \file FillTree.h
 *
 * \ingroup BeamTiming
 * 
 * \brief Class def header for a class FillTree
 *
 * @author david
 */

/** \addtogroup BeamTiming

    @{*/

#ifndef LARLITE_FILLTREE_H
#define LARLITE_FILLTREE_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class FillTree
     User custom analysis class made by SHELL_USER_NAME
   */
  class FillTree : public ana_base{
  
  public:

    /// Default constructor
    FillTree();

    /// Default destructor
    virtual ~FillTree(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setFlashProducer(std::string s) { _flashProducer = s; }
    void setMinPE(double pe) { _pe_min = pe; }
    void setMinT(double t) { _minT = t; }
    void setMaxT(double t) { _maxT = t; }

  protected:

    TTree *_tree;

    int _trig_word;
    int _run, _subrun, _event;
    double _pe_total;
    double _dt;
    double _flash_time, _flash_abstime;
    double _trig_time;

    std::string _flashProducer;
    double _pe_min;
    double _minT, _maxT;
    
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
