/**
 * \file ValidateNIM.h
 *
 * \ingroup Paddles
 * 
 * \brief Class def header for a class ValidateNIM
 *
 * @author david
 */

/** \addtogroup Paddles

    @{*/

#ifndef LARLITE_VALIDATENIM_H
#define LARLITE_VALIDATENIM_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class ValidateNIM
     User custom analysis class made by SHELL_USER_NAME
   */
  class ValidateNIM : public ana_base{
  
  public:

    /// Default constructor
    ValidateNIM();

    /// Default destructor
    virtual ~ValidateNIM(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setVerbose(bool on) { _verbose = on; }

  protected:

    void resetTriggers();

    TTree *_tree;
    int _trig_num;
    double _trig_time;
    double _beam_time;
    int _trig_bits;

    int _rwm;
    int _strb1;
    int _strb2;
    int _numi;
    int _bnb;
    int _led;
    int _flash;
    int _pdls;

    bool _verbose;
    
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
