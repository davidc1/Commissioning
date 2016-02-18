/**
 * \file SignalProcAreaNormalization.h
 *
 * \ingroup HitCalibration
 * 
 * \brief Class def header for a class SignalProcAreaNormalization
 *
 * @author david caratelli
 */

/** \addtogroup HitCalibration

    @{*/

#ifndef LARLITE_SIGNALPROCAREANORMALIZATION_H
#define LARLITE_SIGNALPROCAREANORMALIZATION_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class SignalProcAreaNormalization
     User custom analysis class made by SHELL_USER_NAME
   */
  class SignalProcAreaNormalization : public ana_base{
  
  public:

    /// Default constructor
    SignalProcAreaNormalization(){ _name="SignalProcAreaNormalization"; _fout=0;}

    /// Default destructor
    virtual ~SignalProcAreaNormalization(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    /// set hit producer
    void setHitProducer(std::string s) { _hit_producer = s; }

  protected:

    /// hit producer
    std::string _hit_producer;

    /// TTree
    TTree* _tree;
    double _reco_area, _raw_area;
    int _tick;
    int _trk_size;
    int _wire, _chan;
    int _hit_multiplicity;
    
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
