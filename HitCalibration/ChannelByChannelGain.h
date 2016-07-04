/**
 * \file ChannelByChannelGain.h
 *
 * \ingroup HitCalibration
 * 
 * \brief Class def header for a class ChannelByChannelGain
 *
 * @author davidc1
 */

/** \addtogroup HitCalibration

    @{*/

#ifndef LARLITE_CHANNELBYCHANNELGAIN_H
#define LARLITE_CHANNELBYCHANNELGAIN_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class ChannelByChannelGain
     User custom analysis class made by SHELL_USER_NAME
   */
  class ChannelByChannelGain : public ana_base{
  
  public:

    /// Default constructor
  ChannelByChannelGain() : _tree(nullptr) { _name="ChannelByChannelGain"; _fout=0;}

    /// Default destructor
    virtual ~ChannelByChannelGain(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();
    
    void SetHitProducer(std::string s) { _hit_producer = s; }
    
  protected:
    
    std::string _hit_producer;
    
    TTree *_tree;
    int    _pl;
    double _area;
    double _chi;
    int    _ch;
    int    _run;
    int    _subrun;
    int    _event;
    
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
