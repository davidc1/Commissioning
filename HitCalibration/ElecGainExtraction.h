/**
 * \file ElecGainExtraction.h
 %*
 * \ingroup HitCalibration
 * 
 * \brief Class def header for a class ElecGainExtraction
 *
 * @author david caratelli
 */

/** \addtogroup HitCalibration

    @{*/

#ifndef LARLITE_ELECGAINEXTRACTION_H
#define LARLITE_ELECGAINEXTRACTION_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class ElecGainExtraction
     User custom analysis class made by SHELL_USER_NAME
   */
  class ElecGainExtraction : public ana_base{
  
  public:

    /// Default constructor
    ElecGainExtraction(){ _name="ElecGainExtraction"; _fout=0;}

    /// Default destructor
    virtual ~ElecGainExtraction(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    /// set hit producer
    void setHitProducer(std::string s) { _hit_producer = s; }
    void setClusProducer(std::string s) { _clus_producer = s; }

  protected:

    /// hit producer
    std::string _hit_producer;
    std::string _clus_producer;

    /// TTree
    TTree* _tree;
    double _reco_area, _reco_ampl;
    double _q; // charge form simch
    int _tick;
    int _pl;
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
