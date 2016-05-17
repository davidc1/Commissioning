/**
 * \file TrkMIPCalibration.h
 *
 * \ingroup HitCalibration
 * 
 * \brief Class def header for a class TrkMIPCalibration
 *
 * @author david
 */

/** \addtogroup HitCalibration

    @{*/

#ifndef LARLITE_TRKMIPCALIBRATION_H
#define LARLITE_TRKMIPCALIBRATION_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoTrajectory.h"
#include "AnalysisAlg/CalorimetryAlg.h"
#include "TTree.h"

namespace larlite {
  /**
     \class TrkMIPCalibration
     User custom analysis class made by SHELL_USER_NAME
   */
  class TrkMIPCalibration : public ana_base{
  
  public:

    /// Default constructor
    TrkMIPCalibration();

    /// Default destructor
    virtual ~TrkMIPCalibration(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetTrkProducer(std::string s) { _trk_producer = s; }

  protected:

    double _recomb_factor;

    std::string _trk_producer;

    geoalgo::Vector_t Get3DDir(const geoalgo::Trajectory_t& trj);

    // calorimetry service
    calo::CalorimetryAlg _caloAlg;

    // output tree
    TTree *_tree;
    double _trk_len;
    int _plane;
    double _pitch;
    std::vector<double> _amp_v;
    std::vector<double> _dqdx_v;
    std::vector<double> _dedx_v;
    
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
