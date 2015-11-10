/**
 * \file PaddleTrackOpflash.h
 *
 * \ingroup Paddles
 * 
 * \brief Class def header for a class PaddleTrackOpflash
 *
 * @author Rui
 */

/** \addtogroup Paddles

    @{*/

#ifndef LARLITE_PADDLETRACKOPFLASH_H
#define LARLITE_PADDLETRACKOPFLASH_H

#include "Analysis/ana_base.h"
#include "DataFormat/track.h"
#include "DataFormat/ophit.h"
#include "DataFormat/opflash.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/Geometry.h"
#include "TTree.h"
#include "PaddleTrackAna.h"

namespace larlite {
  /**
     \class PaddleTrackOpflash
     User custom analysis class made by SHELL_USER_NAME
   */
  class PaddleTrackOpflash : public ana_base{
  
  public:

    /// Default constructor
    PaddleTrackOpflash()
      :_tree(nullptr)
      { _name="PaddleTrackOpflash";
	_fout=0;
	_fout=0;
      }

    /// Default destructor
    virtual ~PaddleTrackOpflash(){}

    /** IMPLEMENT in PaddleTrackOpflash.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in PaddleTrackOpflash.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in PaddleTrackOpflash.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:
    
    TTree* _tree;
    TH2F* _x_z_ints_top;
    TH2F* _x_z_ints_bottom;
    bool _save_histos;

    size_t _n_ev_reco;

    int _n_evt;
    int _n_tracks;
    int _n_intersections_FV;
    int _n_intersections_mucs_top;
    int _n_intersections_mucs_bottom;
    
    std::vector<double> _t_opflash;
    
    //muon intersection w/ MuCS
    double _MuCS_ints_x_top;
    double _MuCS_ints_z_top;
    double _MuCS_ints_x_bottom;
    double _MuCS_ints_z_bottom;
    
    double _length_xfiducial;
    double _length_yfiducial;
    double _length_zfiducial;
    
    ::geoalgo::GeoAlgo _geoAlgo;
    ::geoalgo::AABox _vfiducial;
    ::geoalgo::AABox _vmucs_top;
    ::geoalgo::AABox _vmucs_bottom;
    std::vector<::geoalgo::Trajectory> _trj;
    //tracks contained in TPC
    std::vector<::geoalgo::Trajectory> _trj_con;
    //tracks pass MuCS
    std::vector<::geoalgo::Trajectory> _trj_mucs;
    std::vector<::geoalgo::HalfLine> _trj_prj;
    std::vector<::geoalgo::LineSegment> _prj_lineseg;
    
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
