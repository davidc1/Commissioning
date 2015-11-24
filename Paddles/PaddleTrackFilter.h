/**
 * \file PaddleTrackFilter.h
 *
 * \ingroup Paddles
 * 
 * \brief Class def header for a class PaddleTrackFilter
 *
 * @author Rui
 */

/** \addtogroup Paddles

    @{*/

#ifndef LARLITE_PADDLETRACKFILTER_H
#define LARLITE_PADDLETRACKFILTER_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/Geometry.h"
#include "TTree.h"


namespace larlite {
  /**
     \class PaddleTrackFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class PaddleTrackFilter : public ana_base{
  
  public:

    /// Default constructor
    PaddleTrackFilter()
      : _tree(nullptr)
      { _name="PaddleTrackFilter"; _fout=0;}

    /// Default destructor
    virtual ~PaddleTrackFilter(){}

    /** IMPLEMENT in PaddleTrackFilter.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in PaddleTrackFilter.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in PaddleTrackFilter.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    // getter function for TPC volume
    ::geoalgo::AABox getFiducialVol() { return _vfiducial; }
    // getter function for top MCS
    ::geoalgo::AABox getTopMuCS() { return _vmucs_top; }
    // getter function for bottom MCS
    ::geoalgo::AABox getBottomMuCS() { return _vmucs_bottom; }
    

    // getter function for tracks in an event
    std::vector<::geoalgo::Trajectory> getTrj()      {return _trj;}
    // getter function for tracks in an event contained by TPCFV
    std::vector<::geoalgo::Trajectory> getTrj_con()  {return _trj_con;}
    // getter function for tracks in an event pass MuCS
    std::vector<::geoalgo::Trajectory> getTrj_mucs() {return _trj_mucs;}
    // getter function for back-projected half line
    std::vector<::geoalgo::HalfLine> getHalfLine()   {return _trj_prj;}
    // getter function for linesegment (proj w/ first and end point)
    std::vector<::geoalgo::LineSegment> getLineseg() {return _prj_lineseg;}
    // getter functions for event information
    int getRun()    {return _run;}
    int getSubrun() {return _subrun;}
    int getEvent()  {return _event;}
    int getTrackID(){return _trk_id;}
    
  protected:
    
    TTree* _tree;
    
    size_t _n_ev_reco;

    int _run;
    int _subrun;
    int _event;
    int _trk_id;

    
    int _n_evt;
    int _n_tracks;
    int _n_intersections_FV;
    int _n_intersections_mucs_top;
    int _n_intersections_mucs_bottom;
    
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
