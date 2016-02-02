/**
 * \file MuCSTagger.h
 *
 * \ingroup Paddles
 * 
 * \brief Class def header for a class MuCSTagger
 *
 * @author kazuhiro
 */

/** \addtogroup Paddles

    @{*/

#ifndef LARLITE_MUCSTAGGER_H
#define LARLITE_MUCSTAGGER_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include <TH1D.h>
#include <string>
namespace larlite {
  /**
     \class MuCSTagger
     User custom analysis class made by SHELL_USER_NAME
   */
  class MuCSTagger : public ana_base{
  
  public:

    /// Default constructor
    MuCSTagger();

    /// Default destructor
    virtual ~MuCSTagger(){}

    /** IMPLEMENT in MuCSTagger.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MuCSTagger.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MuCSTagger.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void configure(const std::string config_file);

    // Getter functions
    const std::string& producer() const { return _producer; }

    const std::vector<geoalgo::Trajectory> matched_trajectory() const { return _matched_trj_v; }

    const std::vector<geoalgo::HalfLine> matched_dir() const { return _matched_dir_v; }

    const ::geoalgo::AABox upper_box() const { return _mucs_upper_box; }

    const ::geoalgo::AABox lower_box() const { return _mucs_lower_box; }

    const std::vector<::geoalgo::Vector> upper_pt() const { return _upper_pt;       }

    const std::vector<::geoalgo::Vector> lower_pt() const { return _lower_pt;       }

    const float get_ctag_score()const{return _ctag_score;}

    std::vector<geoalgo::HalfLine> get_prj_start(){return _temps1;}
    std::vector<geoalgo::HalfLine> get_prj_end(){return _temps2;}
    
    
  protected:

    //bool Intersect(const TVector3& start, const TVector3& end);
    bool Intersect(const ::geoalgo::HalfLine_t&);
    bool IntersectDumb(const TVector3& start, const TVector3& end);

    bool _configured;
    ::geoalgo::AABox _mucs_upper_box;
    ::geoalgo::AABox _mucs_lower_box;
    std::vector<::geoalgo::Vector> _upper_pt;
    std::vector<::geoalgo::Vector> _lower_pt;
    ::geoalgo::AABox _tpc_av;
    ::geoalgo::HalfLine _mucs_dir;
    
    std::vector<geoalgo::HalfLine> _matched_dir_v;
    std::vector<geoalgo::Trajectory> _matched_trj_v;

    std::vector<geoalgo::HalfLine> _temps1;
    std::vector<geoalgo::HalfLine> _temps2;
    
    std::string _producer;
    double _xmin;
    double _xmax;

    bool   _store_match;
    bool   _hit_both_box;
    double _min_track_length;
    double _segment_length;
    double _scan_length;
    bool   _allow_flip_direction;
    TH1D*  _hNumTracks;
    TH1D*  _hNumTagged;

    float _ctag_score;
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
