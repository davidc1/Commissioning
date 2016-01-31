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

    const std::string& producer() const { return _producer; }

  protected:

    bool Intersect(const TVector3& start, const TVector3& end);

    bool _configured;
    ::geoalgo::AABox _mucs_upper_box;
    ::geoalgo::AABox _mucs_lower_box;
    ::geoalgo::HalfLine _mucs_dir;

    std::string _producer;
    double _xmin;
    double _xmax;
    
    bool   _hit_upper_box;
    bool   _hit_lower_box;
    double _min_track_length;
    double _segment_length;
    double _scan_length;
    bool   _allow_flip_direction;
    TH1D*  _hNumTracks;
    TH1D*  _hNumTagged;
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
