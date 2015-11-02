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

namespace larlite {
  /**
     \class PaddleTrackFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class PaddleTrackFilter : public ana_base{
  
  public:

    /// Default constructor
    PaddleTrackFilter(){ _name="PaddleTrackFilter"; _fout=0;}

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

  protected:
    
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
