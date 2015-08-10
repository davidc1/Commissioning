/**
 * \file Overlapevents.h
 *
 * \ingroup OverlapEvents
 * 
 * \brief Class def header for a class OverlapEvents
 *
 * @author davidc1
 */

/** \addtogroup OverlapEvents

    @{*/

#ifndef LARLITE_OVERLAPEVENTS_H
#define LARLITE_OVERLAPEVENTS_H

#include "Analysis/ana_base.h"
#include "DataFormat/rawdigit.h"

namespace larlite {
  /**
     \class OverlapEvents
     User custom analysis class made by SHELL_USER_NAME
   */
  class OverlapEvents : public ana_base{
  
  public:

    /// Default constructor
    OverlapEvents(){ _name="OverlapEvents"; _fout=0;}

    /// Default destructor
    virtual ~OverlapEvents(){}

    /** IMPLEMENT in OverlapEvents.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in OverlapEvents.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in OverlapEvents.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    std::vector<std::vector<float> > _laser_chans;

    size_t _evt;
    
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
