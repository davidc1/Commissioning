/**
 * \file LaserFinder.h
 *
 * \ingroup LaserFinder
 * 
 * \brief Class def header for a class LaserFinder
 *
 * @author cadams
 */

/** \addtogroup LaserFinder

    @{*/

#ifndef LARLITE_LASERFINDER_H
#define LARLITE_LASERFINDER_H

#include "Analysis/ana_base.h"
#include "TVector3.h"

namespace larlite {
  /**
     \class LaserFinder
     User custom analysis class made by SHELL_USER_NAME
   */
  class LaserFinder : public ana_base{
  
  public:

    /// Default constructor
    LaserFinder(){ _name="LaserFinder"; _fout=0;}

    /// Default destructor
    virtual ~LaserFinder(){}

    /** IMPLEMENT in LaserFinder.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in LaserFinder.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in LaserFinder.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();


    float DistanceToLine3D( const TVector3 & pointOnLine, const TVector3 & directionOfLine, 
                            const TVector3 & targetPoint) const;

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
