/**
 * \file MakeHits.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class MakeHits
 *
 * @author david
 */

/** \addtogroup Clusterer

    @{*/

#ifndef LARLITE_MAKEHITS_H
#define LARLITE_MAKEHITS_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class MakeHits
     User custom analysis class made by SHELL_USER_NAME
   */
  class MakeHits : public ana_base{
  
  public:

    /// Default constructor
    MakeHits(){ _name="MakeHits"; _fout=0;}

    /// Default destructor
    virtual ~MakeHits(){}

    /** IMPLEMENT in MakeHits.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MakeHits.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MakeHits.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void setHitProducer(std::string p) { _producer = p; }

  protected:

    std::string _producer;
    
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
