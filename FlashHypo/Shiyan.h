/**
 * \file Shiyan.h
 *
 * \ingroup FlashHypo
 * 
 * \brief Class def header for a class Shiyan
 *
 * @author Rui
 */

/** \addtogroup FlashHypo

    @{*/

#ifndef LARLITE_SHIYAN_H
#define LARLITE_SHIYAN_H

#include "Analysis/ana_base.h"
#include <iostream>
#include "OpT0Finder/PhotonLibrary/PhotonVisibilityService.h"
#include <numeric>
#include "Analysis/ana_base.h"
#include "DataFormat/track.h"
#include "DataFormat/ophit.h"
#include "DataFormat/opflash.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/Geometry.h"
#include <functional>
#include <algorithm>
#include "TTree.h"

namespace larlite {
  /**
     \class Shiyan
     User custom analysis class made by SHELL_USER_NAME
   */
  class Shiyan : public ana_base{
  
  public:

    /// Default constructor
    Shiyan()
      :_tree(nullptr)
      { _name="Shiyan"; _fout=0;}
    
        
    /// Default destructor
    virtual ~Shiyan(){}

    /** IMPLEMENT in Shiyan.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in Shiyan.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in Shiyan.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    std::vector<double> getPE(){return _pe;}

  protected:
    std::vector<double> _pe;
    TTree* _tree;
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
