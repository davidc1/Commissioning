/**
 * \file SimpleClusterer.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class SimpleClusterer
 *
 * @author david caratelli
 */

/** \addtogroup Clusterer

    @{*/

#ifndef LARLITE_SIMPLECLUSTERER_H
#define LARLITE_SIMPLECLUSTERER_H

#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"
#include <map>

namespace larlite {
  /**
     \class SimpleClusterer
     User custom analysis class made by SHELL_USER_NAME
   */
  class SimpleClusterer : public ana_base{
  
  public:

    /// Default constructor
    SimpleClusterer();

    /// Default destructor
    virtual ~SimpleClusterer(){}

    /** IMPLEMENT in SimpleClusterer.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SimpleClusterer.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SimpleClusterer.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    /// Set the size of each cell for hit-map
    void setCellSize(double d) { _cellSize = d; }
    /// Set the radius around which to search for hits
    /// if two hits are within this distance of each other
    /// then they go into the same cluster
    void setRadius(double d) { _radius = d; }
    /// Set which plane to select hits from
    void setPlane(int pl) { _plane = pl; }
    /// Verbosity setter
    void setVerbose(bool on) { _verbose = on; }
    /// Set Hit Producer
    void setHitProducer(std::string s) { _hitProducer = s; }

  protected:

    /// size of each cell [cm]
    double _cellSize;

    /// radius to count charge around [cm]
    double _radius;
    
    /// plane to select hits from
    int _plane;

    /// verbosity flag
    bool _verbose;

    /// conversion factors for hits
    double _wire2cm, _time2cm;

    /// Hit producer name
    std::string _hitProducer;

    /// Map making function
    void MakeHitMap(const event_hit* hitlist, int plane);

    /// Function to get neighboring hits (from self + neighoring cells)
    std::vector<size_t> getNeighboringHits(const std::pair<int,int> pair);

    /// map connecting coordinate index (i,j) to [h1,h2,h3] (hit index list)
    std::map<std::pair<int,int>, std::vector<size_t> > _hitMap;


    
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
