/**
 * \file LinearHitRemoval.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class LinearHitRemoval
 *
 * @author david caratelli
 */

/** \addtogroup Clusterer

    @{*/

#ifndef LARLITE_LINEARHITREMOVAL_H
#define LARLITE_LINEARHITREMOVAL_H

#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"
#include <map>

namespace larlite {
  /**
     \class LinearHitRemoval
     User custom analysis class made by SHELL_USER_NAME
   */
  class LinearHitRemoval : public ana_base{
  
  public:

    /// Default constructor
    LinearHitRemoval();

    /// Default destructor
    virtual ~LinearHitRemoval(){}

    /** IMPLEMENT in LinearHitRemoval.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in LinearHitRemoval.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in LinearHitRemoval.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    /// Set the size of each cell for hit-map
    void setCellSize(double d) { _cellSize = d; }
    /// Set the radius around which to search for hits
    /// if two hits are within this distance of each other
    /// then they go into the same cluster
    void setRadius(double d) { _radius = d; }
    /// set maximum linearity allowed for this
    void setMaxLinearity(double l) { _max_lin = l; }
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

    /// maximum linearity for hits
    double _max_lin;
    
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

    /// Functions to decide if two hits should belong to the same cluster or not
    bool HitsCompatible(const hit& h1, const hit& h2);

    /// Function to get neighboring hits (from self + neighoring cells)
    void getNeighboringHits(const std::pair<int,int>& pair, std::vector<size_t>& hitIndices);

    /// map connecting coordinate index (i,j) to [h1,h2,h3] (hit index list)
    std::map<std::pair<int,int>, std::vector<size_t> > _hitMap;

    /// covariance, standard deviation, mean
    double cov (const std::vector<double>& data1,
		const std::vector<double>& data2) const;
    double stdev(const std::vector<double>& data) const;
    double mean (const std::vector<double>& data) const;
    
    /// maximum i'th and j'th
    int _maxI;
    int _maxJ;

    
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
