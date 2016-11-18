/**
 * \file PhotonClusterer.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class PhotonClusterer
 *
 * @author david
 */

/** \addtogroup Clusterer

    @{*/

#ifndef LARLITE_PHOTONCLUSTERER_H
#define LARLITE_PHOTONCLUSTERER_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class PhotonClusterer
     User custom analysis class made by SHELL_USER_NAME
   */
  class PhotonClusterer : public ana_base{
  
  public:

    /// Default constructor
    PhotonClusterer();

    /// Default destructor
    virtual ~PhotonClusterer(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    /// Set Producers
    void setClusterProducer(std::string s) { _clusterProducer = s; }

    /// set maximum linearity allowed for this
    void setMaxLinearity(double l) { _max_lin_v.push_back( l ); }
    
    /// Verbosity setter
    void setVerbose(bool on) { _verbose = on; }

    /// set minimum number of hits
    void setMinNHits(int n) { _min_n_hits_v.push_back( n ); }

    void setNMin(int   n) { _minHits = n; }
    void setNMax(int   n) { _maxHits = n; }
    void setQMin(float q) { _minQ    = q; }
    
  protected:
    
    bool _verbose;
    std::string _clusterProducer;
    int _maxHits, _minHits;
    float _minQ;
    std::vector<int> _min_n_hits_v;
    std::vector<float> _max_lin_v;

    double _wire2cm, _time2cm;
    
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
