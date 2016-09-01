/**
 * \file LinearClusterLocalRemoval.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class LinearClusterLocalRemoval
 *
 * @author david caratelli
 */

/** \addtogroup Clusterer

    @{*/

#ifndef LARLITE_LINEARCLUSTERLOCALREMOVAL_H
#define LARLITE_LINEARCLUSTERLOCALREMOVAL_H

#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"
#include <map>
#include "TTree.h"

namespace larlite {
  /**
     \class LinearClusterLocalRemoval
     User custom analysis class made by SHELL_USER_NAME
   */
  class LinearClusterLocalRemoval : public ana_base{
  
  public:

    /// Default constructor
    LinearClusterLocalRemoval();

    /// Default destructor
    virtual ~LinearClusterLocalRemoval(){}

    /** IMPLEMENT in LinearClusterLocalRemoval.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in LinearClusterLocalRemoval.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in LinearClusterLocalRemoval.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    /// set maximum linearity allowed for this
    void setMaxLinearity(double l) { _max_lin = l; }
    /// Verbosity setter
    void setVerbose(bool on) { _verbose = on; }
    /// Set Hit Producer
    void setClusterProducer(std::string s) { _clusterProducer = s; }
    /// set max radius
    void setMaxRadius(double r) { _radius = r; }
    /// set minimum number of hits in cluster
    void setMinNHits(int n) { _min_n_hits = n; }

  protected:

    /// maximum linearity for hits
    double _max_lin;

    /// minimum number of his
    int _min_n_hits;

    /// verbosity flag
    bool _verbose;

    /// conversion factors for hits
    double _wire2cm, _time2cm;

    /// max radius for local hit cluster
    double _radius;

    /// Hit producer name
    std::string _clusterProducer;

    /// covariance, standard deviation, mean
    double cov (const std::vector<double>& data1,
		const std::vector<double>& data2) const;
    double stdev(const std::vector<double>& data) const;
    double mean (const std::vector<double>& data) const;
    double linearity(const std::vector<double>& data1,
		     const std::vector<double>& data2) const;

    /// get local neighbordhood of hits
    void getNeighboringHits(const unsigned int& hit_idx,
			    const std::vector<unsigned int>& hit_idx_v,
			    std::vector<unsigned int>& out_hit_v,
			    const event_hit* ev_hit) const;

    TTree *_tree;
    double _l;
    int    _n_hits;
    int    _pl;
    
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
