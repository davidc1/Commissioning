/**
 * \file ClusterFilter.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class ClusterFilter
 *
 * @author david
 */

/** \addtogroup Clusterer

    @{*/

#ifndef LARLITE_CLUSTERFILTER_H
#define LARLITE_CLUSTERFILTER_H

#include "Analysis/ana_base.h"

#include "Linearity.h"

#include "Linearity.h"

namespace larlite {
  /**
     \class ClusterFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class ClusterFilter : public ana_base{
  
  public:

    /// Default constructor
    ClusterFilter();

    /// Default destructor
    virtual ~ClusterFilter(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setMaxNHits(int n)    { _max_n_hits = n; }
    void setMaxArea (double a) { _Amax = a;       }
    void setMaxDist (double d) { _d_max = d;      }

    void setClusProducer(std::string s) { _clusProducer = s; }
    void setVtxProducer (std::string s) { _vtxProducer  = s; }
    void setOutHitProducer(std::string s) { _outhitProducer = s; }
    void setOutClusterProducer(std::string s) { _outclusterProducer = s; }

  protected:

    std::string _clusProducer;
    std::string _vtxProducer;
    std::string _outhitProducer;
    std::string _outclusterProducer;

    /// vertex coordinates
    std::vector<double> _vtx_w_cm;
    std::vector<double> _vtx_t_cm;

    // maximum area for a cluster
    double _Amax;
    // maximum distance for a cluster from the vtx
    double _d_max;
    // maximum number of hits for a cluster
    int _max_n_hits;

    double _wire2cm, _time2cm;
    
    void getClusterPoints(const std::vector<unsigned int>& hit_idx_v,
			  larlite::event_hit* evt_hit,
			  std::vector<double>& w_v,
			  std::vector<double>& t_v) const;

    void getClusterBounds(const std::vector<double>& w_v,
			  const std::vector<double>& t_v,
			  const int& pl,
			  std::pair<double,double>& w_bounds,
			  std::pair<double,double>& t_bounds,
			  double& d_max, double& d_min) const;
      
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
