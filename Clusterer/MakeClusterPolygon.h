/**
 * \file MakeClusterPolygon.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class MakeClusterPolygon
 *
 * @author david caratelli
 */

/** \addtogroup Clusterer

    @{*/

#ifndef MAKECLUSTERPOLYGON_H
#define MAKECLUSTERPOLYGON_H

#include "Analysis/ana_base.h"
#include "RecoTool/ClusterRecoUtil/Base/Polygon2D.h"

namespace larlite {
  /**
     \class MakeClusterPolygon
     User custom analysis class made by david caratelli
   */
  class MakeClusterPolygon {
  
  public:

    /// Default constructor
    MakeClusterPolygon(){}

    /// Default destructor
    virtual ~MakeClusterPolygon(){}

    /**
     * @brief return a Polygon2D object given a list of larlite hits
     */
    Polygon2D MakePolygon(const std::vector<larlite::hit>& hits,
			  const double& frac=1.);

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
