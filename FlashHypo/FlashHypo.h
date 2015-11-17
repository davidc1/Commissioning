/**
 * \file FlashHypo.h
 *
 * \ingroup FlashHypo
 * 
 * \brief Class def header for a class FlashHypo
 *
 * @author Rui
 */

/** \addtogroup FlashHypo

    @{*/
#ifndef FLASHHYPO_H
#define FLASHHYPO_H

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

/**
   \class FlashHypo
   User defined class FlashHypo ... these comments are used to generate
   doxygen documentation!
 */
class FlashHypo{

 public:

  /// Default constructor
  FlashHypo(){}
  
  /// Default destructor
  ~FlashHypo(){}
  
  // Setter function
  bool TrackStart( bool a) { _start =a; return _start;}
  bool TrackEnd  ( bool b) { _end   =b; return _end;}
  double Set_Gap (double x){ _gap   =x; return _gap;}
  
  // Getter function
  std::vector<double> FlashHypothesis(::geoalgo::Trajectory trj) const;
  
  // Calculation fuction
  std::vector<double> PhotonLibrary(::geoalgo::Vector pt_1, ::geoalgo::Vector pt_2, std::vector<double> pe) const;
  
 protected:
  bool   _start = true;
  bool   _end = true;
  double _gap = 0.5;
};

#endif
/** @} */ // end of doxygen group 

