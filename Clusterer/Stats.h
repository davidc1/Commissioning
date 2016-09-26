/**
 * \file Stats.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class Stats
 *
 * @author david
 */

/** \addtogroup Clusterer

    @{*/
#ifndef STATS_H
#define STATS_H

#include <iostream>
#include <cmath>

/**
   \class Stats
   User defined class Stats ... these comments are used to generate
   doxygen documentation!
 */
class Stats{

public:

  /// Default constructor
  Stats(){}

  /// Default destructor
  ~Stats(){}


 private:

  double gamma(const double& Z) const;

  double igf(double S,double Z) const;

  double chisqr(const int& Dof, const double& Cv) const;

};

#endif
/** @} */ // end of doxygen group 

