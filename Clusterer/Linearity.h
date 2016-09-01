/**
 * \file Linearity.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class Linearity
 *
 * @author david
 */

/** \addtogroup Clusterer

    @{*/
#ifndef LINEARITY_H
#define LINEARITY_H

#include <iostream>
#include <vector>
#include <cmath>

/**
   \class Linearity
   User defined class Linearity ... these comments are used to generate
   doxygen documentation!
 */
class Linearity{

public:

  /// Default constructor
  Linearity(){}

  Linearity(const std::vector<double>& x_v, const std::vector<double>& y_v);

  /// Default destructor
  ~Linearity(){}

  /// covariance, standard deviation, mean
  double cov (const std::vector<double>& data1,
	      const std::vector<double>& data2);
  
  double stdev(const std::vector<double>& data);
  
  double mean (const std::vector<double>& data);
  
  double linearity(const std::vector<double>& data1,
		   const std::vector<double>& data2);

  double _slope;
  double _covxx;
  double _covyy;
  double _covxy;
  double _lin;
  double _meanx;
  double _meany;
  double _stdx;
  double _stdy;
  int    _dof;

};

#endif
/** @} */ // end of doxygen group 

