#ifndef LINEARITY_CXX
#define LINEARITY_CXX

#include "Linearity.h"

Linearity::Linearity(const std::vector<double>& x_v, const std::vector<double>& y_v)
{

  linearity(x_v,y_v);
  
}

// covariance
double Linearity::cov (const std::vector<double>& datax,
		       const std::vector<double>& datay)
{
  if(datax.size() == 0) { std::cout << "zero-vector!" << std::endl; return 0; }
  if(datay.size() == 0) { std::cout << "zero-vector!" << std::endl; return 0; }
  
  double result = 0.0;
  _meanx = mean(datax);
  _meany = mean(datay);
  
  for(size_t i = 0; i < datax.size(); ++i)
    result += (datax[i] - _meanx)*(datay[i] - _meany);
  
  return result/((double)datax.size());
  
}

double Linearity::stdev(const std::vector<double>& data)
{
  if(data.size() == 0) { std::cout << "zero-vector!" << std::endl; return 0; }
  
  double result = 0.0;
  auto    avg   = mean(data);
  for(const auto& d: data)
    result += (d - avg)*(d - avg);
  
  return sqrt(result/((double)data.size()));
}

double Linearity::mean(const std::vector<double>& data)
{
  if(data.size() == 0) { std::cout << "zero-vector!" << std::endl; return 0; }
  
  double result = 0.0;
  
  for(const auto& d : data) 
    result += d;
  
  return (result / ((double)data.size()));
}

double Linearity::linearity(const std::vector<double>& datax,
			    const std::vector<double>& datay)
{
  
  if (datax.size() < 2)
    return 1.;
  
  _covxy = cov(datax,datay);
  _covxx = cov(datax,datax);
  _covyy = cov(datay,datay);
  
  double r_num = _covxy;
  double r_den = sqrt( _covxx * _covyy );
  double r = 0.;
  
  if (r_den == 0)
    r = 0.;
  else
    r = r_num / r_den;
  if (r > 1.) r = 1.;
  if (r < -1) r = -1;
  
  _dof = datax.size() - 2;
  
  _slope = r_num / _covxx;
  
  double lin = sqrt( (1-r*r) * _covyy / _covxx / _dof );
  
  return lin;
  
}

#endif
