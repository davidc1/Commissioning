#ifndef STATS_CXX
#define STATS_CXX

#include "Stats.h"

double Stats::gamma(const double& Z) const
{

  const double RECIP_E = 0.36787944117144232159552377016147;  // RECIP_E = (E^-1) = (1.0 / E)
  const double TWOPI = 6.283185307179586476925286766559;  // TWOPI = 2.0 * PI

  double D = 1.0 / (10.0 * Z);
  D = 1.0 / ((12 * Z) - D);
  D = (D + Z) * RECIP_E;
  D = pow(D, Z);
  D *= sqrt(TWOPI / Z);

  return D;
}

double Stats::igf(double S, double Z) const
{
  if(Z < 0.0)
    {
      return 0.0;
    }
  double Sc = (1.0 / S);
  Sc *= pow(Z, S);
  Sc *= exp(-Z);

  double Sum = 1.0;
  double Nom = 1.0;
  double Denom = 1.0;

  for(int I = 0; I < 200; I++)
    {
      Nom *= Z;
      S++;
      Denom *= S;
      Sum += (Nom / Denom);
    }

  return Sum * Sc;
}

double Stats::chisqr(const int& Dof, const double& Cv) const
{
  if(Cv < 0 || Dof < 1)
    {
      return 0.0;
    }
  double K = ((double)Dof) * 0.5;
  double X = Cv * 0.5;
  if(Dof == 2)
    {
      return exp(-1.0 * X);
    }

  double PValue = igf(K, X);
  if(isnan(PValue) || isinf(PValue) || PValue <= 1e-8)
    {
      return 1e-14;
    }

  PValue /= gamma(K);

  return (1.0 - PValue);
}

#endif
