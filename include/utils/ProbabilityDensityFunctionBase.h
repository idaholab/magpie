#ifndef PROBABILITYDENSITYFUNCTIONBASE_H
#define PROBABILITYDENSITYFUNCTIONBASE_H

#include "Moose.h"
#include "RandomInterface.h"

template<class T>
class ProbabilityDensityFunctionBase
{
public:
  ProbabilityDensityFunctionBase(Real magnitude = 1.0) :
    _magnitude(magnitude)
  {
  }

  /**
   * Override this function for implementing the
   * desired sampling behavior
   */
  virtual T drawSample() = 0;

protected:
  /// Magnitude to allow scaling with other PDFs
  Real _magnitude;
};

#endif //PROBABILITYDENSITYFUNCTIONBASE_H
