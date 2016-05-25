#ifndef DiscretePKAPDFBaseBASE_H
#define DiscretePKAPDFBaseBASE_H

#include "MultiIndex.h"

/**
 * Implements a discrete PDF for sampling
 * PKAs
 */
class DiscretePKAPDFBase
{
public:
  DiscretePKAPDFBase(Real magnitude);

  /**
   * A struct storing the inital state of a primary knock-on atom
   */
  struct initialPKAState
  {
    unsigned int _Z;
    Real _mass;
    Real _energy;
    Point _direction;
  };

  /// Uses the discrete probabilities for sampling the initial pka state
  virtual void drawSample(initialPKAState & initial_state) = 0;

protected:
  /// magnitude for correct scaling with potential other DiscretePKAPDFBase objects
  Real _magnitude;
};

#endif // DiscretePKAPDFBase
