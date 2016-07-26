#ifndef PKAGENERATORNEUTRONICSBASE_H
#define PKAGENERATORNEUTRONICSBASE_H

#include "PKAGeneratorBase.h"
#include "MultiIndex.h"
#include "MyTRIMRasterizer.h"
#include "mytrim/ion.h"
#include "DiscreteFissionPKAPDF.h"

class PKAGeneratorNeutronicsBase;

template<>
InputParameters validParams<PKAGeneratorNeutronicsBase>();

/**
 * A PKAGenerator class that uses neutronics computed reaction rates for
 * sampling the inital state of PKAs.
 */
class PKAGeneratorNeutronicsBase : public PKAGeneratorBase
{
public:
  PKAGeneratorNeutronicsBase(const InputParameters & parameters);

  /// helper function to set pdf based on neutronics information from Mammoth
  virtual void setPDF(Real magnitude, const std::vector<unsigned int> & ZAID, const std::vector<Real> & energies, const MultiIndex<Real> & probabilities) = 0;
};

#endif // PKAGENERATORNEUTRONICSBASE_H
