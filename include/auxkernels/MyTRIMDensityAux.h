#ifndef MYTRIMDENSITYAUX_H
#define MYTRIMDENSITYAUX_H

#include "AuxKernel.h"
#include "mytrim/simconf.h"

// forward declarations
class MyTRIMRasterizer;
class MyTRIMDensityAux;

template<>
InputParameters validParams<MyTRIMDensityAux>();

class MyTRIMDensityAux : public AuxKernel
{
public:
  MyTRIMDensityAux(const InputParameters & params);
  virtual ~MyTRIMDensityAux() {}

  virtual Real computeValue();

protected:
  const MyTRIMRasterizer & _rasterizer;

  /// number of elements used in the problem
  unsigned int _nvars;

  /// Element masses
  const std::vector<Real> & _trim_mass;

private:
  /// internal TRIM simulation status object
  MyTRIM_NS::SimconfType _simconf;

  /// calculate values only for qp 0 and cache them here
  Real _value_cache;
};

#endif //MYTRIMDENSITYAUX_H
