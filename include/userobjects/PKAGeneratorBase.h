#ifndef PKAGENERATORBASE_H
#define PKAGENERATORBASE_H

#include "DiscreteElementUserObject.h"
#include "MyTRIMRasterizer.h"
#include "mytrim/ion.h"

class PKAGeneratorBase;

template<>
InputParameters validParams<PKAGeneratorBase>();

/**
 * Abstract base class for PKA calculation UOs that plug into MyTRIMRasterizer
 * to generate a set of PKAs for the current element
 */
class PKAGeneratorBase : public DiscreteElementUserObject
{
public:
  PKAGeneratorBase(const InputParameters & parameters);

  /**
   * Append the ions for the current element and time window dt.
   * The element volume is passed in as it is computed in the MyTRIMRasterizer anyways.
   */
  virtual void appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list, Real dt, Real vol, const MyTRIMRasterizer::AveragedData &) const = 0;

  virtual void initialize() {}

protected:
  /// helper function to set the ion position to a random location in the current element
  void setPosition(MyTRIM_NS::IonBase & ion) const;

  /// helper function to set the ion direction to a random direction
  void setRandomDirection(MyTRIM_NS::IonBase & ion) const;

  /// Return a point with random uniformly distributed coordinates in the unit cube (temp variables are required to ensure execution order!)
  Point getRandomPoint() const { const Real X = getRandomReal(), Y = getRandomReal(), Z = getRandomReal(); return Point(X,Y,Z); }
};

#endif // PKAGENERATORBASE_H
