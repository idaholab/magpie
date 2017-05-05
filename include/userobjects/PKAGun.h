#ifndef PKAGUN_H
#define PKAGUN_H

#include "PKAFixedPointGenerator.h"

class PKAGun;

template<>
InputParameters validParams<PKAGun>();

/**
 * Starts PKAs at a fixed point in a fixed direction
 */
class PKAGun : public PKAFixedPointGenerator
{
public:
  PKAGun(const InputParameters & parameters);

protected:
  /// provides a mean to override the angular distribution of the PKAs in derived class
  virtual void setDirection(MyTRIM_NS::IonBase & ion) const;

  /// the direction along which the PKAs move
  const Point _direction;
};

#endif // PKAGUN_H
