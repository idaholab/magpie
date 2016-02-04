#ifndef MOOSEMYTRIMMATERIAL_H
#define MOOSEMYTRIMMATERIAL_H

#include "mytrim/material.h"
#include "mytrim/element.h"

/**
 * MyTRIM material class that stores element indices as tags
 */
class MooseMyTRIMMaterial : public MyTRIM_NS::MaterialBase
{
public:
  MooseMyTRIMMaterial(MyTRIM_NS::SimconfType * simconf);

  virtual MyTRIM_NS::ElementBase * getElement(unsigned int nn);

  /**
   * Calculate the density rho by multiplying the fractions t with the atomic weight
   * and dividing by the site volume.
   * This needs to be called before prepare, as prepare renormalizes teh fractions
   * to sum up to one. However we are assuming the fraction difference to one to be
   * the fraction of unoccupied sites.
   */
  virtual void calculateDensity(Real site_volume);
};

#endif //MOOSEMYTRIMMATERIAL_H
