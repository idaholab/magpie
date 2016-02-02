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
  MooseMyTRIMMaterial(MyTRIM_NS::SimconfType * simconf, double rho);

  virtual MyTRIM_NS::ElementBase * getElement(unsigned int nn);
};

#endif //MOOSEMYTRIMMATERIAL_H
