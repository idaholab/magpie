#include "MooseMyTRIMMaterial.h"

MooseMyTRIMMaterial::MooseMyTRIMMaterial(MyTRIM_NS::SimconfType * simconf, double rho) :
    MyTRIM_NS::MaterialBase(simconf, rho)
{
}

MyTRIM_NS::ElementBase *
MooseMyTRIMMaterial::getElement(unsigned int nn)
{
  // store the recoil element index in the tag field
  tag = nn;
  return element[nn];
}
