#include "MooseMyTRIMMaterial.h"

MooseMyTRIMMaterial::MooseMyTRIMMaterial(double rho) :
    MyTRIM_NS::materialBase(rho)
{
}

MyTRIM_NS::elementBase *
MooseMyTRIMMaterial::getElement(unsigned int nn)
{
  // store the recoil element index in the tag field
  tag = nn;
  return element[nn];
}
