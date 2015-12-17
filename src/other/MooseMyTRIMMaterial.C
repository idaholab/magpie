#include "MooseMyTRIMMaterial.h"

MooseMyTRIMMaterial::MooseMyTRIMMaterial(MyTRIM_NS::simconfType * simconf_, double rho_) :
    MyTRIM_NS::materialBase(simconf_, rho_)
{
}

MyTRIM_NS::elementBase *
MooseMyTRIMMaterial::getElement(unsigned int nn)
{
  // store the recoil element index in the tag field
  tag = nn;
  return element[nn];
}
