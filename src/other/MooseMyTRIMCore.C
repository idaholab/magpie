/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MooseMyTRIMCore.h"
#include "MooseMyTRIMSample.h"

MooseMyTRIMCore::MooseMyTRIMCore(MyTRIM_NS::SimconfType * simconf, MooseMyTRIMSample * sample) :
    MyTRIM_NS::TrimBase(simconf, sample),
    _dim(sample->getDim())
{
}

void
MooseMyTRIMCore::vacancyCreation()
{
  // if MyTRIM was to discard this recoil we'd miss this recoil
  if (_recoil->_state == MyTRIM_NS::IonBase::DELETE)
  {
    // instead we set its energy to 0
    _recoil->_E = 0.0;

    // and set the state to VACANCY below. It well be deleted in the recoil loop instead
  }

  // called if both atoms in the recoil have sufficient energy to leave the site
  _recoil->_state = MyTRIM_NS::IonBase::VACANCY;
}

void
MooseMyTRIMCore::replacementCollision()
{
  // if MyTRIM was to discard this recoil we'd miss this recoil
  if (_recoil->_state == MyTRIM_NS::IonBase::DELETE)
  {
    // instead we set its energy to 0
    _recoil->_E = 0.0;

    // and set the state to REPLACEMENT below. It well be deleted in the recoil loop instead
  }

  // called if only recoil atom has sufficient energy to leave the site
  _recoil->_state = MyTRIM_NS::IonBase::REPLACEMENT;
}
