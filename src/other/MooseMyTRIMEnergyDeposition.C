/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MooseMyTRIMEnergyDeposition.h"
#include "MooseMyTRIMSample.h"

MooseMyTRIMEnergyDeposition::MooseMyTRIMEnergyDeposition(MyTRIM_NS::SimconfType * simconf, MooseMyTRIMSample * sample,
                                                         std::list<std::pair<Point, unsigned int> > & vac_list,
                                                         std::list<std::pair<Point, Real> > & edep_list) :
    MooseMyTRIMCore(simconf, sample, vac_list),
    _edep_list(edep_list)
{
}

void
MooseMyTRIMEnergyDeposition::checkPKAState()
{
  switch (_pka->_state)
  {
    // PKA is gone, and with it all energy
    case MyTRIM_NS::IonBase::LOST:
      return;

    // only deposit electronic stopping (PKA is moving on)
    case MyTRIM_NS::IonBase::MOVING:
      depositEnergy(_pka, _dee);
      return;

    case MyTRIM_NS::IonBase::INTERSTITIAL:
      // deposit residual energy of the stopped PKA and electronic stopping
      depositEnergy(_pka, _pka->_E + _dee);
      return;

    case MyTRIM_NS::IonBase::REPLACEMENT:
    case MyTRIM_NS::IonBase::SUBSTITUTIONAL:
      // deposit residual energy of the stopped PKA, electronic stopping, and binding energy to the new lattice site
      depositEnergy(_pka, _pka->_E + _element->_Elbind + _dee);
      return ;
  }
}

void
MooseMyTRIMEnergyDeposition::dissipateRecoilEnergy()
{
  // new recoil is not leaving its lattice site, reimburse binding energy
  depositEnergy(_recoil, _recoil->_E + _element->_Elbind);
}

void
MooseMyTRIMEnergyDeposition::vacancyCreation()
{
  // TODO: lattice relaxation around the vacancy?
  MooseMyTRIMCore::vacancyCreation();
}

bool
MooseMyTRIMEnergyDeposition::followRecoil()
{
  // TODO: if we ever return false here we need to deposit the discarded recoil energy
  return true;
}

void
MooseMyTRIMEnergyDeposition::depositEnergy(MyTRIM_NS::IonBase * ion, Real E)
{
  _edep_list.push_back(std::make_pair(Point(ion->_pos(0), ion->_pos(1), _dim == 2 ? 0.0 :  ion->_pos(2)), E));
}
