/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "MooseMyTRIMCore.h"
#include "mytrim/ion.h"

class MooseMyTRIMSample;

/**
 * MyTRIM simulation with energy deposition tracking
 */
class MooseMyTRIMEnergyDeposition : public MooseMyTRIMCore
{
public:
  MooseMyTRIMEnergyDeposition(MyTRIM_NS::SimconfType * simconf,
                              MooseMyTRIMSample * sample,
                              std::list<std::pair<Point, Real>> & edep_list);

  virtual void checkPKAState();
  virtual void dissipateRecoilEnergy();
  virtual bool followRecoil();

protected:
  void depositEnergy(MyTRIM_NS::IonBase * ion, Real E);

  /// list of energy deposition points generated during the recoil
  std::list<std::pair<Point, Real>> & _edep_list;
};
