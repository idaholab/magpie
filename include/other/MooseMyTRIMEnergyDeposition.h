#ifndef MOOSEMYTRIMENERGYDEPOSITION_H
#define MOOSEMYTRIMENERGYDEPOSITION_H

#include "MooseMyTRIMCore.h"
#include "mytrim/ion.h"

class MooseMyTRIMSample;

/**
 * MyTRIM simulation with energy deposition tracking
 */
class MooseMyTRIMEnergyDeposition : public MooseMyTRIMCore
{
public:
  MooseMyTRIMEnergyDeposition(MyTRIM_NS::SimconfType * simconf, MooseMyTRIMSample * sample,
                              std::vector<std::pair<Point, unsigned int> > & vac,
                              std::vector<std::pair<Point, Real> > & edep);

  virtual void checkPKAState();
  virtual void dissipateRecoilEnergy();
  virtual bool followRecoil();

protected:
  void depositEnergy(MyTRIM_NS::IonBase * ion, Real E);

  /// list of energy deposition points generated during the recoil
  std::vector<std::pair<Point, Real> > & _edep;
};

#endif //MOOSEMYTRIMENERGYDEPOSITION_H
