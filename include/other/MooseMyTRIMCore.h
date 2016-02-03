#ifndef MOOSEMYTRIMCORE_H
#define MOOSEMYTRIMCORE_H

#include "mytrim/trim.h"
#include "MooseTypes.h"

class MooseMyTRIMSample;

/**
 * MyTRIM simulation core class
 */
class MooseMyTRIMCore : public MyTRIM_NS::TrimBase
{
public:
  MooseMyTRIMCore(MyTRIM_NS::SimconfType * simconf, MooseMyTRIMSample * sample, std::vector<std::pair<Point, unsigned int> > & vac);

  virtual void vacancyCreation();
  virtual void checkPKAState();
  virtual void dissipateRecoilEnergy();

protected:
  /// list of vacancies generated during the recoil
  std::vector<std::pair<Point, unsigned int> > & _vac;

  // dimension of the mesh
  const unsigned int _dim;
};

#endif //MOOSEMYTRIMCORE_H
